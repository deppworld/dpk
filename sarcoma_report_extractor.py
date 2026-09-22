"""
sarcoma_report_extractor.py
===========================
Extract Epignostix SARCOMA Methylation Classifier HTML reports into an Excel
summary for LabVantage ingestion.

Differences from the CNS extractor
----------------------------------
  * No Family / Superfamily levels          -> removed
  * No MGMT prediction block                -> removed
  * Therefore NO OCR: pytesseract, opencv-python and numpy are NOT required.
    Tesseract-OCR does not need to be installed on the scheduler machine.
  * The sarcoma report presents results as a proper matrix table
    (Level | Prediction | Calibrated Score | Interpretation | Evidence),
    so parsing is by column header rather than by label-in-first-cell.

Dependencies
------------
    pip install pandas beautifulsoup4 lxml openpyxl

Usage
-----
    # 1. Testing: point at the folder that holds the HTML files
    python sarcoma_report_extractor.py --folder "C:\\Temp\\sarcoma_test"

    # 2. Production (Task Scheduler): edit WORKING_FOLDER below, then just
    python sarcoma_report_extractor.py

    # Other options
    --mode auto|flat|runs   how to look for reports      (default: auto)
    --out  <path>           override the output .xlsx location
    --reprocess             ignore the processed log and redo everything
    --quiet                 suppress console output (Task Scheduler)

Author: Deepak  <depp@live.in>
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from datetime import datetime

import pandas as pd
from bs4 import BeautifulSoup

# ===========================================================================
#  CONFIG  —  edit these two for production
# ===========================================================================

# Testing now:  a local folder containing the .html reports.
# Later:        the network path, e.g.
#   r"\\hlm\data\apps\ClinicalLabs\ADL Data Storage\iScan_Archive\Methylation_Clinical\Output-Sarcoma"
WORKING_FOLDER = r"C:\Temp\sarcoma_test"

# "auto" = decide automatically (recommended)
#   flat : HTML files sit directly in WORKING_FOLDER (or one level down)
#   runs : WORKING_FOLDER\<run_folder>\report\*.html   (the CNS layout)
FOLDER_MODE = "auto"

# Name of the Excel file written into each folder that contains reports.
OUTPUT_FILENAME = "Sarcoma_Report_Summary.xlsx"

# Log of what has already been done, so the scheduler does not redo work.
# Kept beside WORKING_FOLDER by default; move it to a local drive if the
# network share is read-only for the service account.
PROCESSED_LOG_NAME = "processed_sarcoma.txt"
PROCESSED_LOG_DIR = None          # None -> use WORKING_FOLDER

# Write scores as numbers (True) or as text (False).
# Numbers sort and filter correctly in Excel; switch to False only if the
# LabVantage import mapping insists on strings.
SCORES_AS_NUMBERS = True

# Lab interpretation thresholds. NOTE: Epignostix itself reports a simple
# binary (match >= 0.9, no match < 0.9) and that value is captured separately
# in the *_Interpretation columns, so the two never have to be reconciled.
THRESHOLD_MATCHED = 0.90
THRESHOLD_SUGGESTIVE = 0.50

# Column order written to Excel.
COLUMNS = [
    "LabVantage", "MRN", "Sample_Identifier", "Sentrix_ID",
    "Array_Type", "Material_Type", "Biological_Sex",
    "Class", "Class_Score", "Class_Report", "Class_Interpretation", "Class_Evidence",
    "Subclass", "Subclass_Score", "Subclass_Report", "Subclass_Interpretation", "Subclass_Evidence",
    "Classifier_Version", "Report_Generated",
    "Description", "File_Name", "Source_Folder", "Parsed_At",
]

log = logging.getLogger("sarcoma")


# ===========================================================================
#  HELPERS
# ===========================================================================

def get_report(score) -> str:
    """Lab three-tier interpretation of a calibrated score."""
    try:
        score = float(score)
    except (TypeError, ValueError):
        return ""
    if score >= THRESHOLD_MATCHED:
        return "Matched"
    if THRESHOLD_SUGGESTIVE <= score < THRESHOLD_MATCHED:
        return "Suggestive: correlate clinically"
    return "No-match"


def to_score(text):
    """Parse a calibrated score cell. Returns float, str or '' per config."""
    if text is None:
        return ""
    cleaned = str(text).strip().replace(",", ".")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", cleaned)
    if not m:
        return ""
    try:
        value = float(m.group(0))
    except ValueError:
        return ""
    return value if SCORES_AS_NUMBERS else f"{value}"


def clean_interpretation(text: str) -> str:
    """'no match | ?'  ->  'no match'   (drops the symbol glyph)."""
    if not text:
        return ""
    return str(text).split("|")[0].strip()


def setup_logging(log_path: str | None, quiet: bool) -> None:
    log.handlers.clear()
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s")

    if not quiet:
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(fmt)
        log.addHandler(console)

    if log_path:
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            fh = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            fh.setFormatter(fmt)
            log.addHandler(fh)
        except OSError as exc:                       # read-only share, etc.
            log.warning("Could not open log file %s (%s) - console only", log_path, exc)


# ===========================================================================
#  PROCESSED LOG
# ===========================================================================

def processed_log_path(root: str) -> str:
    folder = PROCESSED_LOG_DIR or root
    return os.path.join(folder, PROCESSED_LOG_NAME)


def load_processed(root: str) -> set[str]:
    path = processed_log_path(root)
    if not os.path.exists(path):
        return set()
    done = set()
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                key = line.strip().split("\t")[0]
                if key:
                    done.add(key)
    except OSError as exc:
        log.warning("Could not read processed log (%s) - treating all as new", exc)
    return done


def save_processed(root: str, key: str) -> None:
    path = processed_log_path(root)
    try:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{key}\t{datetime.now():%Y-%m-%d %H:%M:%S}\n")
    except OSError as exc:
        log.warning("Could not write processed log (%s)", exc)


# ===========================================================================
#  DISCOVERY
# ===========================================================================

def _html_in(folder: str) -> list[str]:
    try:
        return sorted(
            os.path.join(folder, f)
            for f in os.listdir(folder)
            if f.lower().endswith((".html", ".htm"))
        )
    except OSError:
        return []


def discover(root: str, mode: str) -> dict[str, list[str]]:
    """Return {folder_containing_reports: [html files]}."""
    if not os.path.isdir(root):
        log.error("Working folder does not exist: %s", root)
        return {}

    groups: dict[str, list[str]] = {}

    # --- runs layout:  root\<run>\report\*.html -----------------------------
    if mode in ("auto", "runs"):
        for name in sorted(os.listdir(root)):
            run = os.path.join(root, name)
            if not os.path.isdir(run):
                continue
            report_dir = os.path.join(run, "report")
            if os.path.isdir(report_dir):
                files = _html_in(report_dir)
                if files:
                    groups[report_dir] = files
        if mode == "runs":
            return groups

    # --- flat layout: HTML directly in root, or one level down --------------
    if mode in ("auto", "flat"):
        files = _html_in(root)
        if files:
            groups.setdefault(root, files)
        for name in sorted(os.listdir(root)):
            sub = os.path.join(root, name)
            if not os.path.isdir(sub) or sub in groups:
                continue
            if os.path.basename(sub).lower() == "report":
                continue
            files = _html_in(sub)
            if files:
                groups.setdefault(sub, files)

    return groups


# ===========================================================================
#  PARSING
# ===========================================================================

def _blank_row() -> dict:
    return {c: "" for c in COLUMNS}


def parse_sample_table(soup, row: dict) -> None:
    """Sample Identifier | Sentrix ID | Array type | Material type | Sex."""
    for table in soup.find_all("table"):
        headers = [th.get_text(" ", strip=True) for th in table.find_all("th")]
        if not any("sample identifier" in h.lower() for h in headers):
            continue

        body = table.find("tbody") or table
        for tr in body.find_all("tr"):
            values = [td.get_text(" ", strip=True) for td in tr.find_all("td")]
            if not values:
                continue

            pairs = dict(zip([h.lower() for h in headers], values))
            sample = pairs.get("sample identifier", values[0] if values else "")

            row["Sample_Identifier"] = sample
            row["Sentrix_ID"] = pairs.get("sentrix id", "")
            row["Array_Type"] = pairs.get("array type", "")
            row["Material_Type"] = pairs.get("material type", "")
            row["Biological_Sex"] = pairs.get("biological sex", "")

            # LabVantage-MRN split. If there is no hyphen the whole identifier
            # goes into LabVantage rather than being silently dropped.
            if sample and "-" in sample:
                lv, mrn = sample.rsplit("-", 1)
                row["LabVantage"] = lv.strip()
                row["MRN"] = mrn.strip()
            else:
                row["LabVantage"] = sample
                row["MRN"] = ""
            return


def parse_classifier_table(soup, row: dict) -> None:
    """Hierarchical results: Level | Prediction | Calibrated Score | ... ."""
    for table in soup.find_all("table"):
        headers = [th.get_text(" ", strip=True).lower() for th in table.find_all("th")]
        if not ("level" in headers and "prediction" in headers):
            continue

        idx = {name: i for i, name in enumerate(headers)}
        i_level = idx.get("level", 0)
        i_pred = idx.get("prediction", 1)
        i_score = idx.get("calibrated score", 2)
        i_interp = idx.get("interpretation")
        i_evid = idx.get("evidence")

        body = table.find("tbody") or table
        for tr in body.find_all("tr"):
            cells = [td.get_text(" ", strip=True) for td in tr.find_all("td")]
            if len(cells) <= max(i_level, i_pred):
                continue

            level = cells[i_level].strip().lower()
            pick = lambda i: cells[i] if (i is not None and i < len(cells)) else ""

            # Check subclass FIRST — "subclass" also contains "class".
            # Family / Superfamily rows, if a future version emits them,
            # fall through and are ignored by design.
            if "subclass" in level:
                prefix = "Subclass"
            elif "class" in level:
                prefix = "Class"
            else:
                continue

            row[prefix] = pick(i_pred)
            row[f"{prefix}_Score"] = to_score(pick(i_score))
            row[f"{prefix}_Interpretation"] = clean_interpretation(pick(i_interp))
            row[f"{prefix}_Evidence"] = pick(i_evid)
        return


def parse_description(lines: list[str], row: dict) -> None:
    captured, capturing = [], False
    stop = ("Copy Number Variation Profile", "Quality Control",
            "MGMT Prediction", "Disclaimer", "Run information")
    for line in lines:
        if line == "Description" or "Methylation Subclass Description" in line:
            capturing = True
            continue
        if capturing and any(line.startswith(s) for s in stop):
            break
        if capturing:
            captured.append(line)
    row["Description"] = " ".join(captured).strip()


def parse_run_info(text: str, row: dict) -> None:
    m = re.search(r"Classifier\s+version\s*(v[\d.]+)", text, re.I) \
        or re.search(r"Classifier\s+(v[\d.]+)", text, re.I)
    if m:
        row["Classifier_Version"] = m.group(1)

    m = re.search(r"Report generated on:\s*(.+)", text)
    if m:
        row["Report_Generated"] = m.group(1).strip()


def parse_html(file_path: str) -> dict:
    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        soup = BeautifulSoup(fh, "html.parser")

    for tag in soup(["script", "style"]):
        tag.decompose()

    text = soup.get_text("\n")
    lines = [l.strip() for l in text.splitlines() if l.strip()]

    row = _blank_row()
    row["File_Name"] = os.path.basename(file_path)
    row["Source_Folder"] = os.path.dirname(file_path)
    row["Parsed_At"] = f"{datetime.now():%Y-%m-%d %H:%M:%S}"

    parse_sample_table(soup, row)
    parse_classifier_table(soup, row)
    parse_description(lines, row)
    parse_run_info(text, row)

    row["Class_Report"] = get_report(row["Class_Score"])
    row["Subclass_Report"] = get_report(row["Subclass_Score"])

    if not row["Class"]:
        log.warning("No Class prediction found in %s", row["File_Name"])

    return row


# ===========================================================================
#  EXCEL OUTPUT
# ===========================================================================

def write_excel(rows: list[dict], out_file: str) -> str | None:
    df = pd.DataFrame(rows, columns=COLUMNS)
    try:
        df.to_excel(out_file, index=False, sheet_name="Sarcoma")
        return out_file
    except PermissionError:
        # File open in Excel, or locked by another scheduler run.
        alt = out_file.replace(".xlsx", f"_{datetime.now():%Y%m%d_%H%M%S}.xlsx")
        log.warning("%s is locked - writing %s instead", out_file, os.path.basename(alt))
        df.to_excel(alt, index=False, sheet_name="Sarcoma")
        return alt
    except OSError as exc:
        log.error("Could not write %s: %s", out_file, exc)
        return None


# ===========================================================================
#  MAIN
# ===========================================================================

def process_group(folder: str, files: list[str], out_override: str | None) -> int:
    rows = []
    for path in files:
        log.info("Parsing %s", os.path.basename(path))
        try:
            rows.append(parse_html(path))
        except Exception as exc:                          # noqa: BLE001
            log.exception("Failed to parse %s: %s", path, exc)

    if not rows:
        return 0

    out_file = out_override or os.path.join(folder, OUTPUT_FILENAME)
    written = write_excel(rows, out_file)
    if written:
        log.info("Wrote %d row(s) -> %s", len(rows), written)
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="Epignostix SARCOMA report -> Excel")
    ap.add_argument("--folder", default=WORKING_FOLDER,
                    help="folder to scan (local path or UNC network path)")
    ap.add_argument("--mode", default=FOLDER_MODE, choices=["auto", "flat", "runs"])
    ap.add_argument("--out", default=None, help="explicit output .xlsx path")
    ap.add_argument("--reprocess", action="store_true",
                    help="ignore the processed log and redo everything")
    ap.add_argument("--quiet", action="store_true", help="no console output")
    args = ap.parse_args()

    root = os.path.abspath(args.folder) if not args.folder.startswith("\\\\") else args.folder

    setup_logging(os.path.join(root, "sarcoma_extractor.log"), args.quiet)
    log.info("=" * 66)
    log.info("Sarcoma extractor starting | folder=%s | mode=%s", root, args.mode)

    groups = discover(root, args.mode)
    if not groups:
        log.warning("No HTML reports found under %s", root)
        return 1

    processed = set() if args.reprocess else load_processed(root)
    total_rows = 0
    total_groups = 0

    for folder, files in groups.items():
        key = os.path.relpath(folder, root) if folder != root else "."
        if key in processed:
            log.info("Skipping %s (already processed)", key)
            continue

        log.info("Processing %s (%d file(s))", key, len(files))
        try:
            n = process_group(folder, files, args.out)
        except Exception as exc:                          # noqa: BLE001
            log.exception("Error in %s: %s", folder, exc)
            continue

        if n:
            total_rows += n
            total_groups += 1
            if not args.reprocess:
                save_processed(root, key)

    log.info("Done. %d folder(s), %d row(s).", total_groups, total_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
