#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

export SUNPY_CONFIGDIR="${SUNPY_CONFIGDIR:-/tmp/sunpy}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-codex}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp}"

PYTHON_BIN="${PYTHON_BIN:-python}"

log() {
  printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Required command not found: $1" >&2
    exit 1
  }
}

prompt_yes_no() {
  local prompt="$1"
  local default="${2:-N}"
  local answer
  if [[ "$default" == "Y" ]]; then
    read -r -p "$prompt [Y/n] " answer
    answer="${answer:-Y}"
  else
    read -r -p "$prompt [y/N] " answer
    answer="${answer:-N}"
  fi
  [[ "$answer" =~ ^[Yy]$ ]]
}

prompt_choice() {
  local prompt="$1"
  shift
  local choices=("$@")
  local answer
  echo "$prompt"
  for choice in "${choices[@]}"; do
    echo "  $choice"
  done
  read -r -p "Select an option: " answer
  echo "$answer"
}

mtime_or_zero() {
  local path="$1"
  if [[ -e "$path" ]]; then
    stat -f "%m" "$path"
  else
    echo 0
  fi
}

check_required_files() {
  log "Checking input file integrity"
  "$PYTHON_BIN" - <<'PY'
from pathlib import Path
from astropy.io import fits
from alignment_common import DATA_DIR, HMI_PATH, IRIS_SOURCE_PATH, SST_WB_SOURCE_PATH, SST_NB_SOURCE_PATH, file_is_truncated_primary_hdu

required = [
    ("HMI", HMI_PATH),
    ("IRIS", IRIS_SOURCE_PATH),
    ("SST WB", SST_WB_SOURCE_PATH),
]

for label, path in required:
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"Missing required file: {label}: {path}")
    fits.getheader(path, 0)
    print(f"OK: {label}: {path}")

nb_path = Path(SST_NB_SOURCE_PATH)
if nb_path.exists():
    fits.getheader(nb_path, 0)
    truncated = file_is_truncated_primary_hdu(nb_path)
    print(f"NB source: {nb_path} | truncated={truncated}")
else:
    print(f"NB source missing: {nb_path}")
PY
}

clear_old_outputs() {
  log "Removing old output products"
  mkdir -p "$REPO_DIR/outputs"
  rm -f "$REPO_DIR"/outputs/*.png
  rm -f "$REPO_DIR"/outputs/*.gif
  rm -f "$REPO_DIR"/outputs/*.json
  rm -f "$REPO_DIR"/outputs/*_aligned.fits
  rm -f "$REPO_DIR"/outputs/alignment_check_sunpy_frame_*.png
}

run_initial_alignment() {
  log "Running initial coarse alignment"
  "$PYTHON_BIN" "$REPO_DIR/align_iris_sst.py" --reset
}

verify_iris_save() {
  local report_mtime_before="$1"
  local fits_mtime_before="$2"
  local report_path="$REPO_DIR/outputs/alignment_report.json"
  local iris_path="$REPO_DIR/outputs/iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits"
  local report_mtime_after fits_mtime_after
  report_mtime_after="$(mtime_or_zero "$report_path")"
  fits_mtime_after="$(mtime_or_zero "$iris_path")"

  log "Checking saved IRIS alignment values"
  "$PYTHON_BIN" - <<'PY'
import json
from astropy.io import fits
from pathlib import Path

report_path = Path("outputs/alignment_report.json")
iris_path = Path("outputs/iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits")

report = json.loads(report_path.read_text())
hdr = fits.getheader(iris_path, 0)
print("IRIS report shift:", report["iris"]["world_shift_arcsec"])
print("IRIS manual shift:", report.get("iris_manual_adjustment_arcsec", {}))
print("IRIS CRVAL1/CRVAL2:", hdr.get("CRVAL1"), hdr.get("CRVAL2"))
print("IRIS ALNIRX/ALNIRY:", hdr.get("ALNIRX"), hdr.get("ALNIRY"))
PY

  if [[ "$report_mtime_after" -le "$report_mtime_before" && "$fits_mtime_after" -le "$fits_mtime_before" ]]; then
    echo "No new IRIS save was detected." >&2
    prompt_yes_no "Continue anyway?" || exit 1
  fi
}

verify_sst_save() {
  local report_mtime_before="$1"
  local fits_mtime_before="$2"
  local manual_mtime_before="$3"
  local report_path="$REPO_DIR/outputs/alignment_report.json"
  local sst_path="$REPO_DIR/outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits"
  local manual_path="$REPO_DIR/outputs/sst_manual_adjustment.json"
  local report_mtime_after fits_mtime_after manual_mtime_after
  report_mtime_after="$(mtime_or_zero "$report_path")"
  fits_mtime_after="$(mtime_or_zero "$sst_path")"
  manual_mtime_after="$(mtime_or_zero "$manual_path")"

  log "Checking saved SST alignment values"
  "$PYTHON_BIN" - <<'PY'
import json
from astropy.io import fits
from pathlib import Path

report_path = Path("outputs/alignment_report.json")
sst_path = Path("outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits")
manual_path = Path("outputs/sst_manual_adjustment.json")

report = json.loads(report_path.read_text())
hdr = fits.getheader(sst_path, 0)
print("SST report total shift:", report.get("current_total_sst_shift_arcsec", {}))
print("SST manual shift:", report.get("manual_adjustment_arcsec", {}))
print("SST header ALNWX/ALNWY:", hdr.get("ALNWX"), hdr.get("ALNWY"))
print("SST header ALNMANX/ALNMANY:", hdr.get("ALNMANX"), hdr.get("ALNMANY"))
if manual_path.exists():
    print("SST manual save json:", json.loads(manual_path.read_text()))
PY

  if [[ "$report_mtime_after" -le "$report_mtime_before" && "$fits_mtime_after" -le "$fits_mtime_before" && "$manual_mtime_after" -le "$manual_mtime_before" ]]; then
    echo "No new SST save was detected." >&2
    prompt_yes_no "Continue anyway?" || exit 1
  fi
}

run_iris_app() {
  local report_path="$REPO_DIR/outputs/alignment_report.json"
  local iris_path="$REPO_DIR/outputs/iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits"
  local report_mtime_before fits_mtime_before
  report_mtime_before="$(mtime_or_zero "$report_path")"
  fits_mtime_before="$(mtime_or_zero "$iris_path")"

  log "Launching IRIS manual alignment app"
  echo "Use the app, press Save, then Quit And Close."
  "$PYTHON_BIN" "$REPO_DIR/interactive_iris_hmi_align.py"

  verify_iris_save "$report_mtime_before" "$fits_mtime_before"
}

run_sst_app() {
  local report_path="$REPO_DIR/outputs/alignment_report.json"
  local sst_path="$REPO_DIR/outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits"
  local manual_path="$REPO_DIR/outputs/sst_manual_adjustment.json"
  local report_mtime_before fits_mtime_before manual_mtime_before
  report_mtime_before="$(mtime_or_zero "$report_path")"
  fits_mtime_before="$(mtime_or_zero "$sst_path")"
  manual_mtime_before="$(mtime_or_zero "$manual_path")"

  log "Launching SST manual alignment app"
  echo "Use the app, press Save, then Quit And Close."
  "$PYTHON_BIN" "$REPO_DIR/interactive_sst_manual_align.py"

  verify_sst_save "$report_mtime_before" "$fits_mtime_before" "$manual_mtime_before"
}

regenerate_wb_products() {
  log "Regenerating WB comparison products"
  "$PYTHON_BIN" "$REPO_DIR/plot_vs_hmi.py" \
    --sst-kind wb \
    --num-samples 3 \
    --output "$REPO_DIR/outputs/compare_to_hmi_wb.png" \
    --match-output "$REPO_DIR/outputs/compare_to_hmi_wb_matches.json"

  "$PYTHON_BIN" "$REPO_DIR/plot_sst_iris_time_overlap.py" \
    --sst-kind wb \
    --num-samples 3 \
    --output "$REPO_DIR/outputs/sst_iris_time_overlap_wb.png" \
    --match-output "$REPO_DIR/outputs/sst_iris_time_overlap_matches_wb.json"

  "$PYTHON_BIN" "$REPO_DIR/sunpy_alignment_check.py"

  echo "Generated WB products:"
  echo "  outputs/compare_to_hmi_wb.png"
  echo "  outputs/sst_iris_time_overlap_wb.png"
  echo "  outputs/alignment_check_sunpy.gif"
}

review_wb_loop() {
  while true; do
    echo
    echo "Please inspect the WB comparison products now."
    local choice
    choice="$(prompt_choice "What do you want to do next?" \
      "1) Continue to the remaining products" \
      "2) Rerun IRIS interactive only" \
      "3) Rerun SST interactive only" \
      "4) Rerun both interactive steps" \
      "5) Regenerate WB plots only" \
      "6) Exit")"

    case "$choice" in
      1)
        return 0
        ;;
      2)
        run_iris_app
        regenerate_wb_products
        ;;
      3)
        run_sst_app
        regenerate_wb_products
        ;;
      4)
        run_iris_app
        run_sst_app
        regenerate_wb_products
        ;;
      5)
        regenerate_wb_products
        ;;
      6)
        echo "Stopping after WB products."
        exit 0
        ;;
      *)
        echo "Unrecognized choice: $choice"
        ;;
    esac
  done
}

regenerate_nb_products() {
  log "Regenerating NB comparison plot"
  "$PYTHON_BIN" "$REPO_DIR/plot_vs_hmi.py" \
    --sst-kind nb \
    --num-samples 3 \
    --output "$REPO_DIR/outputs/compare_to_hmi_nb.png" \
    --match-output "$REPO_DIR/outputs/compare_to_hmi_nb_matches.json"

  "$PYTHON_BIN" "$REPO_DIR/plot_sst_iris_time_overlap.py" \
    --sst-kind nb \
    --num-samples 3 \
    --output "$REPO_DIR/outputs/sst_iris_time_overlap_nb.png" \
    --match-output "$REPO_DIR/outputs/sst_iris_time_overlap_matches_nb.json"
}

build_nb_cube() {
  log "Building aligned NB cube"
  "$PYTHON_BIN" "$REPO_DIR/update_sst_nb_wcs.py"
}

main() {
  require_command "$PYTHON_BIN"
  check_required_files
  clear_old_outputs
  run_initial_alignment

  run_iris_app
  run_sst_app

  regenerate_wb_products
  review_wb_loop

  regenerate_nb_products

  if prompt_yes_no "Build the aligned NB cube as well? This can duplicate ~78 GB temporarily and may hit storage limits?" "N"; then
    build_nb_cube
  else
    echo "Skipping aligned NB cube rewrite."
  fi

  log "Pipeline complete"
}

main "$@"
