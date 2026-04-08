from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/sunpy")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from alignment_common import (
    OUTPUT_DIR,
    REPORT_PATH,
    SST_NB_SOURCE_PATH,
    load_hmi_map,
    make_iris_map,
    make_nb_wing_map,
    make_wb_map,
    normalize_for_display,
    project_map_to_extent,
    read_iris_times,
    read_wb_corners,
    read_wb_times,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare the current IRIS solution and the current SST WB/NB solution separately against HMI."
    )
    parser.add_argument("--sst-kind", choices=["wb", "nb"], default="wb", help="Use WB or NB wing images for the SST column.")
    parser.add_argument("--nb-path", default=str(SST_NB_SOURCE_PATH), help="NB cube to use when --sst-kind nb.")
    parser.add_argument("--scan-indices", default="", help="Comma-separated SST scan indices. If omitted, sample automatically.")
    parser.add_argument("--num-samples", type=int, default=3, help="Number of SST scans to sample when scan indices are omitted.")
    parser.add_argument(
        "--output",
        default=str(OUTPUT_DIR / "compare_to_hmi.png"),
        help="Output PNG path.",
    )
    parser.add_argument(
        "--match-output",
        default=str(OUTPUT_DIR / "compare_to_hmi_matches.json"),
        help="Output JSON path with the SST/IRIS frame pairing used in the figure.",
    )
    return parser.parse_args()


def choose_scan_indices(report: dict, num_samples: int) -> list[int]:
    wb_times = read_wb_times(Path(report["sst_wb_output"]))
    if num_samples <= 1:
        return [0]
    return np.linspace(0, len(wb_times) - 1, num_samples, dtype=int).tolist()


def nearest_iris_index(iris_times, target_time):
    offsets = np.abs((iris_times - target_time).sec)
    idx = int(np.argmin(offsets))
    return idx, float(offsets[idx])


def build_sst_map(hmi_map, report: dict, sst_kind: str, scan_index: int, nb_path: Path):
    wb_path = Path(report["sst_wb_output"])
    if sst_kind == "wb":
        sst_map, (corners, _) = make_wb_map(hmi_map, scan_index, wb_path=wb_path)
        label = "WB"
    else:
        sst_map, wing_indices = make_nb_wing_map(hmi_map, scan_index, nb_path=nb_path, wb_path=wb_path)
        corners = read_wb_corners(wb_path)[scan_index]
        label = f"NB wings {wing_indices[0]} / {wing_indices[1]}"
    return sst_map, corners, label


def extent_from_footprint(corners: np.ndarray, margin_arcsec: float = 80.0) -> list[float]:
    x0 = float(corners[:, 0, 0].mean())
    x1 = float(corners[:, 1, 0].mean())
    y0 = float(corners[0, :, 1].mean())
    y1 = float(corners[1, :, 1].mean())
    cx = 0.5 * (x0 + x1)
    cy = 0.5 * (y0 + y1)
    half = 0.5 * max(abs(x1 - x0), abs(y1 - y0)) + margin_arcsec
    return [
        cx - half,
        cx + half,
        cy - half,
        cy + half,
    ]


def plot_overlay(
    ax,
    hmi_map,
    overlay_map,
    extent: list[float],
    title: str,
    shape: tuple[int, int] = (800, 800),
    overlay_alpha: float = 0.55,
    hmi_alpha: float = 0.65,
) -> None:
    hmi_image = normalize_for_display(project_map_to_extent(hmi_map, extent, shape))
    overlay_image = normalize_for_display(project_map_to_extent(overlay_map, extent, shape))
    ax.imshow(hmi_image, extent=extent, origin="lower", cmap="binary_r", alpha=hmi_alpha, aspect="equal")
    ax.imshow(overlay_image, extent=extent, origin="lower", cmap="gray", alpha=overlay_alpha, aspect="equal")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("Helioprojective Longitude (Solar-X)")
    ax.set_ylabel("Helioprojective Latitude (Solar-Y)")
    ax.set_title(title)


def main() -> None:
    args = parse_args()
    with open(REPORT_PATH) as fh:
        report = json.load(fh)
    hmi_map = load_hmi_map()
    wb_times = read_wb_times(Path(report["sst_wb_output"]))
    iris_times = read_iris_times(Path(report["iris_output"]))

    if args.scan_indices.strip():
        scan_indices = [int(value.strip()) for value in args.scan_indices.split(",") if value.strip()]
    else:
        scan_indices = choose_scan_indices(report, args.num_samples)

    fig = plt.figure(figsize=(12, max(4, 4 * len(scan_indices))), dpi=160)
    matches = []

    for row, scan_index in enumerate(scan_indices, start=1):
        sst_time = wb_times[scan_index]
        iris_idx, abs_dt = nearest_iris_index(iris_times, sst_time)
        iris_map = make_iris_map(hmi_map, iris_idx, iris_path=Path(report["iris_output"]))
        try:
            sst_map, sst_corners, sst_label = build_sst_map(hmi_map, report, args.sst_kind, scan_index, Path(args.nb_path))
        except RuntimeError as exc:
            print(exc, file=sys.stderr)
            raise SystemExit(1)

        with plt.rc_context({"axes.grid": False}):
            extent = extent_from_footprint(sst_corners)
            ax_left = fig.add_subplot(len(scan_indices), 2, 2 * row - 1)
            plot_overlay(
                ax_left,
                hmi_map,
                iris_map,
                extent,
                f"IRIS vs HMI  |  IRIS frame {iris_idx}  |  {iris_times[iris_idx].isot}",
                overlay_alpha=0.12,
                hmi_alpha=0.78,
            )

            ax_right = fig.add_subplot(len(scan_indices), 2, 2 * row)
            plot_overlay(
                ax_right,
                hmi_map,
                sst_map,
                extent,
                f"SST {sst_label} vs HMI  |  scan {scan_index}  |  {sst_time.isot}",
                overlay_alpha=0.50,
                hmi_alpha=0.72,
            )

        matches.append(
            {
                "scan_index": int(scan_index),
                "sst_kind": args.sst_kind,
                "sst_time": sst_time.isot,
                "iris_frame_index": int(iris_idx),
                "iris_time": iris_times[iris_idx].isot,
                "time_offset_seconds": abs_dt,
            }
        )

    fig.tight_layout()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    Path(args.match_output).write_text(json.dumps(matches, indent=2))
    print(output_path)


if __name__ == "__main__":
    main()
