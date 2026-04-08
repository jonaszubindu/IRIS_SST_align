from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/sunpy")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

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
        description="Plot IRIS/SST overlap checks for several time steps across the observation overlap."
    )
    parser.add_argument("--sst-kind", choices=["wb", "nb"], default="nb", help="Use WB or NB wing images for the SST panel.")
    parser.add_argument("--nb-path", default=str(SST_NB_SOURCE_PATH), help="NB cube to use when --sst-kind nb.")
    parser.add_argument("--scan-indices", default="", help="Comma-separated SST scan indices. If omitted, sample the overlap automatically.")
    parser.add_argument("--num-samples", type=int, default=4, help="Number of automatically sampled scan indices.")
    parser.add_argument(
        "--output",
        default=str(OUTPUT_DIR / "sst_iris_time_overlap.png"),
        help="Output PNG path.",
    )
    parser.add_argument(
        "--match-output",
        default=str(OUTPUT_DIR / "sst_iris_time_overlap_matches.json"),
        help="Output JSON file with SST/IRIS time matches.",
    )
    return parser.parse_args()


def choose_scan_indices(num_samples: int) -> list[int]:
    with open(REPORT_PATH) as fh:
        report = json.load(fh)
    wb_times = read_wb_times(Path(report["sst_wb_output"]))
    if num_samples <= 1:
        return [0]
    return np.linspace(0, len(wb_times) - 1, num_samples, dtype=int).tolist()


def nearest_iris_index(iris_times, target_time):
    offsets = np.abs((iris_times - target_time).sec)
    idx = int(np.argmin(offsets))
    return idx, float(offsets[idx])


def build_sst_map(hmi_map, sst_kind: str, scan_index: int, nb_path: Path, wb_path: Path):
    if sst_kind == "wb":
        sst_map, (corners, _) = make_wb_map(hmi_map, scan_index)
        label = "WB"
    else:
        sst_map, wing_indices = make_nb_wing_map(hmi_map, scan_index, nb_path=nb_path)
        corners = read_wb_corners(wb_path)[scan_index]
        label = f"NB wings {wing_indices[0]} / {wing_indices[1]}"
    return sst_map, corners, label


def extent_from_corners(corners: np.ndarray, margin_arcsec: float = 6.0) -> list[float]:
    x0 = float(corners[:, 0, 0].mean())
    x1 = float(corners[:, 1, 0].mean())
    y0 = float(corners[0, :, 1].mean())
    y1 = float(corners[1, :, 1].mean())
    return [
        min(x0, x1) - margin_arcsec,
        max(x0, x1) + margin_arcsec,
        min(y0, y1) - margin_arcsec,
        max(y0, y1) + margin_arcsec,
    ]


def main() -> None:
    args = parse_args()
    hmi_map = load_hmi_map()
    with open(REPORT_PATH) as fh:
        report = json.load(fh)
    wb_times = read_wb_times(Path(report["sst_wb_output"]))
    iris_times = read_iris_times(Path(report["iris_output"]))

    if args.scan_indices.strip():
        scan_indices = [int(value.strip()) for value in args.scan_indices.split(",") if value.strip()]
    else:
        scan_indices = choose_scan_indices(args.num_samples)

    n_panels = len(scan_indices)
    fig = plt.figure(figsize=(8, max(4, 4 * n_panels)), dpi=160)
    matches = []

    for row, scan_index in enumerate(scan_indices, start=1):
        sst_time = wb_times[scan_index]
        iris_idx, abs_dt = nearest_iris_index(iris_times, sst_time)
        iris_map = make_iris_map(hmi_map, iris_idx)
        try:
            sst_map, corners, sst_label = build_sst_map(
                hmi_map,
                args.sst_kind,
                scan_index,
                Path(args.nb_path),
                Path(report["sst_wb_output"]),
            )
        except RuntimeError as exc:
            print(exc, file=sys.stderr)
            raise SystemExit(1)

        extent = extent_from_corners(corners)
        shape = (800, 800)
        iris_image = normalize_for_display(project_map_to_extent(iris_map, extent, shape))
        sst_image = normalize_for_display(project_map_to_extent(sst_map, extent, shape))

        ax = fig.add_subplot(n_panels, 1, row)
        ax.imshow(iris_image, extent=extent, origin="lower", cmap="binary_r", alpha=0.5, aspect="equal")
        ax.imshow(sst_image, extent=extent, origin="lower", cmap="gray", alpha=0.5, aspect="equal")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("Helioprojective Longitude (Solar-X)")
        ax.set_ylabel("Helioprojective Latitude (Solar-Y)")
        ax.set_title(
            f"Scan {scan_index} ({sst_label})  |  SST {sst_time.isot}  |  IRIS frame {iris_idx} {iris_times[iris_idx].isot}  |  |dt|={abs_dt:.1f}s"
        )

        matches.append(
            {
                "scan_index": int(scan_index),
                "sst_time": sst_time.isot,
                "iris_frame_index": int(iris_idx),
                "iris_time": iris_times[iris_idx].isot,
                "time_offset_seconds": abs_dt,
                "sst_kind": args.sst_kind,
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
