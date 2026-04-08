from __future__ import annotations

import argparse

from alignment_common import (
    DEFAULT_IRIS_ROTATION_DEG,
    DEFAULT_IRIS_SHIFT_ARCSEC,
    DEFAULT_SST_SHIFT_ARCSEC,
    IRIS_ALIGNED_PATH,
    REPORT_PATH,
    SST_WB_ALIGNED_PATH,
    ensure_initial_alignment,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the simple initial IRIS/SST alignment products used by the interactive app."
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Rebuild the aligned outputs from the original source FITS files.",
    )
    parser.add_argument("--sst-dx", type=float, default=DEFAULT_SST_SHIFT_ARCSEC[0], help="Initial SST X shift in arcsec.")
    parser.add_argument("--sst-dy", type=float, default=DEFAULT_SST_SHIFT_ARCSEC[1], help="Initial SST Y shift in arcsec.")
    parser.add_argument(
        "--iris-dx",
        type=float,
        default=DEFAULT_IRIS_SHIFT_ARCSEC[0],
        help="Initial IRIS X shift in arcsec.",
    )
    parser.add_argument(
        "--iris-dy",
        type=float,
        default=DEFAULT_IRIS_SHIFT_ARCSEC[1],
        help="Initial IRIS Y shift in arcsec.",
    )
    parser.add_argument(
        "--iris-rotation",
        type=float,
        default=DEFAULT_IRIS_ROTATION_DEG,
        help="Initial IRIS rotation in degrees.",
    )
    parser.add_argument(
        "--no-solve-iris-hmi",
        action="store_true",
        help="Skip the simple automatic IRIS-to-HMI solve and use the provided IRIS shift/rotation values directly.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = ensure_initial_alignment(
        overwrite=args.reset,
        sst_shift_arcsec=(args.sst_dx, args.sst_dy),
        iris_shift_arcsec=(args.iris_dx, args.iris_dy),
        iris_rotation_deg=args.iris_rotation,
        solve_iris_against_hmi=not args.no_solve_iris_hmi,
    )
    print("Initial alignment products are ready.")
    print(f"SST aligned FITS: {SST_WB_ALIGNED_PATH}")
    print(f"IRIS aligned FITS: {IRIS_ALIGNED_PATH}")
    print(f"Report: {REPORT_PATH}")
    print(
        "Initial shifts: "
        f"SST ({report['current_total_sst_shift_arcsec']['dx']:+.3f}, "
        f"{report['current_total_sst_shift_arcsec']['dy']:+.3f}) arcsec, "
        f"IRIS rotation {report['iris']['rotation_deg']:+.3f} deg."
    )


if __name__ == "__main__":
    main()
