from __future__ import annotations

import argparse

from alignment_common import DEFAULT_JSOC_EMAIL, DEFAULT_NOAA_AR, HMI_PATH, ensure_hmi_fits


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ensure the HMI reference FITS exists locally, downloading it from JSOC through drms if needed."
    )
    parser.add_argument(
        "--email",
        default=DEFAULT_JSOC_EMAIL,
        help="Registered JSOC email address. Defaults to IRIS_SST_ALIGN_JSOC_EMAIL or JSOC_EMAIL.",
    )
    parser.add_argument(
        "--noaa-ar",
        default=DEFAULT_NOAA_AR,
        help="NOAA active-region number used to build the HMI JSOC query.",
    )
    parser.add_argument(
        "--output",
        default=str(HMI_PATH),
        help="Output path for the downloaded HMI FITS.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hmi_path = ensure_hmi_fits(output_path=args.output, email=args.email, noaa_ar=args.noaa_ar)
    print(f"HMI reference FITS is available at {hmi_path}")


if __name__ == "__main__":
    main()
