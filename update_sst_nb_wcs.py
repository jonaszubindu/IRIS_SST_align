from __future__ import annotations

import argparse
import sys

from alignment_common import SST_NB_ALIGNED_PATH, SST_NB_SOURCE_PATH, SST_WB_ALIGNED_PATH, update_nb_cube_wcs_from_wb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Propagate the aligned SST wideband spatial WCS into an SST narrowband cube."
    )
    parser.add_argument("--nb-path", default=str(SST_NB_SOURCE_PATH), help="Input NB cube.")
    parser.add_argument("--wb-path", default=str(SST_WB_ALIGNED_PATH), help="Aligned WB cube used as the WCS source.")
    parser.add_argument("--output", default=str(SST_NB_ALIGNED_PATH), help="Output aligned NB cube.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        output_path = update_nb_cube_wcs_from_wb(args.nb_path, args.wb_path, args.output)
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        raise SystemExit(1)
    print(f"Aligned NB cube written to {output_path}")


if __name__ == "__main__":
    main()
