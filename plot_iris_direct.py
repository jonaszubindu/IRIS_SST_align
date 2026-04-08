from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/sunpy")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from alignment_common import OUTPUT_DIR, load_hmi_map, load_report, make_iris_map, normalize_for_display, project_map_to_extent


def main() -> None:
    out = OUTPUT_DIR / "iris_map_direct.png"
    report = load_report()
    hmi_map = load_hmi_map()
    iris_map = make_iris_map(hmi_map, int(report["iris"]["frame_index"]), Path(report["iris_output"]))
    center = iris_map.center
    extent = [
        float(iris_map.bottom_left_coord.Tx.arcsec),
        float(iris_map.top_right_coord.Tx.arcsec),
        float(iris_map.bottom_left_coord.Ty.arcsec),
        float(iris_map.top_right_coord.Ty.arcsec),
    ]
    image = normalize_for_display(project_map_to_extent(iris_map, extent, iris_map.data.shape))

    fig = plt.figure(figsize=(7, 7), dpi=180)
    ax = fig.add_subplot(111)
    ax.imshow(image, extent=extent, origin="lower", cmap="gray", aspect="equal")
    ax.plot(center.Tx.arcsec, center.Ty.arcsec, marker="+", color="red", markersize=18, markeredgewidth=2.2)
    ax.set_xlabel("Helioprojective Longitude (Solar-X)")
    ax.set_ylabel("Helioprojective Latitude (Solar-Y)")
    ax.set_title(f"IRIS direct SunPy map | frame {report['iris']['frame_index']}")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(out)
    print("center", center.Tx.arcsec, center.Ty.arcsec)


if __name__ == "__main__":
    main()
