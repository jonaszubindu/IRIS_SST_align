from __future__ import annotations

import os

os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/sunpy")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from alignment_common import (
    OUTPUT_DIR,
    REPORT_PATH,
    load_hmi_map,
    load_report,
    make_iris_map,
    make_wb_map,
    normalize_for_display,
    project_map_to_extent,
)


GIF_PATH = OUTPUT_DIR / "alignment_check_sunpy.gif"


def main() -> None:
    report = load_report()
    hmi_map = load_hmi_map()
    iris_idx = int(report["iris"]["frame_index"])
    sst_idx = int(report["sst_wb"]["frame_index"])
    sji_map = make_iris_map(hmi_map, iris_idx)
    sst_map, (corners, _) = make_wb_map(hmi_map, sst_idx)

    x0 = float(corners[:, 0, 0].mean())
    x1 = float(corners[:, 1, 0].mean())
    y0 = float(corners[0, :, 1].mean())
    y1 = float(corners[1, :, 1].mean())
    margin = 6.0
    extent = [
        min(x0, x1) - margin,
        max(x0, x1) + margin,
        min(y0, y1) - margin,
        max(y0, y1) + margin,
    ]
    shape = (900, 900)
    iris_image = normalize_for_display(project_map_to_extent(sji_map, extent, shape))
    sst_image = normalize_for_display(project_map_to_extent(sst_map, extent, shape))

    fig = plt.figure(figsize=(8, 8), dpi=120)
    ax = fig.add_subplot(111)

    iris_im = ax.imshow(iris_image, extent=extent, origin="lower", cmap="binary_r", alpha=1.0, aspect="equal")
    sst_im = ax.imshow(sst_image, extent=extent, origin="lower", cmap="gray", alpha=0.0, aspect="equal")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("Helioprojective Longitude (Solar-X)")
    ax.set_ylabel("Helioprojective Latitude (Solar-Y)")
    ax.set_title("Blinking alignment check")

    frame_paths = []
    for idx, (a1, a2) in enumerate([(0.3, 0.7), (0.7, 0.3)]):
        iris_im.set_alpha(a1)
        sst_im.set_alpha(a2)
        frame_path = OUTPUT_DIR / f"alignment_check_sunpy_frame_{idx}.png"
        fig.savefig(frame_path, dpi=120)
        frame_paths.append(frame_path)

    plt.close(fig)
    frames = []
    for path in frame_paths:
        with Image.open(path) as image:
            frames.append(image.convert("P"))
    frames[0].save(
        GIF_PATH,
        save_all=True,
        append_images=[frames[1], frames[0], frames[1]],
        duration=1000,
        loop=0,
    )
    print(GIF_PATH)


if __name__ == "__main__":
    main()
