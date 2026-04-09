from __future__ import annotations

import json
import math
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta


WORKDIR = Path(__file__).resolve().parent
OUTPUT_DIR = WORKDIR / "outputs"
DATA_HMI_DIR = WORKDIR / "data/hmi"
DEFAULT_DATA_DIR = Path("/Users/jonaszbinden/Desktop/Align_IRIS_SST_proj")
DATA_DIR = Path(os.environ.get("IRIS_SST_ALIGN_DATA_DIR", str(DEFAULT_DATA_DIR))).expanduser()

DEFAULT_HMI_FILENAME = "hmi.sharp_720s.reference.continuum.fits"
EXPLICIT_HMI_PATH = os.environ.get("IRIS_SST_ALIGN_HMI_PATH", "").strip()
HMI_PATH = Path(EXPLICIT_HMI_PATH).expanduser() if EXPLICIT_HMI_PATH else (DATA_HMI_DIR / DEFAULT_HMI_FILENAME)
IRIS_SOURCE_PATH = DATA_DIR / "iris_l2_20250619_072925_3660106834_SJI_2832_t000.fits"
SST_WB_SOURCE_PATH = DATA_DIR / "wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im.fits"
SST_NB_SOURCE_PATH = DATA_DIR / "nb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_cmapcorr_im.fits"
STIC_UTILS_PATH = DATA_DIR / "STIC_pyfun_utils.py"

DEFAULT_JSOC_EMAIL = os.environ.get("IRIS_SST_ALIGN_JSOC_EMAIL", os.environ.get("JSOC_EMAIL", "")).strip()
DEFAULT_NOAA_AR = os.environ.get("IRIS_SST_ALIGN_NOAA_AR", "").strip()
DEFAULT_HMI_SERIES_TEMPLATE = "hmi.sharp_720s[][{time}_TAI/1m][?(NOAA_AR={noaa_ar})?]"
DEFAULT_HMI_SEGMENT = "continuum"

IRIS_ALIGNED_PATH = OUTPUT_DIR / "iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits"
SST_WB_ALIGNED_PATH = OUTPUT_DIR / "wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits"
SST_NB_ALIGNED_PATH = OUTPUT_DIR / "nb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_cmapcorr_im_aligned.fits"
REPORT_PATH = OUTPUT_DIR / "alignment_report.json"
MANUAL_SAVE_PATH = OUTPUT_DIR / "sst_manual_adjustment.json"

DEFAULT_SST_SHIFT_ARCSEC = (3.602630034836494, -3.5941129840652195)
DEFAULT_IRIS_SHIFT_ARCSEC = (0.0, 0.0)
DEFAULT_IRIS_ROTATION_DEG = 0.0


def ensure_output_dir() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _format_gib(num_bytes: int) -> str:
    return f"{num_bytes / (1024 ** 3):.1f} GiB"


def copy_file_with_progress(
    source_path: Path,
    output_path: Path,
    chunk_size: int = 64 * 1024 * 1024,
) -> None:
    source_path = Path(source_path)
    output_path = Path(output_path)
    total_bytes = source_path.stat().st_size
    copied = 0
    last_percent = -1

    print(f"Copying NB cube: {_format_gib(total_bytes)} total", file=sys.stderr, flush=True)
    with source_path.open("rb") as src, output_path.open("wb") as dst:
        while True:
            chunk = src.read(chunk_size)
            if not chunk:
                break
            dst.write(chunk)
            copied += len(chunk)
            percent = int((100 * copied) / total_bytes) if total_bytes else 100
            if percent != last_percent:
                filled = min(40, int(percent * 40 / 100))
                bar = "#" * filled + "-" * (40 - filled)
                print(
                    f"\r[{bar}] {percent:3d}%  {_format_gib(copied)} / {_format_gib(total_bytes)}",
                    end="",
                    file=sys.stderr,
                    flush=True,
                )
                last_percent = percent
    shutil.copystat(source_path, output_path)
    print(file=sys.stderr, flush=True)


def ensure_hmi_dir() -> None:
    HMI_PATH.parent.mkdir(parents=True, exist_ok=True)


def infer_reference_time(
    iris_path: Path = IRIS_SOURCE_PATH,
    wb_path: Path = SST_WB_SOURCE_PATH,
) -> Time:
    wb_path = Path(wb_path)
    if wb_path.exists():
        return read_wb_times(wb_path)[0]
    iris_path = Path(iris_path)
    if iris_path.exists():
        return read_iris_times(iris_path)[0]
    raise RuntimeError(
        "Could not infer the HMI reference time because neither the SST WB source nor the IRIS source file exists."
    )


def format_jsoc_tai_time(reference_time: Time) -> str:
    return reference_time.tai.strftime("%Y.%m.%d_%H:%M:%S")


def build_hmi_export_query(
    noaa_ar: str | int,
    reference_time: Time,
    segment: str = DEFAULT_HMI_SEGMENT,
) -> str:
    noaa_ar_str = str(noaa_ar).strip()
    if not noaa_ar_str:
        raise RuntimeError(
            "A NOAA active-region number is required to auto-download HMI. "
            "Set IRIS_SST_ALIGN_NOAA_AR or pass --noaa-ar."
        )
    series = DEFAULT_HMI_SERIES_TEMPLATE.format(
        time=format_jsoc_tai_time(reference_time),
        noaa_ar=noaa_ar_str,
    )
    return f"{series}{{{segment}}}"


def default_hmi_output_path(
    noaa_ar: str | int,
    reference_time: Time,
    segment: str = DEFAULT_HMI_SEGMENT,
) -> Path:
    timestamp = reference_time.tai.strftime("%Y%m%d_%H%M%S")
    return DATA_HMI_DIR / f"hmi.sharp_720s.noaa{str(noaa_ar).strip()}.{timestamp}_TAI.{segment}.fits"


def export_hmi_reference(
    email: str,
    output_path: Path | None = None,
    noaa_ar: str | int | None = None,
    export_query: str | None = None,
    reference_time: Time | None = None,
    segment: str = DEFAULT_HMI_SEGMENT,
) -> Path:
    try:
        import drms
    except ImportError as exc:  # pragma: no cover - dependency/runtime specific
        raise RuntimeError(
            "The HMI reference file is missing and the 'drms' package is not installed. "
            "Install the package dependencies or place the HMI FITS at the configured HMI path."
        ) from exc

    if not email:
        raise RuntimeError(
            "The HMI reference file is missing and no JSOC email is configured. "
            "Set IRIS_SST_ALIGN_JSOC_EMAIL or JSOC_EMAIL, or place the HMI FITS at the configured HMI path."
        )

    resolved_noaa_ar = str(noaa_ar or DEFAULT_NOAA_AR).strip()
    if export_query is None:
        if reference_time is None:
            reference_time = infer_reference_time()
        export_query = build_hmi_export_query(resolved_noaa_ar, reference_time, segment=segment)
    if output_path is None:
        if EXPLICIT_HMI_PATH:
            output_path = HMI_PATH
        else:
            if reference_time is None:
                reference_time = infer_reference_time()
            output_path = default_hmi_output_path(resolved_noaa_ar, reference_time, segment=segment)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    client = drms.Client(email=email)
    export_request = client.export(export_query, method="url", protocol="fits")
    if hasattr(export_request, "wait"):
        export_request.wait()

    if hasattr(export_request, "download"):
        export_request.download(output_path.parent)
    else:  # pragma: no cover - defensive API fallback
        raise RuntimeError("DRMS export request does not support direct download in this environment.")

    downloaded = sorted(output_path.parent.glob("*.fits"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not downloaded:
        raise RuntimeError(f"JSOC export succeeded but no FITS file was downloaded into {output_path.parent}.")

    candidate = downloaded[0]
    if candidate.resolve() != output_path.resolve():
        if output_path.exists():
            output_path.unlink()
        candidate.replace(output_path)
    return output_path


def ensure_hmi_fits(
    output_path: Path | None = None,
    email: str | None = None,
    noaa_ar: str | int | None = None,
    export_query: str | None = None,
    reference_time: Time | None = None,
    segment: str = DEFAULT_HMI_SEGMENT,
) -> Path:
    if output_path is None:
        if EXPLICIT_HMI_PATH:
            output_path = HMI_PATH
        elif HMI_PATH.exists():
            output_path = HMI_PATH
        elif REPORT_PATH.exists():
            report = load_report()
            report_hmi_value = str(report.get("reference_hmi", "")).strip()
            report_hmi = Path(report_hmi_value).expanduser() if report_hmi_value else None
            if report_hmi and report_hmi.exists():
                output_path = report_hmi
        elif DATA_HMI_DIR.exists():
            existing_hmi = sorted(DATA_HMI_DIR.glob("*.fits"), key=lambda path: path.stat().st_mtime, reverse=True)
            if existing_hmi:
                output_path = existing_hmi[0]
        else:
            resolved_noaa_ar = str(noaa_ar or DEFAULT_NOAA_AR).strip()
            if resolved_noaa_ar:
                if reference_time is None:
                    reference_time = infer_reference_time()
                output_path = default_hmi_output_path(resolved_noaa_ar, reference_time, segment=segment)
            else:
                output_path = HMI_PATH
    output_path = Path(output_path)
    if output_path.exists():
        return output_path
    resolved_email = (email or DEFAULT_JSOC_EMAIL).strip()
    return export_hmi_reference(
        resolved_email,
        output_path=output_path,
        noaa_ar=noaa_ar,
        export_query=export_query,
        reference_time=reference_time,
        segment=segment,
    )


def load_report() -> dict:
    if REPORT_PATH.exists():
        return json.loads(REPORT_PATH.read_text())
    return {}


def save_report(report: dict) -> None:
    ensure_output_dir()
    REPORT_PATH.write_text(json.dumps(report, indent=2))


def file_is_truncated_primary_hdu(path: Path) -> bool:
    path = Path(path)
    hdr = fits.getheader(path, 0)
    header_bytes = len(hdr.tostring(sep=""))
    bitpix = int(hdr["BITPIX"])
    n_dim = int(hdr["NAXIS"])
    data_bytes = abs(bitpix) // 8
    for axis in range(1, n_dim + 1):
        data_bytes *= int(hdr[f"NAXIS{axis}"])
    padded_data_bytes = int(math.ceil(data_bytes / 2880.0) * 2880) if data_bytes else 0
    expected_primary_end = header_bytes + padded_data_bytes
    return path.stat().st_size < expected_primary_end


def read_wb_tab(path: Path) -> np.ndarray:
    with fits.open(path, memmap=False) as hdul:
        return np.asarray(hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0], dtype=np.float64)


def read_wb_corners(path: Path) -> np.ndarray:
    return read_wb_tab(path)[..., :2]


def read_wb_times_seconds(path: Path) -> np.ndarray:
    return read_wb_tab(path)[:, 0, 0, 2]


def double2time(var: float) -> str:
    hh = int(np.floor(var / 3600))
    rest = var - hh * 3600
    mm = int(np.floor(rest / 60))
    rest -= mm * 60
    ss = var - hh * 3600 - mm * 60
    return ":".join(("{0:02d}".format(hh), "{0:02d}".format(mm), "{0:06.3f}".format(ss)))


def read_wb_times(path: Path) -> Time:
    with fits.open(path, memmap=False) as hdul:
        obs_date = hdul[0].header["DATE-OBS"].split("T")[0]
    sod_seconds = read_wb_times_seconds(path)
    isot = [f"{obs_date}T{double2time(float(sec))}" for sec in sod_seconds]
    return Time(isot, format="isot")


def read_iris_times(path: Path) -> Time:
    with fits.open(path, memmap=False) as hdul:
        base_time = Time(hdul[0].header["DATE_OBS"])
        aux = np.asarray(hdul[1].data, dtype=np.float64)
    return base_time + TimeDelta(aux[:, 0], format="sec")


def closest_iris_frame(target_time: Time, iris_path: Path = IRIS_ALIGNED_PATH) -> tuple[int, Time]:
    iris_times = read_iris_times(iris_path)
    dt = np.abs((iris_times - target_time).sec)
    idx = int(np.argmin(dt))
    return idx, iris_times[idx]


def normalize_alignment_image(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(image)
    if not finite.any():
        return np.zeros_like(image, dtype=np.float32)
    fill = np.nanmedian(image[finite])
    image = np.where(finite, image, fill)
    p1, p99 = np.percentile(image, [1, 99])
    image = np.clip(image, p1, p99)
    image = image - image.mean()
    std = image.std()
    return image / std if std > 0 else image


def project_data_wcs_to_extent(data: np.ndarray, wcs_in: WCS, extent: list[float], shape: tuple[int, int]) -> np.ndarray:
    from scipy.ndimage import map_coordinates

    ny, nx = shape
    y, x = np.indices((ny, nx), dtype=np.float64)
    lon_arcsec = extent[0] + (extent[1] - extent[0]) * x / max(nx - 1, 1)
    lat_arcsec = extent[2] + (extent[3] - extent[2]) * y / max(ny - 1, 1)
    cunit1 = u.Unit(wcs_in.wcs.cunit[0] or "deg")
    cunit2 = u.Unit(wcs_in.wcs.cunit[1] or "deg")
    lon = (lon_arcsec * u.arcsec).to_value(cunit1)
    lat = (lat_arcsec * u.arcsec).to_value(cunit2)
    sx, sy = wcs_in.world_to_pixel_values(lon, lat)
    sampled = map_coordinates(
        np.asarray(data, dtype=np.float32),
        [sy, sx],
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    return sampled.astype(np.float32)


def project_map_to_extent(map_obj, extent: list[float], shape: tuple[int, int]) -> np.ndarray:
    return project_data_wcs_to_extent(
        np.asarray(map_obj.data, dtype=np.float32),
        map_obj.wcs.celestial,
        extent,
        shape,
    )


def world_bounds_from_header(header: fits.Header) -> tuple[float, float, float, float]:
    wcs = WCS(header).celestial
    nx = int(header["NAXIS1"])
    ny = int(header["NAXIS2"])
    px = np.array([0.0, nx - 1.0, 0.0, nx - 1.0], dtype=np.float64)
    py = np.array([0.0, 0.0, ny - 1.0, ny - 1.0], dtype=np.float64)
    lon, lat = wcs.pixel_to_world_values(px, py)
    return float(np.nanmin(lon)), float(np.nanmax(lon)), float(np.nanmin(lat)), float(np.nanmax(lat))


def world_bounds_from_corners(corners: np.ndarray) -> tuple[float, float, float, float]:
    lon = np.asarray(corners[..., 0], dtype=np.float64)
    lat = np.asarray(corners[..., 1], dtype=np.float64)
    return float(np.nanmin(lon)), float(np.nanmax(lon)), float(np.nanmin(lat)), float(np.nanmax(lat))


def square_extent_from_bounds(bounds: tuple[float, float, float, float], pad_arcsec: float = 12.0) -> list[float]:
    xmin, xmax, ymin, ymax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    half = 0.5 * max(xmax - xmin, ymax - ymin) + pad_arcsec
    return [cx - half, cx + half, cy - half, cy + half]


def project_hmi_to_iris_header(hmi_map, iris_header: fits.Header) -> np.ndarray:
    from scipy.ndimage import map_coordinates

    ny = int(iris_header["NAXIS2"])
    nx = int(iris_header["NAXIS1"])
    y, x = np.indices((ny, nx), dtype=np.float64)
    iris_wcs = WCS(iris_header).celestial
    lon, lat = iris_wcs.pixel_to_world_values(x, y)
    hx, hy = hmi_map.wcs.world_to_pixel_values(lon, lat)
    sampled = map_coordinates(
        np.asarray(hmi_map.data, dtype=np.float32),
        [hy, hx],
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    return sampled.astype(np.float32)


def solve_iris_to_hmi(
    hmi_map,
    iris_path: Path = IRIS_ALIGNED_PATH,
    frame_index: int | None = None,
    rotation_grid_deg: np.ndarray | None = None,
) -> dict:
    from scipy.ndimage import rotate as ndi_rotate, shift as ndi_shift
    from skimage.registration import phase_cross_correlation

    report = load_report()
    if frame_index is None:
        frame_index = int(report["iris"]["frame_index"])

    with fits.open(iris_path, memmap=False) as hdul:
        iris_data = np.asarray(hdul[0].data[frame_index], dtype=np.float32)
        iris_header = hdul[0].header.copy()

    if rotation_grid_deg is None:
        rotation_grid_deg = np.arange(-0.6, 0.61, 0.05, dtype=np.float64)

    hmi_on_iris = project_hmi_to_iris_header(hmi_map, iris_header)
    row_min, row_max, col_min, col_max = active_region_bounds(-normalize_alignment_image(hmi_on_iris), margin=70)
    ref = normalize_alignment_image(hmi_on_iris[row_min : row_max + 1, col_min : col_max + 1])
    best: dict | None = None

    for theta in rotation_grid_deg:
        rotated = ndi_rotate(iris_data, float(theta), reshape=False, order=1, mode="constant", cval=np.nan, prefilter=False)
        mov_full = normalize_alignment_image(rotated)
        mov = mov_full[row_min : row_max + 1, col_min : col_max + 1]
        shift, error, _ = phase_cross_correlation(
            np.nan_to_num(ref, nan=0.0),
            np.nan_to_num(mov, nan=0.0),
            upsample_factor=20,
        )
        coarse_y = np.arange(shift[0] - 2.0, shift[0] + 2.01, 0.5)
        coarse_x = np.arange(shift[1] - 2.0, shift[1] + 2.01, 0.5)
        local_best: dict | None = None
        for dy in coarse_y:
            for dx in coarse_x:
                shifted = ndi_shift(mov, shift=(dy, dx), order=1, mode="constant", cval=np.nan, prefilter=False)
                mask = np.isfinite(ref) & np.isfinite(shifted)
                if mask.sum() < 5000:
                    continue
                score = float(np.corrcoef(ref[mask].ravel(), shifted[mask].ravel())[0, 1])
                candidate = {
                    "rotation_deg": float(theta),
                    "pixel_shift": [float(dx), float(dy)],
                    "world_shift_arcsec": [
                        float(dx * iris_header["CDELT1"]),
                        float(dy * iris_header["CDELT2"]),
                    ],
                    "score": score,
                    "error": float(error),
                    "frame_index": int(frame_index),
                }
                if local_best is None or candidate["score"] > local_best["score"]:
                    local_best = candidate
        if local_best and (best is None or local_best["score"] > best["score"]):
            best = local_best

    assert best is not None
    return best


def build_rotation_pc(angle_deg: float) -> tuple[float, float, float, float]:
    theta = math.radians(angle_deg)
    return math.cos(theta), -math.sin(theta), math.sin(theta), math.cos(theta)


def apply_iris_shift_and_rotation(
    source_path: Path,
    output_path: Path,
    shift_arcsec: tuple[float, float],
    rotation_deg: float,
    overwrite: bool,
) -> None:
    if output_path.exists() and not overwrite:
        return
    shutil.copy2(source_path, output_path)
    with fits.open(output_path, mode="update", memmap=False) as hdul:
        hdr = hdul[0].header
        hdr["CRVAL1"] = float(hdr["CRVAL1"]) + float(shift_arcsec[0])
        hdr["CRVAL2"] = float(hdr["CRVAL2"]) + float(shift_arcsec[1])
        pc11, pc12, pc21, pc22 = build_rotation_pc(rotation_deg)
        hdr["PC1_1"] = pc11
        hdr["PC1_2"] = pc12
        hdr["PC2_1"] = pc21
        hdr["PC2_2"] = pc22
        hdr["ALNROTI"] = float(rotation_deg)
        hdr["ALNIRX"] = float(shift_arcsec[0])
        hdr["ALNIRY"] = float(shift_arcsec[1])
        hdul.flush()


def apply_sst_shift(
    source_path: Path,
    output_path: Path,
    shift_arcsec: tuple[float, float],
    overwrite: bool,
) -> None:
    if output_path.exists() and not overwrite:
        return
    shutil.copy2(source_path, output_path)
    with fits.open(output_path, mode="update", memmap=False) as hdul:
        tab = np.asarray(hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0], dtype=np.float64)
        tab[..., 0] += float(shift_arcsec[0])
        tab[..., 1] += float(shift_arcsec[1])
        hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0] = tab
        hdr = hdul[0].header
        hdr["ALNWX"] = float(shift_arcsec[0])
        hdr["ALNWY"] = float(shift_arcsec[1])
        hdr["ALNMANX"] = 0.0
        hdr["ALNMANY"] = 0.0
        hdul.flush()


def create_initial_report(
    sst_shift_arcsec: tuple[float, float] = DEFAULT_SST_SHIFT_ARCSEC,
    iris_shift_arcsec: tuple[float, float] = DEFAULT_IRIS_SHIFT_ARCSEC,
    iris_rotation_deg: float = DEFAULT_IRIS_ROTATION_DEG,
    iris_frame_index: int | None = None,
    reference_hmi_path: Path | None = None,
    noaa_ar: str | int | None = None,
) -> dict:
    sst_times = read_wb_times(SST_WB_ALIGNED_PATH)
    sst_frame_index = 0
    sst_time = sst_times[sst_frame_index]
    iris_time_path = IRIS_ALIGNED_PATH if IRIS_ALIGNED_PATH.exists() else IRIS_SOURCE_PATH
    if iris_frame_index is None:
        iris_frame_index, iris_time = closest_iris_frame(sst_time, iris_time_path)
    else:
        iris_time = read_iris_times(iris_time_path)[int(iris_frame_index)]
    return {
        "reference_hmi": str(reference_hmi_path or HMI_PATH),
        "noaa_active_region": str(noaa_ar or DEFAULT_NOAA_AR),
        "iris_output": str(IRIS_ALIGNED_PATH),
        "sst_wb_output": str(SST_WB_ALIGNED_PATH),
        "sst_nb_output": str(SST_NB_ALIGNED_PATH),
        "sst_nb_source": str(SST_NB_SOURCE_PATH),
        "initial_sst_shift_arcsec": {"dx": float(sst_shift_arcsec[0]), "dy": float(sst_shift_arcsec[1])},
        "manual_adjustment_arcsec": {"dx": 0.0, "dy": 0.0},
        "current_total_sst_shift_arcsec": {"dx": float(sst_shift_arcsec[0]), "dy": float(sst_shift_arcsec[1])},
        "iris": {
            "frame_index": int(iris_frame_index),
            "observation_time": iris_time.isot,
            "world_shift_arcsec": [float(iris_shift_arcsec[0]), float(iris_shift_arcsec[1])],
            "rotation_deg": float(iris_rotation_deg),
        },
        "sst_wb": {
            "frame_index": int(sst_frame_index),
            "observation_time": sst_time.isot,
            "world_shift_arcsec": [float(sst_shift_arcsec[0]), float(sst_shift_arcsec[1])],
        },
        "notes": [
            "Initial alignment uses a single global SST shift and optional IRIS rotation.",
            "Further alignment is intended to be done manually with the interactive app.",
        ],
    }


def ensure_initial_alignment(
    overwrite: bool = False,
    sst_shift_arcsec: tuple[float, float] = DEFAULT_SST_SHIFT_ARCSEC,
    iris_shift_arcsec: tuple[float, float] = DEFAULT_IRIS_SHIFT_ARCSEC,
    iris_rotation_deg: float = DEFAULT_IRIS_ROTATION_DEG,
    solve_iris_against_hmi: bool = True,
    noaa_ar: str | int | None = None,
) -> dict:
    ensure_output_dir()
    apply_sst_shift(SST_WB_SOURCE_PATH, SST_WB_ALIGNED_PATH, sst_shift_arcsec, overwrite=overwrite)
    solved_iris_frame = None
    hmi_path = ensure_hmi_fits(noaa_ar=noaa_ar) if solve_iris_against_hmi else HMI_PATH
    if solve_iris_against_hmi:
        tmp_header_report = create_initial_report(
            sst_shift_arcsec=sst_shift_arcsec,
            iris_shift_arcsec=(0.0, 0.0),
            iris_rotation_deg=0.0,
            reference_hmi_path=hmi_path,
            noaa_ar=noaa_ar,
        )
        save_report(tmp_header_report)
        hmi_map = load_hmi_map(noaa_ar=noaa_ar, hmi_path=hmi_path)
        iris_solution = solve_iris_to_hmi(hmi_map, iris_path=IRIS_SOURCE_PATH, frame_index=tmp_header_report["iris"]["frame_index"])
        iris_shift_arcsec = tuple(iris_solution["world_shift_arcsec"])
        iris_rotation_deg = float(iris_solution["rotation_deg"])
        solved_iris_frame = int(iris_solution["frame_index"])
    apply_iris_shift_and_rotation(
        IRIS_SOURCE_PATH,
        IRIS_ALIGNED_PATH,
        iris_shift_arcsec,
        iris_rotation_deg,
        overwrite=overwrite,
    )
    report = create_initial_report(
        sst_shift_arcsec=sst_shift_arcsec,
        iris_shift_arcsec=iris_shift_arcsec,
        iris_rotation_deg=iris_rotation_deg,
        iris_frame_index=solved_iris_frame,
        reference_hmi_path=hmi_path,
        noaa_ar=noaa_ar,
    )
    if solve_iris_against_hmi:
        report["iris"]["solved_against_hmi"] = True
    save_report(report)
    if overwrite and MANUAL_SAVE_PATH.exists():
        MANUAL_SAVE_PATH.unlink()
    return report


def load_hmi_map(noaa_ar: str | int | None = None, hmi_path: Path | None = None):
    import sunpy.map

    resolved_path = ensure_hmi_fits(output_path=hmi_path, noaa_ar=noaa_ar)
    return sunpy.map.Map(resolved_path)


def make_iris_map(hmi_map, frame_index: int, iris_path: Path = IRIS_ALIGNED_PATH):
    import sunpy.map

    with fits.open(iris_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data[frame_index], dtype=np.float32)
        src = hdul[0].header
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = data.shape[1]
    header["NAXIS2"] = data.shape[0]
    for key in [
        "CTYPE1",
        "CTYPE2",
        "CUNIT1",
        "CUNIT2",
        "CRPIX1",
        "CRPIX2",
        "CRVAL1",
        "CRVAL2",
        "CDELT1",
        "CDELT2",
        "PC1_1",
        "PC1_2",
        "PC2_1",
        "PC2_2",
    ]:
        if key in src:
            header[key] = src[key]
    header["DATE-OBS"] = src.get("DATE_OBS", src.get("DATE-OBS"))
    header["TELESCOP"] = src.get("TELESCOP", "IRIS")
    header["INSTRUME"] = src.get("INSTRUME", "SJI")
    header["WAVELNTH"] = 2832
    header["WAVEUNIT"] = "Angstrom"
    header["BUNIT"] = src.get("BUNIT", "Corrected DN")
    for dst_key, src_key in [
        ("RSUN_OBS", "rsun_obs"),
        ("DSUN_OBS", "dsun_obs"),
        ("HGLN_OBS", "crln_obs"),
        ("HGLT_OBS", "crlt_obs"),
    ]:
        header[dst_key] = hmi_map.meta.get(src_key, hmi_map.meta.get(dst_key.lower(), 0.0))
    header["HGln_OBS"] = header["HGLN_OBS"]
    header["HGlt_OBS"] = header["HGLT_OBS"]
    return sunpy.map.Map(data, header)


def corners_to_tan_header(corners: np.ndarray, data_shape: tuple[int, int], src_header: fits.Header, hmi_map) -> fits.Header:
    ny, nx = data_shape
    x0 = float(np.mean(corners[:, 0, 0]))
    x1 = float(np.mean(corners[:, 1, 0]))
    y0 = float(np.mean(corners[0, :, 1]))
    y1 = float(np.mean(corners[1, :, 1]))
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = nx
    header["NAXIS2"] = ny
    header["CTYPE1"] = "HPLN-TAN"
    header["CTYPE2"] = "HPLT-TAN"
    header["CUNIT1"] = "arcsec"
    header["CUNIT2"] = "arcsec"
    header["CRPIX1"] = (nx + 1) / 2.0
    header["CRPIX2"] = (ny + 1) / 2.0
    header["CRVAL1"] = 0.5 * (x0 + x1)
    header["CRVAL2"] = 0.5 * (y0 + y1)
    header["CDELT1"] = (x1 - x0) / (nx - 1)
    header["CDELT2"] = (y1 - y0) / (ny - 1)
    header["DATE-OBS"] = src_header["DATE-OBS"]
    header["TELESCOP"] = src_header.get("TELESCOP", "SST")
    header["INSTRUME"] = src_header.get("INSTRUME", "CHROMIS")
    header["WAVELNTH"] = int(round(float(src_header.get("WAVELNTH", 395.0)) * 10))
    header["WAVEUNIT"] = "Angstrom"
    header["BUNIT"] = "arb. unit"
    for dst_key, src_key in [
        ("RSUN_OBS", "rsun_obs"),
        ("DSUN_OBS", "dsun_obs"),
        ("HGLN_OBS", "crln_obs"),
        ("HGLT_OBS", "crlt_obs"),
    ]:
        header[dst_key] = hmi_map.meta.get(src_key, hmi_map.meta.get(dst_key.lower(), 0.0))
    header["HGln_OBS"] = header["HGLN_OBS"]
    header["HGlt_OBS"] = header["HGLT_OBS"]
    return header


def make_wb_map(hmi_map, frame_index: int, wb_path: Path = SST_WB_ALIGNED_PATH):
    import sunpy.map

    with fits.open(wb_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data[frame_index, 0, 0], dtype=np.float32)
        corners = np.asarray(hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0][frame_index, :, :, :2], dtype=np.float64)
        src = hdul[0].header.copy()
    header = corners_to_tan_header(corners, data.shape, src, hmi_map)
    return sunpy.map.Map(data, header), (corners, src)


def estimate_nb_wavelengths(nb_header: fits.Header) -> np.ndarray:
    n_lambda = int(nb_header["NAXIS3"])
    wave_min = float(nb_header.get("WAVEMIN", nb_header.get("WAVELNTH", 395.0)))
    wave_max = float(nb_header.get("WAVEMAX", nb_header.get("WAVELNTH", 395.0)))
    if n_lambda <= 1:
        return np.array([wave_min], dtype=np.float64)
    return np.linspace(wave_min, wave_max, n_lambda, dtype=np.float64)


def load_nb_wing_image(nb_path: Path, scan_index: int) -> tuple[np.ndarray, fits.Header, tuple[int, int]]:
    if file_is_truncated_primary_hdu(nb_path):
        raise RuntimeError(
            f"{nb_path} is truncated before the end of the primary cube. "
            "A complete NB cube is needed for wing-image plotting."
        )
    with fits.open(nb_path, memmap=True) as hdul:
        src = hdul[0].header.copy()
        try:
            blue = np.asarray(hdul[0].data[scan_index, 0, 0], dtype=np.float32)
            red = np.asarray(hdul[0].data[scan_index, 0, -1], dtype=np.float32)
        except Exception as exc:  # pragma: no cover - runtime path for truncated data
            raise RuntimeError(
                f"Could not read the requested NB scan from {nb_path}. "
                "The local file appears truncated; a complete NB cube is needed for plotting."
            ) from exc
    return 0.5 * (blue + red), src, (0, int(src["NAXIS3"]) - 1)


def make_nb_wing_map(hmi_map, scan_index: int, nb_path: Path = SST_NB_SOURCE_PATH, wb_path: Path = SST_WB_ALIGNED_PATH):
    import sunpy.map

    image, src, wing_indices = load_nb_wing_image(nb_path, scan_index)
    corners = read_wb_corners(wb_path)[scan_index]
    header = corners_to_tan_header(corners, image.shape, src, hmi_map)
    wavelengths = estimate_nb_wavelengths(src)
    header["WAVELNTH"] = float(0.5 * (wavelengths[wing_indices[0]] + wavelengths[wing_indices[1]]) * 10.0)
    return sunpy.map.Map(image, header), wing_indices


def normalize_for_display(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(image)
    if not finite.any():
        return np.zeros_like(image, dtype=np.float32)
    fill = np.nanmedian(image[finite])
    image = np.where(finite, image, fill)
    p1, p99 = np.percentile(image, [1, 99])
    image = np.clip(image, p1, p99)
    if p99 <= p1:
        return np.zeros_like(image, dtype=np.float32)
    return (image - p1) / (p99 - p1)


def sst_world_to_pixel(corners: np.ndarray, shape: tuple[int, int], lon_arcsec: np.ndarray, lat_arcsec: np.ndarray):
    ny, nx = shape
    x_min = 0.5 * (corners[0, 0, 0] + corners[1, 0, 0])
    x_max = 0.5 * (corners[0, 1, 0] + corners[1, 1, 0])
    y_min = 0.5 * (corners[0, 0, 1] + corners[0, 1, 1])
    y_max = 0.5 * (corners[1, 0, 1] + corners[1, 1, 1])
    x = (lon_arcsec - x_min) / (x_max - x_min) * (nx - 1)
    y = (lat_arcsec - y_min) / (y_max - y_min) * (ny - 1)
    return x, y


def project_sst_to_iris(sst_image: np.ndarray, sst_corners: np.ndarray, iris_header: fits.Header):
    from scipy.ndimage import map_coordinates

    ny = iris_header["NAXIS2"]
    nx = iris_header["NAXIS1"]
    y, x = np.indices((ny, nx), dtype=np.float64)
    lon = (x + 1 - iris_header["CRPIX1"]) * iris_header["CDELT1"] + iris_header["CRVAL1"]
    lat = (y + 1 - iris_header["CRPIX2"]) * iris_header["CDELT2"] + iris_header["CRVAL2"]
    sst_x, sst_y = sst_world_to_pixel(sst_corners, sst_image.shape, lon, lat)
    sampled = map_coordinates(
        sst_image,
        [sst_y, sst_x],
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    return sampled.astype(np.float32)


def active_region_bounds(image: np.ndarray, margin: int = 95):
    valid = np.isfinite(image)
    filled = np.where(valid, image, np.nanmedian(image[valid]))
    dark = filled.copy()
    points = []
    for _ in range(6):
        idx = np.unravel_index(np.argmin(dark), dark.shape)
        points.append(idx)
        y, x = idx
        y0 = max(0, y - 60)
        y1 = min(dark.shape[0], y + 61)
        x0 = max(0, x - 60)
        x1 = min(dark.shape[1], x + 61)
        dark[y0:y1, x0:x1] = np.nanmax(dark)
    ys = [p[0] for p in points]
    xs = [p[1] for p in points]
    row_min = max(0, min(ys) - margin)
    row_max = min(image.shape[0] - 1, max(ys) + margin)
    col_min = max(0, min(xs) - margin)
    col_max = min(image.shape[1] - 1, max(xs) + margin)
    return row_min, row_max, col_min, col_max


def iris_extent_arcsec(header: fits.Header, row_min: int, row_max: int, col_min: int, col_max: int) -> list[float]:
    x0 = (col_min + 1 - header["CRPIX1"]) * header["CDELT1"] + header["CRVAL1"]
    x1 = (col_max + 1 - header["CRPIX1"]) * header["CDELT1"] + header["CRVAL1"]
    y0 = (row_min + 1 - header["CRPIX2"]) * header["CDELT2"] + header["CRVAL2"]
    y1 = (row_max + 1 - header["CRPIX2"]) * header["CDELT2"] + header["CRVAL2"]
    return [float(x0), float(x1), float(y0), float(y1)]


def build_nb_wcs_tab_from_aligned_wb(
    wb_path: Path,
    nb_header: fits.Header,
) -> fits.BinTableHDU:
    wb_tab = read_wb_tab(wb_path)
    n_scan = wb_tab.shape[0]
    wavelengths = estimate_nb_wavelengths(nb_header)
    n_lambda = len(wavelengths)
    coords = np.empty((n_scan, n_lambda, 2, 2, 4), dtype=np.float64)
    coords[..., 0] = wb_tab[:, None, :, :, 0]
    coords[..., 1] = wb_tab[:, None, :, :, 1]
    coords[..., 2] = wavelengths[None, :, None, None]
    coords[..., 3] = wb_tab[:, None, :, :, 2]
    hpln_index = np.array([[1.0, 2.0]], dtype=np.float32)
    hplt_index = np.array([[1.0, 2.0]], dtype=np.float32)
    columns = fits.ColDefs(
        [
            fits.Column(
                name="HPLN+HPLT+WAVE+TIME",
                format=f"{coords.size}D",
                dim=f"(4,2,2,{n_lambda},{n_scan})",
                array=[coords],
            ),
            fits.Column(name="HPLN-INDEX", format="2E", array=hpln_index),
            fits.Column(name="HPLT-INDEX", format="2E", array=hplt_index),
        ]
    )
    return fits.BinTableHDU.from_columns(columns, name="WCS-TAB")


def update_nb_cube_wcs_from_wb(
    nb_path: Path,
    wb_aligned_path: Path = SST_WB_ALIGNED_PATH,
    output_path: Path = SST_NB_ALIGNED_PATH,
) -> Path:
    nb_path = Path(nb_path)
    wb_aligned_path = Path(wb_aligned_path)
    output_path = Path(output_path)
    if file_is_truncated_primary_hdu(nb_path):
        raise RuntimeError(
            f"{nb_path} is truncated before the end of the primary cube. "
            "Refusing to copy or rewrite it. Please provide a complete NB cube."
        )
    try:
        with fits.open(nb_path, memmap=True) as hdul:
            nb_header = hdul[0].header.copy()
    except Exception as exc:  # pragma: no cover - depends on local NB file state
        raise RuntimeError(
            f"Could not open {nb_path} as a writable FITS cube. "
            "This local NB file appears truncated; use a complete cube to write an aligned NB product."
        ) from exc

    if output_path.exists():
        output_path.unlink()
    copy_file_with_progress(nb_path, output_path)
    wcs_hdu = build_nb_wcs_tab_from_aligned_wb(wb_aligned_path, nb_header)
    print("Preparing NB WCS-TAB update...", file=sys.stderr, flush=True)
    with fits.open(output_path, mode="update", memmap=True) as hdul:
        hdr = hdul[0].header
        wb_hdr = fits.getheader(wb_aligned_path, 0)
        hdr["ALNWX"] = float(wb_hdr.get("ALNWX", 0.0))
        hdr["ALNWY"] = float(wb_hdr.get("ALNWY", 0.0))
        hdr["ALNMANX"] = float(wb_hdr.get("ALNMANX", 0.0))
        hdr["ALNMANY"] = float(wb_hdr.get("ALNMANY", 0.0))
        hdr["ALNWB"] = str(wb_aligned_path.name)
        if "WCS-TAB" in hdul:
            idx = hdul.index_of("WCS-TAB")
            hdul[idx] = wcs_hdu
        else:
            hdul.append(wcs_hdu)
        print("Writing WCS-TAB to aligned NB cube...", file=sys.stderr, flush=True)
        hdul.flush()
    print("NB cube update complete.", file=sys.stderr, flush=True)
    return output_path
