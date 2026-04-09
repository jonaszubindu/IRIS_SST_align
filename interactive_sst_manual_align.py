from __future__ import annotations

import json
import os
import signal

import dash
import numpy as np
import plotly.graph_objects as go
from astropy.io import fits
from astropy.wcs import WCS
from dash import Input, Output, State, dcc, html

from alignment_common import (
    IRIS_ALIGNED_PATH,
    MANUAL_SAVE_PATH,
    REPORT_PATH,
    SST_WB_ALIGNED_PATH,
    WORKDIR,
    ensure_initial_alignment,
    load_report,
    normalize_for_display,
    project_data_wcs_to_extent,
    project_sst_to_iris,
    square_extent_from_bounds,
    world_bounds_from_corners,
    world_bounds_from_header,
)

APP_PORT = 8054


def load_session() -> dict:
    if not REPORT_PATH.exists():
        ensure_initial_alignment()
    report = load_report()
    with fits.open(IRIS_ALIGNED_PATH, memmap=False) as hdul:
        iris_index = int(report["iris"]["frame_index"])
        iris_image = np.asarray(hdul[0].data[iris_index], dtype=np.float32)
        iris_header = hdul[0].header.copy()
    with fits.open(SST_WB_ALIGNED_PATH, memmap=False) as hdul:
        sst_index = int(report["sst_wb"]["frame_index"])
        sst_image = np.asarray(hdul[0].data[sst_index, 0, 0], dtype=np.float32)
        sst_corners = np.asarray(hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0][sst_index, :, :, :2], dtype=np.float64)
        hdr = hdul[0].header
        base_total_dx = float(hdr.get("ALNWX", 0.0))
        base_total_dy = float(hdr.get("ALNWY", 0.0))

    iris_bounds = world_bounds_from_header(iris_header)
    sst_bounds = world_bounds_from_corners(sst_corners)
    bounds = (
        min(iris_bounds[0], sst_bounds[0]),
        max(iris_bounds[1], sst_bounds[1]),
        min(iris_bounds[2], sst_bounds[2]),
        max(iris_bounds[3], sst_bounds[3]),
    )
    extent = square_extent_from_bounds(bounds, pad_arcsec=12.0)
    shape = (900, 900)
    iris_crop = project_data_wcs_to_extent(iris_image, WCS(iris_header).celestial, extent, shape)
    # Build a temporary TAN header-like grid for the existing SST projector.
    temp_header = fits.Header()
    temp_header["NAXIS1"] = shape[1]
    temp_header["NAXIS2"] = shape[0]
    temp_header["CRPIX1"] = 1.0
    temp_header["CRPIX2"] = 1.0
    temp_header["CRVAL1"] = extent[0]
    temp_header["CRVAL2"] = extent[2]
    temp_header["CDELT1"] = (extent[1] - extent[0]) / (shape[1] - 1)
    temp_header["CDELT2"] = (extent[3] - extent[2]) / (shape[0] - 1)
    sst_crop = project_sst_to_iris(sst_image, sst_corners, temp_header)
    return {
        "iris_crop": normalize_for_display(iris_crop),
        "sst_crop": normalize_for_display(sst_crop),
        "extent": extent,
        "width": extent[1] - extent[0],
        "height": extent[3] - extent[2],
        "base_center_x": 0.5 * (extent[0] + extent[1]),
        "base_center_y": 0.5 * (extent[2] + extent[3]),
        "base_total_dx": base_total_dx,
        "base_total_dy": base_total_dy,
    }


SESSION = load_session()


def heatmap_trace(image: np.ndarray, extent: list[float], colorscale: str, opacity: float):
    ny, nx = image.shape
    return go.Heatmap(
        z=np.flipud(image),
        x0=extent[0],
        dx=(extent[1] - extent[0]) / max(nx - 1, 1),
        y0=extent[3],
        dy=-(extent[3] - extent[2]) / max(ny - 1, 1),
        colorscale=colorscale,
        showscale=False,
        opacity=opacity,
        hovertemplate="x=%{x:.2f} arcsec<br>y=%{y:.2f} arcsec<extra></extra>",
    )


def build_figure(center_x: float, center_y: float) -> go.Figure:
    width = SESSION["width"]
    height = SESSION["height"]
    pad_x = max(12.0, 0.35 * width)
    pad_y = max(12.0, 0.35 * height)
    view_extent = [
        SESSION["extent"][0] - pad_x,
        SESSION["extent"][1] + pad_x,
        SESSION["extent"][2] - pad_y,
        SESSION["extent"][3] + pad_y,
    ]
    sst_extent = [
        center_x - 0.5 * width,
        center_x + 0.5 * width,
        center_y - 0.5 * height,
        center_y + 0.5 * height,
    ]
    fig = go.Figure()
    fig.add_trace(heatmap_trace(SESSION["iris_crop"], SESSION["extent"], "Greys_r", 0.65))
    fig.add_trace(heatmap_trace(SESSION["sst_crop"], sst_extent, "Greys", 0.65))
    fig.add_trace(
        go.Scatter(
            x=[view_extent[0], view_extent[1], view_extent[1], view_extent[0]],
            y=[view_extent[2], view_extent[2], view_extent[3], view_extent[3]],
            mode="markers",
            marker={"size": 18, "opacity": 0.0},
            hovertemplate="Click target<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        clickmode="event",
        xaxis={"title": "Solar X [arcsec]", "scaleanchor": "y", "range": [view_extent[0], view_extent[1]]},
        yaxis={"title": "Solar Y [arcsec]", "range": [view_extent[2], view_extent[3]]},
        margin={"l": 50, "r": 20, "t": 50, "b": 50},
        title="Click to set the new SST center over IRIS. Press 's' or click Save.",
        template="plotly_white",
        uirevision="manual-align",
    )
    return fig


def save_solution(center_x: float, center_y: float):
    extra_dx = center_x - SESSION["base_center_x"]
    extra_dy = center_y - SESSION["base_center_y"]
    total_dx = SESSION["base_total_dx"] + extra_dx
    total_dy = SESSION["base_total_dy"] + extra_dy

    with fits.open(SST_WB_ALIGNED_PATH, mode="update", memmap=False) as hdul:
        tab = np.asarray(hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0], dtype=np.float64)
        tab[..., 0] += extra_dx
        tab[..., 1] += extra_dy
        hdul["WCS-TAB"].data["HPLN+HPLT+TIME"][0] = tab
        hdr = hdul[0].header
        hdr["ALNWX"] = total_dx
        hdr["ALNWY"] = total_dy
        hdr["ALNMANX"] = extra_dx
        hdr["ALNMANY"] = extra_dy
        hdul.flush()

    report = load_report()
    report["manual_adjustment_arcsec"] = {"dx": extra_dx, "dy": extra_dy}
    report["current_total_sst_shift_arcsec"] = {"dx": total_dx, "dy": total_dy}
    REPORT_PATH.write_text(json.dumps(report, indent=2))

    MANUAL_SAVE_PATH.write_text(
        json.dumps(
            {
                "extra_shift_arcsec": {"dx": extra_dx, "dy": extra_dy},
                "total_sst_shift_arcsec": {"dx": total_dx, "dy": total_dy},
                "new_center_arcsec": {"x": center_x, "y": center_y},
                "notes": "Saved from plotly interactive click-to-center aligner.",
            },
            indent=2,
        )
    )
    return extra_dx, extra_dy, total_dx, total_dy


app = dash.Dash(__name__)
app.title = "SST Manual Align"
app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <script>
        document.addEventListener('keydown', function(e) {
            if ((e.key || '').toLowerCase() === 's') {
                const btn = document.getElementById('save-button');
                if (btn) { btn.click(); }
            }
        });
        </script>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

app.layout = html.Div(
    [
        dcc.Store(id="center-store", data={"x": SESSION["base_center_x"], "y": SESSION["base_center_y"]}),
        dcc.Store(id="app-state-store", data={"saved": False, "closing": False}),
        dcc.Interval(id="shutdown-interval", interval=1200, n_intervals=0, disabled=True),
        html.Div(
            "Click anywhere on the plot to place the SST center there. Press 's' or click Save when happy.",
            style={"marginBottom": "10px"},
        ),
        dcc.Graph(
            id="align-graph",
            figure=build_figure(SESSION["base_center_x"], SESSION["base_center_y"]),
            style={"height": "92vh"},
            config={"displayModeBar": True, "scrollZoom": True},
        ),
        html.Div(
            [
                html.Button("Reset", id="reset-button", n_clicks=0),
                dcc.Input(
                    id="x-input",
                    type="number",
                    value=SESSION["base_center_x"],
                    debounce=True,
                    placeholder="Center X [arcsec]",
                    style={"marginLeft": "12px", "width": "150px"},
                ),
                dcc.Input(
                    id="y-input",
                    type="number",
                    value=SESSION["base_center_y"],
                    debounce=True,
                    placeholder="Center Y [arcsec]",
                    style={"marginLeft": "8px", "width": "150px"},
                ),
                html.Button("Shift", id="shift-button", n_clicks=0, style={"marginLeft": "8px"}),
                html.Button("Save", id="save-button", n_clicks=0, style={"marginLeft": "8px"}),
                html.Button("Quit And Close", id="quit-button", n_clicks=0, disabled=False, style={"marginLeft": "8px"}),
                html.Span(id="coord-status", style={"marginLeft": "14px", "whiteSpace": "pre-line"}),
            ],
            style={"marginTop": "10px"},
        ),
        html.Div(
            id="save-status",
            style={
                "marginTop": "10px",
                "fontWeight": "bold",
                "padding": "10px 12px",
                "borderRadius": "8px",
                "display": "none",
            },
        ),
    ],
    style={"maxWidth": "1600px", "margin": "0 auto", "padding": "16px"},
)


@app.callback(
    Output("center-store", "data"),
    Input("align-graph", "clickData"),
    Input("reset-button", "n_clicks"),
    Input("shift-button", "n_clicks"),
    State("center-store", "data"),
    State("x-input", "value"),
    State("y-input", "value"),
    prevent_initial_call=True,
)
def update_center(click_data, _reset_clicks, _shift_clicks, current, x_value, y_value):
    ctx = dash.callback_context
    trigger = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else ""
    if trigger == "reset-button":
        return {"x": SESSION["base_center_x"], "y": SESSION["base_center_y"]}
    if trigger == "shift-button":
        if x_value is None or y_value is None:
            return current
        return {"x": float(x_value), "y": float(y_value)}
    if click_data and click_data.get("points"):
        point = click_data["points"][0]
        return {"x": float(point["x"]), "y": float(point["y"])}
    return current


@app.callback(
    Output("align-graph", "figure"),
    Output("coord-status", "children"),
    Output("x-input", "value"),
    Output("y-input", "value"),
    Input("center-store", "data"),
)
def redraw(center):
    x = float(center["x"])
    y = float(center["y"])
    extra_dx = x - SESSION["base_center_x"]
    extra_dy = y - SESSION["base_center_y"]
    total_dx = SESSION["base_total_dx"] + extra_dx
    total_dy = SESSION["base_total_dy"] + extra_dy
    text = (
        f"clicked center: x={x:+.2f}, y={y:+.2f} arcsec    "
        f"extra shift: dx={extra_dx:+.2f}, dy={extra_dy:+.2f} arcsec    "
        f"total shift: dx={total_dx:+.2f}, dy={total_dy:+.2f} arcsec"
    )
    return build_figure(x, y), text, x, y


@app.callback(
    Output("save-status", "children"),
    Output("app-state-store", "data"),
    Input("save-button", "n_clicks"),
    State("center-store", "data"),
    running=[
        (Output("save-button", "disabled"), True, False),
        (Output("save-button", "children"), "Saving...", "Save"),
        (Output("quit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
def save_current_solution(_n_clicks, center):
    extra_dx, extra_dy, total_dx, total_dy = save_solution(float(center["x"]), float(center["y"]))
    return (
        f"Save complete. extra dx={extra_dx:+.2f}, extra dy={extra_dy:+.2f} arcsec. "
        f"Current total dx={total_dx:+.2f}, dy={total_dy:+.2f} arcsec. "
        f"Plots were not regenerated yet. Close the app first, then we can refresh them separately."
    ), {"saved": True, "closing": False}


@app.callback(
    Output("quit-button", "disabled"),
    Output("quit-button", "children"),
    Input("app-state-store", "data"),
)
def update_quit_button(state):
    state = state or {}
    if state.get("closing"):
        return True, "Closing..."
    return False, "Quit And Close"


@app.callback(
    Output("save-status", "style"),
    Input("app-state-store", "data"),
)
def save_status_style(state):
    state = state or {}
    if state.get("closing"):
        return {
            "marginTop": "10px",
            "fontWeight": "bold",
            "color": "#0b4f8a",
            "backgroundColor": "#e8f1fb",
            "padding": "10px 12px",
            "borderRadius": "8px",
            "display": "block",
        }
    if state.get("saved"):
        return {
            "marginTop": "10px",
            "fontWeight": "bold",
            "color": "#176b2c",
            "backgroundColor": "#eaf7ed",
            "padding": "10px 12px",
            "borderRadius": "8px",
            "display": "block",
        }
    return {
        "marginTop": "10px",
        "fontWeight": "bold",
        "padding": "10px 12px",
        "borderRadius": "8px",
        "display": "none",
    }


@app.callback(
    Output("save-status", "children", allow_duplicate=True),
    Output("app-state-store", "data", allow_duplicate=True),
    Output("shutdown-interval", "disabled"),
    Input("quit-button", "n_clicks"),
    State("app-state-store", "data"),
    prevent_initial_call=True,
)
def begin_quit(n_clicks, state):
    state = state or {}
    if n_clicks:
        return (
            "The app is closing now. You can close this browser tab.",
            {"saved": bool(state.get("saved")), "closing": True},
            False,
        )
    return dash.no_update, dash.no_update, True


@app.callback(
    Output("quit-button", "n_clicks", allow_duplicate=True),
    Input("shutdown-interval", "n_intervals"),
    State("app-state-store", "data"),
    prevent_initial_call=True,
)
def quit_app(_n_intervals, state):
    state = state or {}
    if state.get("closing"):
        os.kill(os.getpid(), signal.SIGTERM)
    return 0


def main() -> None:
    url = f"http://127.0.0.1:{APP_PORT}"
    print(f"Open {url}")
    app.run(debug=False, port=APP_PORT)


if __name__ == "__main__":
    main()
