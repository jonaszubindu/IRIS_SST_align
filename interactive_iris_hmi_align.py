from __future__ import annotations

import base64
import io
import json
import os
import signal
import webbrowser

os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/sunpy")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import dash
import matplotlib
import numpy as np
import plotly.graph_objects as go
from astropy.io import fits
from dash import Input, Output, State, dcc, html

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from alignment_common import (
    IRIS_ALIGNED_PATH,
    REPORT_PATH,
    WORKDIR,
    ensure_initial_alignment,
    load_hmi_map,
    load_report,
    make_iris_map,
    normalize_for_display,
    project_map_to_extent,
)


APP_PORT = 8056


def load_session() -> dict:
    if not REPORT_PATH.exists():
        ensure_initial_alignment()
    report = load_report()
    hmi_map = load_hmi_map()
    with fits.open(IRIS_ALIGNED_PATH, memmap=False) as hdul:
        iris_index = int(report["iris"]["frame_index"])
        iris_image = np.asarray(hdul[0].data[iris_index], dtype=np.float32)
        iris_header = hdul[0].header.copy()
        base_total_dx = float(report["iris"]["world_shift_arcsec"][0])
        base_total_dy = float(report["iris"]["world_shift_arcsec"][1])

    base_crval1 = float(iris_header["CRVAL1"])
    base_crval2 = float(iris_header["CRVAL2"])
    x0 = float(iris_header["CRVAL1"] - (iris_header["CRPIX1"] - 1) * iris_header["CDELT1"])
    x1 = float(iris_header["CRVAL1"] + (iris_header["NAXIS1"] - iris_header["CRPIX1"]) * iris_header["CDELT1"])
    y0 = float(iris_header["CRVAL2"] - (iris_header["CRPIX2"] - 1) * iris_header["CDELT2"])
    y1 = float(iris_header["CRVAL2"] + (iris_header["NAXIS2"] - iris_header["CRPIX2"]) * iris_header["CDELT2"])
    pad = 18.0
    extent = [min(x0, x1) - pad, max(x0, x1) + pad, min(y0, y1) - pad, max(y0, y1) + pad]
    return {
        "hmi_map": hmi_map,
        "iris_index": iris_index,
        "iris_image": iris_image,
        "iris_header": iris_header,
        "extent": extent,
        "base_center_x": base_crval1,
        "base_center_y": base_crval2,
        "base_total_dx": base_total_dx,
        "base_total_dy": base_total_dy,
    }


SESSION = load_session()


def render_overlay_png(center_x: float, center_y: float) -> str:
    import sunpy.map

    iris_header = SESSION["iris_header"].copy()
    iris_header["CRVAL1"] = float(center_x)
    iris_header["CRVAL2"] = float(center_y)
    iris_map = sunpy.map.Map(SESSION["iris_image"], iris_header)
    hmi_map = SESSION["hmi_map"]
    extent = SESSION["extent"]
    shape = SESSION["iris_image"].shape
    hmi_image = normalize_for_display(project_map_to_extent(hmi_map, extent, shape))
    iris_image = normalize_for_display(project_map_to_extent(iris_map, extent, shape))

    fig = plt.figure(figsize=(8.8, 8.8), dpi=160)
    ax = fig.add_subplot(111)
    ax.imshow(hmi_image, extent=extent, origin="lower", cmap="binary_r", alpha=0.72, aspect="equal")
    ax.imshow(iris_image, extent=extent, origin="lower", cmap="gray", alpha=0.38, aspect="equal")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("Helioprojective Longitude (Solar-X)")
    ax.set_ylabel("Helioprojective Latitude (Solar-Y)")
    ax.set_title("Click to set the new IRIS center over HMI. Press 's' or click Save.")
    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=160)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def build_figure(center_x: float, center_y: float) -> go.Figure:
    view_extent = SESSION["extent"]
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=np.zeros((60, 60), dtype=np.float32),
            x=np.linspace(view_extent[0], view_extent[1], 60),
            y=np.linspace(view_extent[2], view_extent[3], 60),
            colorscale=[[0.0, "rgba(0,0,0,0)"], [1.0, "rgba(0,0,0,0)"]],
            showscale=False,
            opacity=0.02,
            hovertemplate="Click target<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>",
            zmin=0.0,
            zmax=1.0,
        )
    )
    fig.update_layout(
        clickmode="event",
        xaxis={"title": "Solar X [arcsec]", "scaleanchor": "y", "range": [view_extent[0], view_extent[1]]},
        yaxis={"title": "Solar Y [arcsec]", "range": [view_extent[2], view_extent[3]]},
        margin={"l": 50, "r": 20, "t": 50, "b": 50},
        template="plotly_white",
        uirevision="manual-iris-hmi",
        images=[
            {
                "source": render_overlay_png(center_x, center_y),
                "xref": "x",
                "yref": "y",
                "x": view_extent[0],
                "y": view_extent[3],
                "sizex": view_extent[1] - view_extent[0],
                "sizey": view_extent[3] - view_extent[2],
                "sizing": "stretch",
                "layer": "below",
                "opacity": 1.0,
            }
        ],
    )
    return fig


def save_solution(center_x: float, center_y: float):
    extra_dx = center_x - SESSION["base_center_x"]
    extra_dy = center_y - SESSION["base_center_y"]
    total_dx = SESSION["base_total_dx"] + extra_dx
    total_dy = SESSION["base_total_dy"] + extra_dy

    with fits.open(IRIS_ALIGNED_PATH, mode="update", memmap=False) as hdul:
        hdr = hdul[0].header
        hdr["CRVAL1"] = float(center_x)
        hdr["CRVAL2"] = float(center_y)
        hdr["ALNIRX"] = total_dx
        hdr["ALNIRY"] = total_dy
        hdr["ALNIRMX"] = extra_dx
        hdr["ALNIRMY"] = extra_dy
        hdul.flush()

    report = load_report()
    report["iris"]["world_shift_arcsec"] = [total_dx, total_dy]
    report["iris_manual_adjustment_arcsec"] = {"dx": extra_dx, "dy": extra_dy}
    REPORT_PATH.write_text(json.dumps(report, indent=2))
    return extra_dx, extra_dy, total_dx, total_dy


app = dash.Dash(__name__)
app.title = "IRIS Manual Align"
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
            "Click anywhere on the plot to place the IRIS center there. Press 's' or click Save when happy.",
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
        f"total IRIS shift: dx={total_dx:+.2f}, dy={total_dy:+.2f} arcsec"
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
        f"Current IRIS total dx={total_dx:+.2f}, dy={total_dy:+.2f} arcsec. "
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
    webbrowser.open(url)
    app.run_server(debug=False, port=APP_PORT)


if __name__ == "__main__":
    main()
