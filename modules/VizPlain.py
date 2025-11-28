import numpy as np
import plotly.graph_objects as go

def VizPlain_plot(centers_list, R, y0_y1=None, display_circle=False, show_path=True):
    CENTERS_arr = np.asarray(centers_list)
    xs = CENTERS_arr[:, 0]
    ys = CENTERS_arr[:, 1]
    show_path = True

    fig = go.Figure()
    if y0_y1:
        fig.add_hrect(y0=y0_y1[0], y1=y0_y1[1], line_width=0, fillcolor="tan", opacity=0.5, layer="below")

    if display_circle:
        for (cx, cy) in CENTERS_arr:
            fig.add_shape(
                type="circle",
                xref="x", yref="y",
                x0=cx - R, x1=cx + R,
                y0=cy - R, y1=cy + R,
                line=dict(width=2),
                fillcolor="rgba(0,0,0,0)",
            )

    if show_path:
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines+markers",
                line=dict(shape="spline", width=5),
                marker=dict(size=10),
                name="centers",
            )
        )

    # --- Compute range ---
    x_min = xs.min() - 2 * R
    x_max = xs.max() + 2 * R
    y_min = ys.min() - 3 * R
    y_max = ys.max() + 3 * R



    # --- final layout ---
    fig.update_xaxes(range=[x_min, x_max], scaleanchor="y", visible=False)
    fig.update_yaxes(range=[y_min, y_max], visible=False)
    fig.update_layout(width=1200, height=900, showlegend=False)

    fig.write_html("protein.html")