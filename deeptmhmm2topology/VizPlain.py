import numpy as np
import plotly.graph_objects as go
import ast




def VizPlain_plot(TopologyDataframe, R,
                  range=(0,1310),
                  y0=None, 
                  y1= None,
                  display_circle=False, 
                  show_path=True,
                  add_NtermAnnotation=True,
                  add_CtermAnnotation=True):
    ##################################################
    # basic lines and dots
    ##################################################
    sub_df = TopologyDataframe.iloc[range[0]:range[-1]]
    centers_list_tmp = sub_df["center"].tolist()
    centers_list = [ast.literal_eval(item) for item in centers_list_tmp]
    CENTERS_arr = np.asarray(centers_list)
    xs = CENTERS_arr[:, 0]
    ys = CENTERS_arr[:, 1]


    fig = go.Figure()
    ## membrane position
    if y0 and y1:
        fig.add_hrect(y0=y0, y1=y1, line_width=0, fillcolor="tan", opacity=0.5, layer="below")

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
                mode="lines+markers+text",
                line=dict(shape="spline", width=1),
                marker=dict(size=10),
                hovertext=sub_df["position"].astype(str) + sub_df["AA"],
                hoverinfo="text"
            )
        )

    if add_NtermAnnotation:
        NtermAA_x, NtermAA_y = centers_list[0][0], centers_list[0][-1]
        fig.add_annotation(x=NtermAA_x, y=NtermAA_y,
                        text="N<sub>2</sub>H-",
                        showarrow=False,
                        xshift=-30)
    if add_CtermAnnotation:
        CtermAA_x, CtermAA_y = centers_list[-1][0], centers_list[-1][-1]
        fig.add_annotation(x=CtermAA_x, y=CtermAA_y,
                        text="-COOH",
                        showarrow=False,
                        xshift=30)


    x_min = xs.min() - 2 * R
    x_max = xs.max() + 2 * R
    y_min = ys.min() - 3 * R
    y_max = ys.max() + 3 * R



    fig.update_xaxes(range=[x_min, x_max], scaleanchor="y", visible=False)
    fig.update_yaxes(range=[y_min, y_max], visible=False)
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", 
        plot_bgcolor="rgba(0,0,0,0)", 
        autosize=True,
        margin=dict(l=20, r=20, t=0, b=0),
        showlegend=False
            )

    fig.write_html("protein_before.html")