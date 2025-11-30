from deeptmhmm2topology.DeepTMHMM2Topology import TopologyCenters
from deeptmhmm2topology.tools import Convert_dff32df, Convert_3line2df
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import ast
from plotly.subplots import make_subplots

#########################################################
# Read results from DeepTMHMM
#########################################################
# deeptmhmm_gff3_fp = ".\examples\Q92508_TMRs.gff3"
# deeptmhmm_3line_fp = ".\examples\Q92508predicted_topologies.3line"
# dff3_df = Convert_dff32df(deeptmhmm_gff3_fp)
# line3_df =Convert_3line2df(deeptmhmm_3line_fp)

#########################################################
# Get topology centers
#########################################################
# Topology = (
#     TopologyCenters(R=0.2, away=4, dff3_df=dff3_df, line3_df=line3_df, max_b=15)
#     .genTMCircleRelativeCenters(membraneThickness=3)
#     .generate()
# )
# TopoDataFrame = Topology.TopoDataFrame
# TopoDataFrame.to_csv(f"./TopologyCenters.csv", encoding='utf-8', index=False)

#########################################################
# Get membrane y0 and y1 axis
#########################################################
# y0 = Topology.membraney0
# y1 = Topology.membraney1
# print(y0, y1) 
# 2.0 4.821227555116749


#########################################################
# Visualization (customized example)
#########################################################
centers_df = pd.read_csv(".\TopologyCenters.csv")
## add more imformation, in this case, mutation records
mut_df = pd.read_csv(".\mutations.csv")
df = pd.merge(centers_df, mut_df, left_on='position', right_on='pos', how='outer')
mutation_list = list(set(df.dropna(subset=['mutation'])["mutation"].tolist()))
# print(mutation_list) =4
# ['deletion', 'termination', 'frameshift', 'insertion']
colors = ["red", "green", "chartreuse", "darkorange"]
mutation_color_map = {
    'deletion': 'red',
    'termination': 'green',
    'frameshift': 'royalblue',
    'insertion': 'darkorange'
}
y0=2
y1=4.821227555116749
fig = make_subplots(rows=2, cols=1, vertical_spacing=0.02)
# my protein is too long, so I split it into 2 parts at around 1310

## subplot row1col1
sub11df = df.iloc[:1310]
centers11_list_tmp = sub11df["center"].tolist()
centers11_list = [ast.literal_eval(item) for item in centers11_list_tmp]
CENTERS_arr11 = np.asarray(centers11_list)
xs11 = CENTERS_arr11[:, 0]
ys11 = CENTERS_arr11[:, 1]
colors11 = [mutation_color_map.get(mut, 'lightgray') for mut in sub11df["mutation"].fillna('no_mutation')]
fig.add_trace(
    go.Scatter(
        x=xs11,
        y=ys11,
        mode="lines+markers+text",
        line=dict(shape="spline", width=1, color="lightgray"),
        marker=dict(size=10, color=colors11),
        hovertext=sub11df["position"].astype(str) + sub11df["AA"],
        hoverinfo="text",
        showlegend=False,
        legendgroup = '1'
    ), row=1, col=1
)
fig.add_hrect(y0=y0, y1=y1, line_width=0, fillcolor="tan", opacity=0.5, layer="below", row=1, col=1)
NtermAA_x, NtermAA_y = xs11[0], ys11[0]
fig.add_annotation(x=NtermAA_x, y=NtermAA_y, text="N<sub>2</sub>H-", showarrow=False, xshift=-30, row=1, col=1)

## subplot row1col1
sub21df = df.iloc[1310:]
centers21_list_tmp = sub21df["center"].tolist()
centers21_list = [ast.literal_eval(item) for item in centers21_list_tmp]
CENTERS_arr21 = np.asarray(centers21_list)
xs21 = CENTERS_arr21[:, 0]
ys21 = CENTERS_arr21[:, 1]
colors21 = [mutation_color_map.get(mut, 'lightgray') for mut in sub21df["mutation"].fillna('no_mutation')]
fig.add_trace(
    go.Scatter(
        x=xs21,
        y=ys21,
        mode="lines+markers+text",
        line=dict(shape="spline", width=1, color="lightgray"),
        marker=dict(size=10, color=colors21),
        hovertext=sub21df["position"].astype(str) + sub21df["AA"],
        hoverinfo="text",
        showlegend=False,
        legendgroup = '2'
    ), row=2, col=1
)
fig.add_hrect(y0=y0, y1=y1, line_width=0, fillcolor="tan", opacity=0.5, layer="below", row=2, col=1)
CtermAA_x, CtermAA_y = xs21[-1], ys21[-1]
fig.add_annotation(x=CtermAA_x, y=CtermAA_y, text="-COOH", showarrow=False, xshift=30, row=2, col=1)



# legend
for mut, color in mutation_color_map.items():
    fig.add_trace(
        go.Scatter(
            x=[None], y=[None], 
            mode='markers', 
            marker=dict(size=10, color=color),
            name=mut,
            showlegend=True,
            legendgroup = '2'
        ), row=2, col=1
    )
fig.update_layout(
    legend=dict(
    x=0,  # Position the legend at the left
    y=1,  # Position it at the top
    xanchor='left',  # Anchor the legend to the left
    yanchor='bottom',   # Anchor the legend to the top
    orientation='v', # Make the legend vertical
    ))
# layout
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)
fig.update_layout(
    paper_bgcolor="rgba(0,0,0,0)", 
    plot_bgcolor="rgba(0,0,0,0)", 
    autosize=True,
    margin=dict(l=20, r=20, t=0, b=0),
    showlegend=True
        )

# fig.write_html("protein_q92508.html")
fig.write_image("protein_q92508.svg")