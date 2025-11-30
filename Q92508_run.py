from deeptmhmm2topology.VizPlain import VizPlain_plot
from deeptmhmm2topology.DeepTMHMM2Topology import TopologyCenters
from deeptmhmm2topology.tools import Convert_dff32df, Convert_3line2df
import pandas as pd

#########################################################
# Read results from DeepTMHMM
#########################################################
deeptmhmm_gff3_fp = ".\examples\Q92508_TMRs.gff3"
deeptmhmm_3line_fp = ".\examples\Q92508predicted_topologies.3line"
dff3_df = Convert_dff32df(deeptmhmm_gff3_fp)
line3_df =Convert_3line2df(deeptmhmm_3line_fp)

#########################################################
# Get topology centers
#########################################################
Topology = (
    TopologyCenters(R=0.2, away=4, dff3_df=dff3_df, line3_df=line3_df, max_b=15)
    .genTMCircleRelativeCenters(membraneThickness=3)
    .generate()
)
TopoDataFrame = Topology.TopoDataFrame
TopoDataFrame.to_csv(f"./TopologyCenters.csv", encoding='utf-8', index=False)

#########################################################
# Get membrane y0 and y1 axis
#########################################################
y0 = Topology.membraney0
y1 = Topology.membraney1
print(y0, y1) 
# 2.0 4.821227555116749


