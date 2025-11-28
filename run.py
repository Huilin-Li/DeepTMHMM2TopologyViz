from modules.VizPlain import VizPlain_plot
from modules.DeepTMHMM2Topology import TopologyCenters
from modules.tools import Convert_dff32df
import numpy as np





deeptmhmm_gff3_list = [".\examples\P1PV6_1TMRs.gff3", ".\examples\Q92508_TMRs.gff3", ".\examples\P7Z0A_1.gff3", ".\examples\P5EH6_1.gff3"]
deeptmhmm_gff3 = deeptmhmm_gff3_list[2]
dff3_df = Convert_dff32df(deeptmhmm_gff3)
print(dff3_df)


# Topology
Topology = (
    TopologyCenters(R=0.25, away=4, dff3_df=dff3_df)
    .genTMCircleRelativeCenters(membraneThickness=3)
    .generate()
)

CENTERS_LIST = Topology.TopoCenters
y0_y1 = (Topology.membraney0, Topology.membraney1)


print(len(CENTERS_LIST), len(set(CENTERS_LIST)))

VizPlain_plot(CENTERS_LIST, R=0.25, y0_y1=y0_y1)

