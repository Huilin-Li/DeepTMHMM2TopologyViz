from deeptmhmm2topology.VizPlain import VizPlain_plot
from deeptmhmm2topology.DeepTMHMM2Topology import TopologyCenters
from deeptmhmm2topology.tools import Convert_dff32df, Convert_3line2df
import pandas as pd




deeptmhmm_gff3_list = [".\examples\P1PV6_1TMRs.gff3", ".\examples\Q92508_TMRs.gff3", ".\examples\P7Z0A_1.gff3", ".\examples\P5EH6_1.gff3"]
Convert_3line2df_list = [".\examples\Q92508predicted_topologies.3line"]
deeptmhmm_gff3 = deeptmhmm_gff3_list[1]
dff3_df = Convert_dff32df(deeptmhmm_gff3)
deeptmhmm_3line_fp = Convert_3line2df_list[0]
line3_df =Convert_3line2df(deeptmhmm_3line_fp)

print(dff3_df)
print(line3_df)


# # Topology
Topology = (
    TopologyCenters(R=0.25, away=4, dff3_df=dff3_df, line3_df=line3_df, max_a=2, max_b=6)
    .genTMCircleRelativeCenters(membraneThickness=3)
    .generate()
)

TopoDataFrame = Topology.TopoDataFrame
print(TopoDataFrame)

y0_y1 = (Topology.membraney0, Topology.membraney1)




# TopoDataFrame.to_csv(f"./Topology_{}_{}", encoding='utf-8', index=False)
print(y0_y1)

# VizPlain_plot(CENTERS_LIST, R=0.25, y0_y1=y0_y1)

f"He said his name is {name}."