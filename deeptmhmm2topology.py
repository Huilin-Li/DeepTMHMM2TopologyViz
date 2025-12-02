import argparse
import os
from deeptmhmm2topology.DeepTMHMM2Topology import TopologyCenters
from deeptmhmm2topology.tools import Convert_dff32df, Convert_3line2df


def main():
    parser = argparse.ArgumentParser(description="DeepTMHMM 2 Topology (Plain version).")
    #########################################################
    # Input arguments
    #########################################################
    parser.add_argument("gff3")
    parser.add_argument("line3")
    parser.add_argument("-r", argument_default=0.2)
    parser.add_argument("-a","--away", argument_default=4)
    parser.add_argument("-maxb", argument_default=10)
    parser.add_argument("-memThk","--membraneThickness", argument_default=3) 
    #########################################################
    # Output arguments
    #########################################################
    parser.add_argument("-oc","--CentersOut", argument_default="TopologyCenters.csv")
    parser.add_argument("-om","--MembraneOut", argument_default="TopologyMembraney0y1.text")
    args = parser.parse_args()

    # load arguments
    dff3_df=Convert_dff32df(args.gff3)
    line3_df=Convert_3line2df(args.line3)
    R = args.r
    away = args.a
    max_b = args.maxb
    membraneThickness = args.memThk

    #########################################################
    # Get topology centers
    #########################################################
    Topology = (
        TopologyCenters(R=R, away=away, dff3_df=dff3_df, line3_df=line3_df, max_b=max_b)
        .genTMCircleRelativeCenters(membraneThickness=membraneThickness)
        .generate()
    )
    TopoDataFrame = Topology.TopoDataFrame
    TopoDataFrame.to_csv(args.o, encoding='utf-8', index=False)

    #########################################################
    # Output membrane y0 and y1
    #########################################################
    y0 = Topology.membraney0
    y1 = Topology.membraney1
    


if __name__ == "__main__":
    main()

