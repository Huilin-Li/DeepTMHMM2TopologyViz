import numpy as np
import operator
from . import tools


class TopologyCenters:
    def __init__(self, R, away, dff3_df):
        self.start_center = (0, 0)
        self.away = int(away)
        self.R = float(R)
        self.df = dff3_df
        self.centers = None
        self.TMCircleCenters_DICT = None

    def generate(self):
        df = self.df
        Nterm_IMO = df.iloc[0]["IMO"]
        if Nterm_IMO == "outside":
            self.OutsideNterm()
        else:
            self.InsideNterm()
        return self

    def OutsideNterm(self):
        df = self.df
        for idx, row in df.iterrows():
            length = row["length"]
            if idx == len(df)-1:
                self.addCterm(length)
                break
            if idx == 0:
                # Nterm
                self.addNterm(length)
                print("idx=", idx, len(self.centers))
            elif idx%4 == 1:
                # download TM
                self.addDownwardingTMCenters(idx)
                print("idx=", idx, len(self.centers))
            elif idx%4 == 2:
                self.addIntracellularNotTMCenters(length) 
                print("idx=", idx, len(self.centers))
            elif idx%4 == 3:
                self.addUpwardingTMCenters(idx)
                print("idx=", idx, len(self.centers))
            else:
                self.addExtracellularNotTMCenters(length)
                print("idx=", idx, len(self.centers))
        return self
    
    def InsideNterm(self):
        df = self.df
        for idx, row in df.iterrows():
            length = row["length"]
            if idx == len(df)-1:
                self.addCterm(length)
                break
            if idx == 0:
                # Nterm
                self.addNterm(length)
                print("idx=", idx, len(self.centers))
            elif idx%4 == 1:
                # download TM
                self.addUpwardingTMCenters(idx) 
                print("idx=", idx, len(self.centers))
            elif idx%4 == 2:
                self.addExtracellularNotTMCenters(length) 
            elif idx%4 == 3:
                self.addDownwardingTMCenters(idx)
            else:
                self.addIntracellularNotTMCenters(length)
        return self

    def genTMCircleRelativeCenters(self, membraneThickness):
        TMUnits_idx_centers = tools.GenerateTMCircleRelativeCenters(df=self.df, membraneThickness=membraneThickness, R=self.R)
        self.TMCircleCenters_DICT = TMUnits_idx_centers
        return self
    
    def genMembraneYUpBottom(self):
        pass

    def addNterm(self, length):
        df = self.df
        CENTERs_bridge_list = [self.start_center]
        Nterm_IMO = df.iloc[0]["IMO"]
        Nterm_centers_bridge = tools.AddNterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Nterm_IMO)
        print("Nterm_centers_bridge", len(Nterm_centers_bridge))
        CENTERs_bridge_list += Nterm_centers_bridge[1:]
        self.centers = CENTERs_bridge_list
        return self
    
    def addCterm(self, length):
        df = self.df
        CENTERs_bridge_list = self.centers
        Cterm_IMO = df.iloc[-1]["IMO"]
        Cterm_centers_bridge = tools.AddCterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Cterm_IMO)
        print("Cterm_centers_bridge", len(Cterm_centers_bridge))
        CENTERs_bridge_list += Cterm_centers_bridge[1:]
        self.centers = CENTERs_bridge_list
        return self

    def addUpwardingTMCenters(self, idx):
        # upwarding
        CENTERs_bridge_list = self.centers
        upTMCenters_relative = self.TMCircleCenters_DICT[idx]
        df = self.df
        pre_IMO = df.iloc[idx-1]["IMO"]
        if pre_IMO == "outside":
            upTMCenters_relative = upTMCenters_relative[::-1]

        upTMCenters = tools.MoveCoords(coords=upTMCenters_relative, new_start=CENTERs_bridge_list[-1])
        CENTERs_bridge_list += upTMCenters
        self.centers = CENTERs_bridge_list
        return self
    
    def addDownwardingTMCenters(self, idx):
        # downwarding
        CENTERs_bridge_list = self.centers
        downTMCenters_relative = self.TMCircleCenters_DICT[idx]
        df = self.df
        pre_IMO = df.iloc[idx-1]["IMO"]
        if pre_IMO == "outside":
            downTMCenters_relative = downTMCenters_relative[::-1]
        
        downTMCenters = tools.MoveCoords(coords=downTMCenters_relative, new_start=CENTERs_bridge_list[-1])
        print("downTMCenters", len(downTMCenters), len(set(downTMCenters)))
        CENTERs_bridge_list += downTMCenters
        self.centers = CENTERs_bridge_list
        return self
    
    def addExtracellularNotTMCenters(self, length):
        CENTERs_bridge_list = self.centers
        ExtracellularNotTMCenters = tools.AddExtracellularNotTMCenters(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R)
        CENTERs_bridge_list += ExtracellularNotTMCenters
        self.centers = CENTERs_bridge_list
        return self

    def addIntracellularNotTMCenters(self, length):
        CENTERs_bridge_list = self.centers
        IntracellularNotTMCenters = tools.AddIntracellularNotTMCenters(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R)
        CENTERs_bridge_list += IntracellularNotTMCenters
        self.centers = CENTERs_bridge_list
        return self
