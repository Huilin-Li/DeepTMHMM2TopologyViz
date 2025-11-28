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
        self.membraney0 = None
        self.membraney1 = None
        self.Height = None
        self.TopoCenters = None

    def generate(self):
        df = self.df
        Nterm_IMO = df.iloc[0]["IMO"]
        if Nterm_IMO == "outside":
            self.OutsideNterm()
        else:
            self.InsideNterm()
        return self
    
    def remove_duplicates(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    
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
                self.membraney1 = self.centers[-1][-1]
                self.membraney0 = self.membraney1 - self.Height 
            elif idx%4 == 1:
                # download TM
                self.addDownwardingTMCenters(idx)
            elif idx%4 == 2:
                self.addIntracellularNotTMCenters(length) 
            elif idx%4 == 3:
                self.addUpwardingTMCenters(idx)
            else:
                self.addExtracellularNotTMCenters(length)

        CENTERS_arr = np.asarray(self.centers)
        CENTERS_arr = np.round(CENTERS_arr, 6)
        CENTERS_arr_list = list(map(tuple, CENTERS_arr.tolist()))
        self.TopoCenters = self.remove_duplicates(CENTERS_arr_list)
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
                self.membraney0 = self.centers[-1][-1]
                self.membraney1 = self.membraney0 + self.Height 
            elif idx%4 == 1:
                # download TM
                self.addUpwardingTMCenters(idx) 
            elif idx%4 == 2:
                self.addExtracellularNotTMCenters(length) 
            elif idx%4 == 3:
                self.addDownwardingTMCenters(idx)
            else:
                self.addIntracellularNotTMCenters(length)
        CENTERS_arr = np.asarray(self.centers)
        CENTERS_arr = np.round(CENTERS_arr, 6)
        CENTERS_arr_list = list(map(tuple, CENTERS_arr.tolist()))
        self.TopoCenters = self.remove_duplicates(CENTERS_arr_list)
        return self

    def genTMCircleRelativeCenters(self, membraneThickness):
        TMUnits_idx_centers = tools.GenerateTMCircleRelativeCenters(df=self.df, membraneThickness=membraneThickness, R=self.R)
        self.TMCircleCenters_DICT = TMUnits_idx_centers
        TMUnits_idx_centers = self.TMCircleCenters_DICT
        Ybottoms = []
        Yups = []
        for k,v in TMUnits_idx_centers.items():
            Ybottoms.append(v[0][-1])
            Yups.append(v[-1][-1])
        bottom = sum(Ybottoms) / float(len(Ybottoms))
        up = sum(Yups) / float(len(Yups))
        self.Height = up - bottom
        return self
    
    def genMembraneYUpBottom(self):
        pass

    def addNterm(self, length):
        df = self.df
        CENTERs_bridge_list = [self.start_center]
        Nterm_IMO = df.iloc[0]["IMO"]
        Nterm_centers_bridge = tools.AddNterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Nterm_IMO)
        CENTERs_bridge_list += Nterm_centers_bridge[1:]
        self.centers = CENTERs_bridge_list
        return self
    
    def addCterm(self, length):
        df = self.df
        CENTERs_bridge_list = self.centers
        Cterm_IMO = df.iloc[-1]["IMO"]
        Cterm_centers_bridge = tools.AddCterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Cterm_IMO)
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
