import numpy as np
import operator
from . import tools

class ProteinTopologyCenters:
    def __init__(self, R, away, dff3_df):
        self.start_center = (0, 0)
        self.away = int(away)
        self.R = float(R)
        self.df = dff3_df

        self.centers = None
        self.TMCircleCenters_DICT = None


    def genTMCircleRelativeCenters(self, membraneThickness):
        TMUnits_idx_centers = tools.GenerateTMCircleRelativeCenters(df=self.df, membraneThickness=membraneThickness, R=self.R)
        self.TMCircleCenters_DICT = TMUnits_idx_centers
        return self
    
    def genNtermCenters(self, length, IMO):
        if IMO == "inside":
            # downwarding
            upward = False
            pre_center = self.start_center
            restCircles = length - self.away # one side
            # get S curve
            centers_list = [pre_center]
            restCircles_comb = tools.CirclesCombinations(restCircles)
            if isinstance(restCircles_comb, list):
                Scurve_centers = tools._gen_horizontalS(pre_center=centers_list[-1], R=self.R, restCircles_comb_is_list=restCircles_comb, upward=upward)
                centers_list += Scurve_centers[1:]
                # one side away
                pre_center = centers_list[-1]
                next_CircleCenters = tools._gen_straight_nCircleCenters(pre_center=pre_center, n=self.away, R=self.R, upward=True)
                centers_list += next_CircleCenters[1:]
                self.centers = centers_list
                return self
                    
            else:
                # short Nterm
                # downwarding straight
                pre_center = centers_list[-1]
                next_CircleCenters = tools._gen_straight_nCircleCenters(pre_center=pre_center, n=restCircles_comb, R=self.R, upward=True)
                centers_list += next_CircleCenters[1:]
                self.centers = centers_list
                return self
        else:
            # outside
            upward = True
            pre_center = self.start_center
            restCircles = length - self.away # one side
            # get S curve
            centers_list = [pre_center]
            restCircles_comb = tools.CirclesCombinations(restCircles)
            if isinstance(restCircles_comb, list):
                Scurve_centers = tools._gen_horizontalS(pre_center=centers_list[-1], R=self.R, restCircles_comb_is_list=restCircles_comb, upward=upward)
                centers_list += Scurve_centers
                # one side away
                pre_center = centers_list[-1]
                next_CircleCenters = tools._gen_straight_nCircleCenters(pre_center=pre_center, n=self.away, R=self.R, upward=False)
                centers_list += next_CircleCenters[1:]
                self.centers = centers_list
                return self
                    
            else:
                # short Nterm
                # downwarding straight
                pre_center = centers_list[-1]
                next_CircleCenters = tools._gen_straight_nCircleCenters(pre_center=pre_center, n=restCircles_comb, R=self.R, upward=False)
                centers_list += next_CircleCenters[1:]
                self.centers = centers_list
                return self
    
    def genTMCenters(self, idx):
        TMCircleCenters_DICT = self.TMCircleCenters_DICT
        centers_list = self.centers
        pre_center = centers_list[-1]
        R = self.R
        upward = True 

        if idx%4 == 1:
            # upwarding
            upTMCenters_relative = TMCircleCenters_DICT[idx]
            upTM_start_center = tools._next_circle_center(pre_center, R, upward, 90)
            upTMCenters = tools.MoveCoords(coords=upTMCenters_relative, new_start=upTM_start_center, first=True)
            centers_list += upTMCenters[1:]
            self.centers = centers_list
            return self
        else:
            # idx%4 == 3
            downTMCenters_relative = TMCircleCenters_DICT[idx]
            downTM_start_center = tools._next_circle_center(pre_center, R, operator.not_(upward), 90)
            centers_list += [downTM_start_center]
            downTMCenters = tools.MoveCoords(coords=downTMCenters_relative, new_start=downTM_start_center, first=False)
            centers_list += downTMCenters[1:]
            self.centers = centers_list
            return self
    
    def genExtracellularCenters(self, length):
        centers_list = self.centers
        pre_center = centers_list[-1]
        extracellular_NotTMCenters = tools._genNotTMCenters(pre_center, length, away=self.away, R=self.R, extracellular=True, Nterm=False, Cterm=False)
        centers_list += extracellular_NotTMCenters
        self.centers = centers_list
        return self
 

