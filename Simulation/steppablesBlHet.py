from PySteppables import *
import CompuCell
import CompuCellSetup
import sys
from PySteppablesExamples import MitosisSteppableBase
from random import *
import math
import time
# from math import *
import numpy
from datetime import datetime
import cProfile
# from scipy import integrate




#STRINGS FOR OUTPUT FILES
#drug='cisplatin'
#size='lgSphere'
#vess='451v_0.99pctA'


#INITIALIZE COUNTING VARIABLES
#numPCancer=0
#numIC50Cancer=0



## TIME FRAMES

# CONSTANTS RELATED TO DIFFUSION COEFFICIENTS
CisGem1Min = 65.678 # 65.678 mcs = 1 min of diffusion time for cells of diameter T24 bladder cancer cell line, for drugs with diffusion coeff. of sodium fluorescein
MCSFractionOfHour = 0.0002537615293 # hours per MCS, based on diffusion time for one T24 cell diameter of sodium fluorescein, proxy for cisplatin and gemcitabine

# CELL TIMINGS
divisionCycleTimeHrs = 30 # average time to division / replication from several cancer cell lines in vitro
phagocytosisEndTime = 24 # dead cells removed at 24 hours
# divisionCycleTimeHrs = 0.001 # TEST average time to division / replication from several cancer cell lines in
# phagocytosisEndTime = 0.1 # TEST dead cells removed at x hours



#IV CISPLATIN
# 70 mg/m^2
cis70EndInfus=180*CisGem1Min # end of 3-h infusion (de Jongh, 2001)

# 60 mg/m^2
drug15Mins=15*CisGem1Min # 18107.745 MCS; subtract from runtime once begin using concentration time course
cIVFirstPoint=(5.742+15)*CisGem1Min # 6931.79 MCS; five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
cEndDataSet=(170.862*CisGem1Min)+cIVFirstPoint
cFirst5Mins=5*CisGem1Min # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
cFirst20Mins=20*CisGem1Min # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
# cFirst5Mins=5*465.189 # five minutes worth of MCSes in normal tissue (Swabb 1974?)
cisplatinIC50=52.71 # FINAL USED FOR OVARIAN MODEL: muM (equiv to (equitoxic) 2h-IC-50 Mistry, 1992
# cisplatinIC50=38.3 # (MATCHES 2H DRUG TIME COURSE): muM (SD=12.6, in SKOV-3, 2h exposure; Table 1, Mistry, 1992)
# cisplatinIC50=126.63 # muM (+/- 12.06 micromols/L, in SKOV3ip1, 48h exposure; Fig. 3, Xu, 2008)
# cisplatinIC50  =3.33 # muM (Nakaro 1997, minimum effective concentration for DNA damage)

# IV GEMCITABINE
drug30Mins=30*CisGem1Min # = 1970.353826 mcs
gemZeroConcTime=240*CisGem1Min # time at final data point

# REGIMEN TIMINGS
# aggressive bladder cancer regimen time frame, NCCN regimen for metastatic bl.canc. (list)
#MCS for each step of Gem-Cis aggressive regimen; MCS at the end of each stage of regimen cycle:day:eventsTotal; e.g. gem infusion cycle 1: day 18: 3 infusions total to this point
#	[MCS start of gem infusion 1:1:1, end of gem infusion 1:1:1,
#        MCS end of week 1,	MCS end of gem infusion 1:8:2,
#        MCS end of week 2,	MCS end of gem infusion 1:15:3
#        MCS end of cycle 1, 21d
# aggressInfusTimesGem = numpy.array([0, 15762.8306,
#                         662038.8854, 677801.716,
#                         1324077.771, 1339840.601,
#                         1986116.656])
aggressInfusTimesGem = numpy.array([0, 15762.8306,
                        94576.98363 ,110339.8142,
                        189153.9673 ,204916.7979,
                        1986116.656])
cycleTime = 21.0*94576.98363 # 21 days * mcs/h
# aggressInfusTimesGem = numpy.array([0,1,2,3,4,5,10.0]) # test
# cycleTime = 10.0 # test
# print 'cycleTime',cycleTime
# print 'aggressInfusTimesGem',aggressInfusTimesGem
# aggressInfusTimesGem = aggressInfusTimesGem + cycleTime
# print 'aggressInfusTimesGem plus 21',aggressInfusTimesGem


# INFUSION SIMULATANEOUSLY WITH GEM
aggressInfusTimeDay1Cis = numpy.array([0.0, 27584.95356])	#MCS, 7h; end of gem infusion 1:1:1 to approximate cis = 0(last data point at 6h, exponential fit should approach 0 between 6 and 7 h) for de Jongh 2001, 70 mg/m^2
# INFUSION DIRECTLY AFTER GEM
# aggressInfusTimeDay1Cis = numpy.array([15762.8306, 43347.78416])	#MCS, 7h; end of gem infusion 1:1:1 to approximate cis = 0(last data point at 6h, exponential fit should approach 0 between 6 and 7 h) for de Jongh 2001, 70 mg/m^2
#aggressInfusTimeDay1Cis = [15762.8306, 28531.83993]	#MCS, end of gem infusion 1:1:1 to approximate cis = 0  for Casper 1984, 60 mg/m^2



# CELLULAR PARAMETERS: IC50, ACCUMULATION OF DRUG
## CELL IDs
# TypeId="4" TypeName="SCSG_BFTC_905"
# TypeId="5" TypeName="SCSG_J82"
# TypeId="6" TypeName="RCRG_RT4"
# TypeId="7" TypeName="RCRG_HT_1197"
# TypeId="8" TypeName="SCRG_SW780"
# TypeId="9" TypeName="SCRG_KU_19_19"
# TypeId="10" TypeName="RCSG_LB831_BLC"
# TypeId="11" TypeName="RCSG_DSH1"
# TypeId="12" TypeName="IC50Cis"
# TypeId="13" TypeName="IC50Gem"
# TypeId="14" TypeName="cisResistant"
# TypeId="15" TypeName="gemResistant"
# TypeId="16" TypeName="dualResist"


## ACCUMULATION RATES
## IC50 cells have same accumulation rate as cell line they came from (rate is carried as part of dictionary)
## Dead cells and LungNormal are set in the SecretionSteppables to be same as middle-of-the-road-least-sensitive line SCRG_SW780
## Dead cells presumedly have active macrophages in their space collecting drug for 24 hours after death -- no fit to a long-tailed distribution, since proximity of macrophages to multiple cells is likely correlated (they probably all get eatn at close to the same time)
## cisplatin, platinum accumulation per cell per time step, based on IC50 of bladder cancer cell line;
##            see paper Table for fits of IC50s to accumulations
## drugAccumFrac_x = concentration removed from voxel: microM/MCS, * siteConcCis(microM) in SecretionSteppable
cispAccumFrac_SCSG_BFTC_905 = 7.98701E-05       # (sens cis and gem)	2.575477619	IC50 microM	cisplatin
cispAccumFrac_SCSG_J82 = 7.69840E-05            # (sens cis and gem)	5.42972235	IC50 microM	cisplatin				
cispAccumFrac_RCRG_RT4 = 5.46716E-05            # (resist cis and gem)	27.49620513	IC50 microM	cisplatin				
cispAccumFrac_RCRG_HT_1197 = 3.80138E-05        # (resist cis and gem)	43.97041406	IC50 microM	cisplatin				
cispAccumFrac_SCRG_SW780 = 6.82909E-05          # (sens cis resist gem)	14.02708355	IC50 microM	cisplatin					
cispAccumFrac_SCRG_KU_19_19 = 7.22571E-05       # (sens cis resist gem)	10.10456195	IC50 microM	cisplatin
cispAccumFrac_RCSG_LB831_BLC = 7.42347E-06      # (resist cis sens gem)	225.1619678	IC50 microM	cisplatin
cispAccumFrac_RCSG_DSH1 = 1.15622E-05           # (resist cis sens gem)	144.5646771	IC50 microM	cisplatin				

## gemcitabine (possibly dFdCtP) accumulation per cell per time step, based on IC50 of bladder cancer cell line;
#  siteConcGem * microM/MCS = (-0.8242 * IC50 + 67.2261) * siteConcGem/50 * 1/1.5E6 * 1/10^9 * 1/$B$6 * $B$9 * 10^6    =	microM gem accumulation / MCS * frac50uMGem
gemAccumFrac_SCSG_BFTC_905 = 4.41575E-04     # (sensitive cis and gem)	5.15E-06	IC50 microM	gemcitabine			
gemAccumFrac_SCSG_J82 = 4.41565E-04          # (sensitive cis and gem)
gemAccumFrac_RCRG_RT4 = 4.22858E-04          # (resistant cis and gem)	13.84278281	IC50 microM	gemcitabine				
gemAccumFrac_RCRG_HT_1197 = 4.39190E-04      # (resistant cis and gem)	1.764248816	IC50 microM	gemcitabine				
gemAccumFrac_SCRG_SW780 = 2.68443E-04        # (sens cis resist gem)	128.0487435	IC50 microM	gemcitabine				
gemAccumFrac_SCRG_KU_19_19 = 4.18843E-04     # (sens cis resist gem)	16.81235445	IC50 microM	gemcitabine				
gemAccumFrac_RCSG_LB831_BLC = 4.41518E-04    # (resist cis sens gem)	0.041854289	IC50 microM	gemcitabine				
gemAccumFrac_RCSG_DSH1 = 4.41445E-04         # (resist cis sens gem)	0.096498675	IC50 microM	gemcitabine

## IC50s
# from the GDSC database, June 2016
cisIC50_SCSG_BFTC_905 = 0.8106177157         # (sens cis and gem)	2.575477619	IC50 microM	cisplatin
cisIC50_SCSG_J82 = 1.647223153               # (sens cis and gem)	5.42972235	IC50 microM	cisplatin				
cisIC50_RCRG_RT4 = 5.923917064               # (resist cis and gem)	27.49620513	IC50 microM	cisplatin				
cisIC50_RCRG_HT_1197 = 6.586828431           # (resist cis and gem)	43.97041406	IC50 microM	cisplatin				
cisIC50_SCRG_SW780 = 3.774888444             # (sens cis resist gem)	14.02708355	IC50 microM	cisplatin					
cisIC50_SCRG_KU_19_19 = 2.877213884          # (sens cis resist gem)	10.10456195	IC50 microM	cisplatin
cisIC50_RCSG_LB831_BLC = 6.586828431         # (resist cis sens gem)	225.1619678	IC50 microM	cisplatin
cisIC50_RCSG_DSH1 = 6.586828431              # (resist cis sens gem)	144.5646771	IC50 microM	cisplatin				

gemIC50_SCSG_BFTC_905 = 0.000017923         # (sensitive cis and gem)	5.15E-06	IC50 microM	gemcitabine			
gemIC50_SCSG_J82 =  0.025787266             # (sensitive cis and gem)
gemIC50_RCRG_RT4 =  46.134163935           # (resistant cis and gem)	13.84278281	IC50 microM	gemcitabine				
gemIC50_RCRG_HT_1197 = 6.106834031         # (resistant cis and gem)	1.764248816	IC50 microM	gemcitabine				
gemIC50_SCRG_SW780 = 270.913928515        # (sens cis resist gem)	128.0487435	IC50 microM	gemcitabine				
gemIC50_SCRG_KU_19_19 = 55.498902953       # (sens cis resist gem)	16.81235445	IC50 microM	gemcitabine				
gemIC50_RCSG_LB831_BLC = 0.145644144        # (resist cis sens gem)	0.041854289	IC50 microM	gemcitabine				
gemIC50_RCSG_DSH1 = 0.335738949             # (resist cis sens gem)	0.096498675	IC50 microM	gemcitabine

## CELL VOLUME PARAMETERS
T24BCCellVol = 1 # bladder cancer cell volume (units = voxels)
normalLambdaVolume = 100.0
cellGrowthLambdaVolume = 90.0 # =90.0, others higher (100.0) to keep dividing cells from replacing pre-existing cells
phagocytosisLambdaVolume = 1000.0
deadLambdaVolume = 100.0

## VASCULARITY
vesselPercentInMetastasis = 0.146 # 0.1460592054 = fraction of vessels per area in bladder cancer metastases, estimated from CLCC ratio of metastatic MVD/primary MVD and bladder cancer primary MVD (microvessel density = MVD)

## use if using percent vessel in total sim area to limit vessel growth; not in use as of 6-30-2016
# totalSimCellsPossible = 20*20*1 #CHANGE WITH SIM DIMENSIONS!
# print "total cells in sim =",totalSimCellsPossible
# maxVesselCellCount = round(vesselPercentInMetastasis*totalSimCellsPossible) # used in mitosis to limit global vessel nums






# PRINT SIMULATION START TIME
now = datetime.now()
print "SIMULATION START TIME =",now




# SETCELLCLOCKS
# SET FOR INITIAL CELL TYPES; IC50 AND RESISTANT CELL TYPES WILL INHERIT BASELINE CHARACTERISTICS (CYCLE TIME, IC50, ACCUMULATION RATE, ETC.)
# SETS UNIQUE INITIAL AGE FROM UNIFORM DISTRIBUTION
# STARTS INTERNAL CLOCK IN CELL
# SETS DIVISION RATE FROM GAUSSIAN DISTRIBUTION
class SetCellDictionaries(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (SetCellDictionaries) is called once before simulation"
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:
            x = gauss(divisionCycleTimeHrs,1/30*divisionCycleTimeHrs) # assumes that div time in final sims is 30h
            y = uniform(0,divisionCycleTimeHrs) # age of cells initialized into simulation
            cell.dict["AgeHrs"]=y
            # cell.dict["AgeHrs"]=0
            cell.dict["cycleHrs"]=x
            cell.dict["generation"]=0
            cell.dict["numDivisions"] = 0
            cell.dict["HrsSinceDeath"]=0
            cell.dict["cisAccum"]=0
            cell.dict["gemAccum"]=0
            cell.dict["cisResistance"]=1
            cell.dict["gemResistance"]=1
            cell.dict["IC50CisOrig"]=0
            cell.dict["IC50GemOrig"]=0
            cell.dict["IC50Cis"]=0
            cell.dict["IC50Gem"]=0
            cell.dict["accumRtCis"]=0
            cell.dict["accumRtGem"]=0
            # ## NO DRUG SYNERGY
            # if cell.type==4:
            #     cell.dict["IC50CisOrig"]=cisIC50_SCSG_BFTC_905
            #     cell.dict["IC50GemOrig"]=gemIC50_SCSG_BFTC_905
            #     cell.dict["IC50Cis"]=cisIC50_SCSG_BFTC_905
            #     cell.dict["IC50Gem"]=gemIC50_SCSG_BFTC_905
            #     cell.dict["accumRtCis"]=cispAccumFrac_SCSG_BFTC_905
            #     cell.dict["accumRtGem"]=gemAccumFrac_SCSG_BFTC_905
            # if cell.type==5:
            #     cell.dict["IC50CisOrig"]=cisIC50_SCSG_J82
            #     cell.dict["IC50GemOrig"]=gemIC50_SCSG_J82
            #     cell.dict["IC50Cis"]=cisIC50_SCSG_J82
            #     cell.dict["IC50Gem"]=gemIC50_SCSG_J82
            #     cell.dict["accumRtCis"]=cispAccumFrac_SCSG_J82
            #     cell.dict["accumRtGem"]=gemAccumFrac_SCSG_J82
            # if cell.type==6:
            #     cell.dict["IC50CisOrig"]=cisIC50_RCRG_RT4
            #     cell.dict["IC50GemOrig"]=gemIC50_RCRG_RT4
            #     cell.dict["IC50Cis"]=cisIC50_RCRG_RT4
            #     cell.dict["IC50Gem"]=gemIC50_RCRG_RT4
            #     cell.dict["accumRtCis"]=cispAccumFrac_RCRG_RT4
            #     cell.dict["accumRtGem"]=gemAccumFrac_RCRG_RT4
            # if cell.type==7:
            #     cell.dict["IC50CisOrig"]=cisIC50_RCRG_HT_1197
            #     cell.dict["IC50GemOrig"]=gemIC50_RCRG_HT_1197
            #     cell.dict["IC50Cis"]=cisIC50_RCRG_HT_1197
            #     cell.dict["IC50Gem"]=gemIC50_RCRG_HT_1197
            #     cell.dict["accumRtCis"]=cispAccumFrac_RCRG_HT_1197
            #     cell.dict["accumRtGem"]=gemAccumFrac_RCRG_HT_1197
            # if cell.type==8:
            #     cell.dict["IC50CisOrig"]=cisIC50_SCRG_SW780
            #     cell.dict["IC50GemOrig"]=gemIC50_SCRG_SW780
            #     cell.dict["IC50Cis"]=cisIC50_SCRG_SW780
            #     cell.dict["IC50Gem"]=gemIC50_SCRG_SW780
            #     cell.dict["accumRtCis"]=cispAccumFrac_SCRG_SW780
            #     cell.dict["accumRtGem"]=gemAccumFrac_SCRG_SW780
            # if cell.type==9:
            #     cell.dict["IC50CisOrig"]=cisIC50_SCRG_KU_19_19
            #     cell.dict["IC50GemOrig"]=gemIC50_SCRG_KU_19_19
            #     cell.dict["IC50Cis"]=cisIC50_SCRG_KU_19_19
            #     cell.dict["IC50Gem"]=gemIC50_SCRG_KU_19_19
            #     cell.dict["accumRtCis"]=cispAccumFrac_SCRG_KU_19_19
            #     cell.dict["accumRtGem"]=gemAccumFrac_SCRG_KU_19_19
            # if cell.type==10:
            #     cell.dict["IC50CisOrig"]=cisIC50_RCSG_LB831_BLC
            #     cell.dict["IC50GemOrig"]=gemIC50_RCSG_LB831_BLC
            #     cell.dict["IC50Cis"]=cisIC50_RCSG_LB831_BLC
            #     cell.dict["IC50Gem"]=gemIC50_RCSG_LB831_BLC
            #     cell.dict["accumRtCis"]=cispAccumFrac_RCSG_LB831_BLC
            #     cell.dict["accumRtGem"]=gemAccumFrac_RCSG_LB831_BLC
            # if cell.type==11:
            #     cell.dict["IC50CisOrig"]=cisIC50_RCSG_DSH1
            #     cell.dict["IC50GemOrig"]=gemIC50_RCSG_DSH1
            #     cell.dict["IC50Cis"]=cisIC50_RCSG_DSH1
            #     cell.dict["IC50Gem"]=gemIC50_RCSG_DSH1
            #     cell.dict["accumRtCis"]=cispAccumFrac_RCSG_DSH1
            #     cell.dict["accumRtGem"]=gemAccumFrac_RCSG_DSH1

            ##DRUG SYNERGY: PRE-TREATMENT AND CO-TREATMENT WITH GEMCITABINE IMPROVES CISPLATIN EFFICACY, ~2.5X (Moufarij, 2003)
            ##remove "*2.5" from any cell line's cisplatin accumulation to make it non-synergistic
            if cell.type==4:
                cell.dict["IC50CisOrig"]=cisIC50_SCSG_BFTC_905
                cell.dict["IC50GemOrig"]=gemIC50_SCSG_BFTC_905
                cell.dict["IC50Cis"]=cisIC50_SCSG_BFTC_905
                cell.dict["IC50Gem"]=gemIC50_SCSG_BFTC_905
                cell.dict["accumRtCis"]=cispAccumFrac_SCSG_BFTC_905*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_SCSG_BFTC_905
            if cell.type==5:
                cell.dict["IC50CisOrig"]=cisIC50_SCSG_J82
                cell.dict["IC50GemOrig"]=gemIC50_SCSG_J82
                cell.dict["IC50Cis"]=cisIC50_SCSG_J82
                cell.dict["IC50Gem"]=gemIC50_SCSG_J82
                cell.dict["accumRtCis"]=cispAccumFrac_SCSG_J82*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_SCSG_J82
            if cell.type==6:
                cell.dict["IC50CisOrig"]=cisIC50_RCRG_RT4
                cell.dict["IC50GemOrig"]=gemIC50_RCRG_RT4
                cell.dict["IC50Cis"]=cisIC50_RCRG_RT4
                cell.dict["IC50Gem"]=gemIC50_RCRG_RT4
                cell.dict["accumRtCis"]=cispAccumFrac_RCRG_RT4*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_RCRG_RT4
            if cell.type==7:
                cell.dict["IC50CisOrig"]=cisIC50_RCRG_HT_1197
                cell.dict["IC50GemOrig"]=gemIC50_RCRG_HT_1197
                cell.dict["IC50Cis"]=cisIC50_RCRG_HT_1197
                cell.dict["IC50Gem"]=gemIC50_RCRG_HT_1197
                cell.dict["accumRtCis"]=cispAccumFrac_RCRG_HT_1197*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_RCRG_HT_1197
            if cell.type==8:
                cell.dict["IC50CisOrig"]=cisIC50_SCRG_SW780
                cell.dict["IC50GemOrig"]=gemIC50_SCRG_SW780
                cell.dict["IC50Cis"]=cisIC50_SCRG_SW780
                cell.dict["IC50Gem"]=gemIC50_SCRG_SW780
                cell.dict["accumRtCis"]=cispAccumFrac_SCRG_SW780*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_SCRG_SW780
            if cell.type==9:
                cell.dict["IC50CisOrig"]=cisIC50_SCRG_KU_19_19
                cell.dict["IC50GemOrig"]=gemIC50_SCRG_KU_19_19
                cell.dict["IC50Cis"]=cisIC50_SCRG_KU_19_19
                cell.dict["IC50Gem"]=gemIC50_SCRG_KU_19_19
                cell.dict["accumRtCis"]=cispAccumFrac_SCRG_KU_19_19*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_SCRG_KU_19_19
            # **** non-synergistic 7-18-2016 (model has since changed, needs to be re-run)
            if cell.type==10:
                cell.dict["IC50CisOrig"]=cisIC50_RCSG_LB831_BLC
                cell.dict["IC50GemOrig"]=gemIC50_RCSG_LB831_BLC
                cell.dict["IC50Cis"]=cisIC50_RCSG_LB831_BLC
                cell.dict["IC50Gem"]=gemIC50_RCSG_LB831_BLC
                cell.dict["accumRtCis"]=cispAccumFrac_RCSG_LB831_BLC*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_RCSG_LB831_BLC
            if cell.type==11:
                cell.dict["IC50CisOrig"]=cisIC50_RCSG_DSH1
                cell.dict["IC50GemOrig"]=gemIC50_RCSG_DSH1
                cell.dict["IC50Cis"]=cisIC50_RCSG_DSH1
                cell.dict["IC50Gem"]=gemIC50_RCSG_DSH1
                cell.dict["accumRtCis"]=cispAccumFrac_RCSG_DSH1*2.5
                cell.dict["accumRtGem"]=gemAccumFrac_RCSG_DSH1

            # print initial dictionary vals for each cell
            # print 'cell.type=',cell.type,'cell.id=',cell.id,'dict=',cell.dict
        # "TO GET ALL CELL ATTRIBUTES" (to see what cell attributes can be accessed/changed in Python):
        print 'Members of SteppableBasePy class'
        print dir(cell)

            # break




# *****************************
# INCREMENT AGE
# INCREMENT TIME SINCE DEATH
# INCREMENT DRUG DELIVERY TIMES IF IN A NEW DRUG CYCLE
class IncrementClocks(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
        self.simulator=_simulator

    def start(self):
        print "This function (IncrementClocks) is called at every MCS"

    def step(self,mcs):
        # # increment cycle time if next cycle has been entered
        if mcs >= aggressInfusTimesGem[6]: # indexing starts at "0"
            print 'last time point gem dosages array',aggressInfusTimesGem[6]
            global aggressInfusTimesGem
            print 'aggressInfusTimesGem',aggressInfusTimesGem
            aggressInfusTimesGem = aggressInfusTimesGem + cycleTime
            print 'aggressInfusTimesGem plus 21',aggressInfusTimesGem
            print 'cycle time mcs added = ', cycleTime
            print 'Gemcitabine 30 min infusion end mcs = ',drug30Mins
        if mcs >= aggressInfusTimeDay1Cis[1]:
            global aggressInfusTimeDay1Cis
            print 'aggressInfusTimesDay1Cis',aggressInfusTimeDay1Cis
            aggressInfusTimeDay1Cis = aggressInfusTimeDay1Cis + cycleTime
            print 'aggressInfusTimesDay1Cis plus 21',aggressInfusTimeDay1Cis
            print 'cycle time mcs added = ', cycleTime
            print 'Cisplatin 3h infusion end mcs = ',cis70EndInfus

        self.cellList=CellList(self.inventory)
        for cell in self.cellList:
            cell.dict["AgeHrs"]+= MCSFractionOfHour
            if cell.type==3:
                cell.dict["HrsSinceDeath"]+= MCSFractionOfHour
            # if cell.type!=1:
            # print 'cell.id=',cell.id,'cell.type=',cell.type,' dict=',cell.dict, 'vol=',cell.targetVolume,'volLambda=',cell.lambdaVolume
            # print cell.dict





# ***************************** ALSO UNCOMMENT CELL DICTIONARY KEYS 1 AND 2 FOR CELL AGE ABOVE
# VOLUMEPARAMSTEPPABLE
# SET CELL TARGET VOLUMES OR COMPRESSIBILITIES
# SET VOLUME CONSTRAINTS
class VolumeParamSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (VolumeParamSteppable) is called once before simulation"
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:
                cell.targetVolume=T24BCCellVol
                cell.lambdaVolume=normalLambdaVolume
                # print 'cell.type=',cell.type,'cell.id=',cell.id,'cell.volume=',cell.targetVolume,'cell.lambdaVolume=',cell.lambdaVolume

    # def step(self,mcs):
    #     for cell in self.cellList:
    #         print "MCS",mcs,'cell.type=',cell.type,'cell.id=',cell.id,'cell.volume=',cell.targetVolume,'cell.lambdaVolume=',cell.lambdaVolume




# *****************************
# GROW CELLS WHEN READY TO DIVIDE
# cells grow any time there is an adjacent space if they have passed a division time;
# this means any spaces that open up in the inner tumor volume will be quickly filled.
# cells that are at IC50 have a chance to die only once each time they reach division cycle time
class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            if cell.type!=1 and cell.type!=2 and cell.type!=3: # all cell types grow and then divide except for Vessel, LungNormal, and Dead, respectively (IC50Cis, and IC50Gem divide)
                # if cell.dict["AgeHrs"]>=divisionCycleTimeHrs:
                # print 'inside growthSteppable for tumor cells, type is',cell.type,'age is', cell.dict["AgeHrs"]
                if cell.dict["AgeHrs"]>=cell.dict["cycleHrs"]:
                    cell.targetVolume=2*T24BCCellVol
                    cell.lambdaVolume=cellGrowthLambdaVolume
                    # print 'I am cell.type',cell.type,'cell.id',cell.id,'targetVolume',cell.targetVolume,'targetLambda',cell.lambdaVolume,'and I want to grow so I can divide.'
        # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        



# *****************************
# MITOSIS; IC50 CELLS DIE WITH 50% CHANCE WHEN MITOSIS IS ATTEMPTED; DEAD, IC50, VESSEL, AND LUNG CELLS DO NOT REPLICATE;
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)

    # def start(self):
    #     print "This function (MitosisSteppable) is called at every MCS"

    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            # cell death from drugs
            # if cells are IC50Cis or IC50Gem
            if cell.type==12 or cell.type==13:
                if cell.dict["cycleHrs"]<cell.dict["AgeHrs"]:
                    deathChance = uniform(0,1)
                    print 'deathChance=',deathChance
                    if deathChance<=0.5:
                        cell.type=3 # cell dies with 50% chance at cell division attempt
                        cell.targetVolume = cell.volume # keep cell same size as when died so it won't expand over other cells
                        cell.lambdaVolume=deadLambdaVolume
                        print 'cell.type', cell.type,'cell.id', cell.id, 'died'
                    else: # become dually resistant if you a. have gained resistance to both drugs or b. you survive being IC50 to both drugs at once (2nd chance at death)
                        print 'cell.type', cell.type,'cell.id', cell.id, 'will increase its resistance'
                        #cisIC50
                        if cell.type==12:
                            if cell.dict["gemResistance"] > 1 or (cell.dict["cisAccum"] > cell.dict["IC50Cis"] and cell.dict["gemAccum"] > cell.dict["IC50Gem"]):
                                if (cell.dict["cisAccum"] > cell.dict["IC50Cis"] and cell.dict["gemAccum"] > cell.dict["IC50Gem"]):
                                    deathChance = uniform(0,1)
                                    print 'deathChance for dual IC50=',deathChance
                                    if deathChance<=0.5:
                                        cell.type=3 # cell dies with 50% chance at cell division attempt
                                        cell.targetVolume = cell.volume # keep cell same size as when died so it won't expand over other cells
                                        cell.lambdaVolume=deadLambdaVolume
                                        print 'cell.type', cell.type,'cell.id', cell.id, 'died'
                                else:
                                    cell.type = 16 # dualResist
                            else:
                                cell.type = 14 # CisResist
                            if cell.dict["cisResistance"] < 29: # Max multiple of IC50 in cell lines gaining resistance to cisplatin within 1-2yrs culturing; Vallo et al., 2015
                                # cell.dict["cisResistance"] += 0 # zero gain
                                # cell.dict["cisResistance"] += 1 # simple gain
                                # cell.dict["cisResistance"] += 0.05 # slow gain (2 y of 30h cell cycles to max)
                                # cell.dict["cisResistance"] += 0.1 # fast gain (1 y of 30h cell cycles to max)
                                gain = uniform(0.05,0.1)
                                print 'cisGain (range 0.05 - 0.1) =',gain
                                cell.dict["cisResistance"] += gain # uniform betw slow and fast gain (1 to 2 y of 30h cell cycles to max)
                                cell.dict["IC50Cis"] = cell.dict["IC50CisOrig"] * cell.dict["cisResistance"]
                            print 'cell.type', cell.type,'cell.id', cell.id, 'increased its cis resistance'
                        # gemIC50
                        if cell.type==13:
                            if cell.dict["cisResistance"] > 1 or (cell.dict["cisAccum"] > cell.dict["IC50Cis"] and cell.dict["gemAccum"] > cell.dict["IC50Gem"]):
                                if (cell.dict["cisAccum"] > cell.dict["IC50Cis"] and cell.dict["gemAccum"] > cell.dict["IC50Gem"]):
                                    deathChance = uniform(0,1)
                                    print 'deathChance for dual IC50=',deathChance
                                    if deathChance<=0.5:
                                        cell.type=3 # cell dies with 50% chance at cell division attempt
                                        cell.targetVolume = cell.volume # keep cell same size as when died so it won't expand over other cells
                                        cell.lambdaVolume=deadLambdaVolume
                                        print 'cell.type', cell.type,'cell.id', cell.id, 'died'
                                else:
                                    cell.type = 16 # dualResist
                            else:
                                cell.type = 15 # GemResist
                            if cell.dict["gemResistance"] < 73: # Max multiple of IC50 in cell lines gaining resistance to gemcitabine within 1-2yrs culturing; Vallo et al., 2015
                                # cell.dict["gemResistance"] += 0 # zero gain
                                # cell.dict["gemResistance"] += 1 # simple gain
                                # cell.dict["gemResistance"] += 0.125 # slow gain (2 y of 30h cell cycles to max)
                                # cell.dict["gemResistance"] += 0.25 # fast gain (1 y of 30h cell cycles to max)
                                gain = uniform(0.125,0.25)
                                print 'gemGain (range 0.125 - 0.25) =',gain
                                cell.dict["gemResistance"] += gain # uniform betw slow and fast gain (1 to 2 y of 30h cell cycles to max)
                                cell.dict["IC50Gem"] = cell.dict["IC50GemOrig"] * cell.dict["gemResistance"]
                                print 'cell.type', cell.type,'cell.id', cell.id, 'increased its gem resistance'
                                
                        cell.dict["AgeHrs"] = 0 # reset cell cycle; cells that haven't grown don't have a chance to try to divide again; would have had to have doubled size as IC50 type, before becoming current resistant type

            if cell.volume==2*T24BCCellVol: # cells only double in size if they have reached their division time, and only divide if they have doubled in size
                cells_to_divide.append(cell) # if cell is already dead but doubled size, it won't divide below
            # print 'celltype',cell.type,'cellid',cell.id,'is dividing at AgeHrs',cell.dict["AgeHrs"]
        for cell in cells_to_divide:
            if cell.type!=1 and cell.type!=2 and cell.type!=3: # all cell types divide except for Vessel, LungNormal, Dead, respectively (IC50Cis, and IC50Gem divide if they are 2x volume and have lived to become resistant type)
                # to change mitosis mode leave one of the below lines uncommented
                # print 'cells to divide increment'
                self.divideCellRandomOrientation(cell)
                # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
                # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
                # self.divideCellAlongMinorAxis(cell)                               # this is a valid option
        # print 'cells to divide == '
        # print cells_to_divide

    def updateAttributes(self):
        #        if self.parentCell.dict["generation"]<=4:  # if cell has already divided once, do not become vessel
        # print "number of cells that are type vessel:", len(self.cellListByType(self.VESSEL))
        # print "total cells:", len(self.cellList)
        vesselCells = float(len(self.cellListByType(self.VESSEL)))
        totalCells = float(len(self.cellList))
        percentVesselCellsInAllCells = (vesselCells / float(totalCells))
        # print 'vessel count = ',vesselCells,'total cell count = ', totalCells, 'percent vessel = ',percentVesselCellsInAllCells,'compared to', vesselPercentInMetastasis
        # cell divides if volume has doubled (condition for GrowthSteppable)
        if percentVesselCellsInAllCells < vesselPercentInMetastasis: # prev used maxVesselCellCount based on sim dimensions
            # if self.parentCell.dict["generation"] > 5:
            chanceToBeVessel = uniform(0,1)
            if chanceToBeVessel <= vesselPercentInMetastasis:
                print 'mitosis, vessel = child'
                self.parentCell.targetVolume /= 2.0 # reduce parent target volume by increasing; = ratio to parent vol
                self.parentCell.lambdaVolume = normalLambdaVolume
                self.parentCell.dict["AgeHrs"] = 0 # re-set cell to keep distribution of vessel more even in sim field -- no more vessel in this region for the time of a cell cycle -- when using % vessel in overall space as control for vascular density
                self.parentCell.dict["numDivisions"] += 1
                self.cloneParent2Child() # copy all parent parameters, then over-write
                self.childCell.type=1 # CHILD IS VESSEL
                self.childCell.targetVolume = 1
                self.childCell.dict["generation"]+=1
                self.childCell.dict["numDivisions"] = 0
                self.childCell.dict["cisAccum"] = 0
                self.childCell.dict["gemAccum"] = 0
                print 'childCell.type=',self.childCell.type, 'childCell.id=',self.childCell.id,' dict=',self.childCell.dict,'childCell.targetVolume=', self.childCell.targetVolume,'childCell.lambdaVolume=', self.childCell.lambdaVolume
                print   'parentCell.type=',self.parentCell.type, 'parentCell.id=',self.parentCell.id,' dict=',self.parentCell.dict,'parentCell.targetVolume=', self.parentCell.targetVolume,'parentCell.lambdaVolume=', self.parentCell.lambdaVolume
            else:
                print 'mitosis, chance at vessel failed'
                self.parentCell.targetVolume /= 2.0 # reduce parent target volume by increasing; = ratio to parent vol
                self.parentCell.lambdaVolume = normalLambdaVolume # make sure parent stays in place
                self.parentCell.dict["AgeHrs"] = 0
                self.parentCell.dict["numDivisions"] += 1
                self.cloneParent2Child() # copy all parent parameters, then over-write
                self.childCell.dict["generation"]+=1
                self.childCell.dict["numDivisions"] = 0
                self.childCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
                self.childCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
                self.parentCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
                self.parentCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
                ## for cell in self.cellList:
                print  'childCell.type=',self.childCell.type,'childCell.id=',self.childCell.id,' dict=',self.childCell.dict,'childCell.targetVolume=', self.childCell.targetVolume,'childCell.lambdaVolume=', self.childCell.lambdaVolume
                print  'parentCell.type=',self.parentCell.type,  'parentCell.id=',self.parentCell.id,' dict=',self.parentCell.dict,'parentCell.targetVolume=', self.parentCell.targetVolume,'parentCell.lambdaVolume=', self.parentCell.lambdaVolume
        else:
            print 'mitosis, no chance at vessel'
            self.parentCell.targetVolume /= 2.0 # reduce parent target volume by increasing; = ratio to parent vol
            self.parentCell.lambdaVolume = normalLambdaVolume # make sure parent stays in place
            self.parentCell.dict["AgeHrs"] = 0
            self.parentCell.dict["numDivisions"] += 1
            self.cloneParent2Child() # copy all parent parameters, then over-write
            self.childCell.dict["generation"]+=1
            self.childCell.dict["numDivisions"] = 0
            self.childCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
            self.childCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
            self.parentCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
            self.parentCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
            ## for cell in self.cellList:
            print  'childCell.type=',self.childCell.type, 'childCell.id=',self.childCell.id,' dict=',self.childCell.dict,'childCell.targetVolume=', self.childCell.targetVolume,'childCell.lambdaVolume=', self.childCell.lambdaVolume
            print 'parentCell.type=',self.parentCell.type, 'parentCell.id=',self.parentCell.id,' dict=',self.parentCell.dict,'parentCell.targetVolume=', self.parentCell.targetVolume,'parentCell.lambdaVolume=', self.parentCell.lambdaVolume


        # uncomment when adding in the generational limitation on mitosis to vessel
        # else:
        #     print 'mitosis, no chance at vessel'
        #     self.parentCell.targetVolume /= 2.0 # reduce parent target volume by increasing; = ratio to parent vol
        #     self.parentCell.lambdaVolume = normalLambdaVolume # make sure parent stays in place
        #     self.parentCell.dict["AgeHrs"] = 0
        #     self.parentCell.dict["numDivisions"] += 1

        #     self.cloneParent2Child() # copy all parent parameters, then over-write
        #     self.childCell.dict["generation"]+=1
        #     self.childCell.dict["numDivisions"] += 0
        #     self.childCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
        #     self.childCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
        #     self.parentCell.dict["cisAccum"] = 0.5 * self.parentCell.dict["cisAccum"]
        #     self.parentCell.dict["gemAccum"] = 0.5 * self.parentCell.dict["gemAccum"]
        #     ## for cell in self.cellList:
        #     print  'childCell.type=',self.childCell.type, 'childCell.id=',self.childCell.id,' dict=',self.childCell.dict,'childCell.targetVolume=', self.childCell.targetVolume,'childCell.lambdaVolume=', self.childCell.lambdaVolume
        #     print 'parentCell.type=',self.parentCell.type, 'parentCell.id=',self.parentCell.id,' dict=',self.parentCell.dict,'parentCell.targetVolume=', self.parentCell.targetVolume,'parentCell.lambdaVolume=', self.parentCell.lambdaVolume

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )

        # if self.parentCell.type==5:
        #     self.childCell.type=6
        # elif self.parentCell.type==6:
        #     self.childCell.type=7
        #     # else:
        # #     self.childCell.type=




        # *****************************
# CELLS DISAPPEAR AFTER DEATH AND PHAGOCYTOSIS
class RemoveDeadCells(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        print "This function (RemoveDeadCells) is called at every MCS"

    def step(self,mcs):
        for cell in self.cellList:
            if cell.type==3:
                # print 'I am dead cell.id', cell.id, 'and I died',cell.dict["HrsSinceDeath"],'hrs ago.'
                if cell.dict["HrsSinceDeath"]>=phagocytosisEndTime:
                    print 'removing dead cell.id', cell.id,'with cell volume',cell.volume,'cell.targetVolume',cell.targetVolume,'and cell.lambdaVolume',cell.lambdaVolume
                    cell.targetVolume=0
                    cell.lambdaVolume=phagocytosisLambdaVolume




                    # <!-- CISPLATIN -->
  # D(VX2 carcinoma for sodium fluorescein, MW376 (Nugent 1984)
  #      = 1207.18273728686 cell diam^2 / 1 min (= 5.64um^2 (voxel edge) / 1/60hr)
  #      = 1 cell diam^2 / 1/1207.183 min (= 5.64um^2 (voxel edge) / 1/1207.183min)
  #      = 72430.9642372114 MCS/hour = 144861.928474423 2 hours-->
  # <!-- CISPLATIN (PT) ACCUMULATION SET BY FIT TO ACCUMULATION PER IC50 IN BLADDER CANCER CELLS
  #      FOR COMPARISON: = 9.4E-15 $\mu$M/cell/min = 1.57E-16 muM/cell/sec = 9.4E-14 $\mu$M/cell/10 min AT 5 $\MU$M IN SKOV3  (Fig.3, Mistry, 1992)
  #      accumulation per SKOV3 cell, with respect to current concentration (uM)
  #      = (- 1.2e-07*x^{3} + 1.3e-05*x^{2} + 0.00058*x + 0.076)/10/465.189 uM, x = current concentration -->
  #      at 100 micromolar for 50min in SKOV3
  #      = 8 +/- 0.8 nmol/mg protein/50 mins: 8 +/- 0.8 nM * 1e-3uM/nM /mg protein * 0.21mg protein/uL * 1uL/10^9um^3 * 179.1 um^3/ SKOV3 cell /50min = 6.02E-12 uM Pt/SKOV3 cell/min = 1.00E-13 uM Pt/SKOV3 cell/s -->
  # <!-- PEAK AVE INTRAPERITONEAL
  #      (90mg/m^2 IP for 4 hour dwell = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; Howell, 1982):
  #      40.5mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol = 134.977503749375 muM = 8098.65022496251 (concentration * 60 sec/MCS)-->
  # <!-- PEAK VALUES IV
  #      (90mg/m^2 IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; 1.6mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol =  5.3324445925679 uM,
  #      = 319.946675554074 (concentration * 60 sec/MCS) -->
  # <!-- MINIMUM EFFECTIVE CONCENTRATION OF CISPLATIN (BASED ON DNA DAMAGING LEVELS in clear cell ovarian carcinoma, Takatori, 2012)
  #      = 1 mug/mL; 1ug/ml*1000ml/L*1g/10^6ug* 1/300.5g/mol*10^6umol/mol = 3.33 muM-->
  # <!-- MINIMUM EFFECTIVE CONCENTRATION OF CISPLATIN (IN SERUM, clinical trials, lung cancer, Nakano, 1997)
  #      = 0.3micrograms/ml; 0.3ug/ml*1000ml/L*1g/10^6ug* 1/300.05g/mol*10^6umol/mol = 0.999833361106482 muM -->
  # <!-- SECRETION = micromolar = umol/
  #      (90mg/m^2IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL 1.6 1.6mug/mL*1g/10^6mug*1000mL/L/Xg/mol cisplatin*10^6umol/mol 5.3324445925679 uM 319.946675554074 (concentration * 60 sec/MCS) -->
  #      EFFLUX  (Fig.4, Mistry, 1992 ) = 3.58E-04 % accumulated Pt effluxed/ 1/465.189min = 3.58E-06 fraction of accumulated Pt effluxed/ 1/465.189min for D in normal tissue OR 1.38E-06 fraction of accumulated Pt effluxed/ 1/1207.183min for D in tumor tissue
  # <!-- DECAY RATE (not used when empirically determined effective diffusion coefficient is used) = ln(2)/half life = 0.23 / 10 min = 0.000049667766485 /(1/465.189min) (Go, Ata, RXMed:Platinol) -->
# Cisplatin intraperitoneal concentration fitted curve (dosage = 270 mg/m^2 ) 451.51e-0.551*t muM, t=   Howell, 1982 (Fig. 2)
# Cisplatin post-IP plasma concentration fitted curve [IP(t)]*-0.000889*t+0.06288 muM, t= Sugarbaker, 1986, 1990, 1997
# Cisplatin intravenous concentration fitted curve (dosage =100mg/m^2) -3.338ln(t) + 16.094 muM, t=min   Himmelstein, 1981 (Fig. )
# Cisplatin post-IV peritoneal concentration fitted curve [IV(tMins)]*( -2E-07*tMins^3 + 0.0002*tMins^2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver); Linear solution:  (0.5431*t - 0.386) muM, t=h; y = [IV(t)]*0.0091t - 0.386, R^2= 0.9509, t=,min; Sugarbaker, 1996


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class DiffusionSolverFESteeringCisplatinIV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
    def start(self):
        pass
    def step(self,mcs):
        # #TEST
        # print 'mcs=',mcs
        # if -1 < mcs < 20: # FLOATS; USE CONDITIONALS WITHOUT "="
        #     tMins= (mcs + aggressInfusTimeDay1Cis[0]) / CisGem1Min # time since injection
        #     IVtMins = 0.3725*tMins # linear fit for 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
        #     print 'pre-20 IVtMins=',IVtMins
        # elif 20 <= mcs < 40: # plateau for 5.7m
        #     IVtMins = 5.59 # constant for ~6 mins mins; highest and first data point after infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
        #     print 'pre-40 IVtMins=',IVtMins
        # elif 40 <= mcs < 60:        # prior to end of IV data set
        #     tMins=((aggressInfusTimeDay1Cis[1] + mcs)/CisGem1Min) # diffusion time for one cell diameter in tumor tissue; take away added infusion and plateau time so fit is correct; use floats
        #     IVtMins = -1.154e-06*tMins**3 + 0.0005737*tMins**2 - 0.09922*tMins + 5.973 # Casper, 1984
        #     print 'pre-60 IVtMins=',IVtMins

        # print 'conditional inner loop', (aggressInfusTimeDay1Cis[0] + cEndDataSet)
        # print 'conditional outer loop', (aggressInfusTimeDay1Cis[1])
        # print 'test IVtMins=',IVtMins
        # print 'test Python rounding calculations that combine floats and ints:', drug30Mins + aggressInfusTimesGem[2]


        ##### DRUG CONCENTRATIONS AFTER IV DELIVERY:
        if aggressInfusTimeDay1Cis[0] < mcs < aggressInfusTimeDay1Cis[1]: # at correct time in regimen
            if aggressInfusTimeDay1Cis[0] <= mcs < aggressInfusTimeDay1Cis[0] + cis70EndInfus:
                tHrs= (mcs - aggressInfusTimeDay1Cis[0]) / (CisGem1Min * 60) # time since injection
                IVtHrs = 0.11*tHrs**3 - 0.83*tHrs**2 + 2.2*tHrs - 2.6e-16 # cubic fit for 3h infusion (de Jongh, 2001; max = 2.11uM)
                print 'infusion cis IVtHrs=',IVtHrs
            elif aggressInfusTimeDay1Cis[0] + cis70EndInfus <= mcs < aggressInfusTimeDay1Cis[1]:
                tHrs=((mcs - aggressInfusTimeDay1Cis[0]) / (CisGem1Min * 60)) # DO NOT SUBTRACT INFUSION TIME
                IVtHrs = 57.4124 * math.exp(-1.0927 * tHrs) # exponential fit for post-3h infusion (de Jongh, 2001; max = 2.11uM)
                print 'decay cis tHrs=',tHrs
                print 'decay cis IVtHrs= ',IVtHrs

            # #  cisplatin 60 mg/m^2: 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            # if aggressInfusTimeDay1Cis[0] < mcs < aggressInfusTimeDay1Cis[0] + drug15Mins:
            #     tMins= (mcs - aggressInfusTimeDay1Cis[0]) / CisGem1Min # time since injection
            #     IVtMins = 0.3725*tMins # linear fit for 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            #     print 'infusion cis IVtMins=',IVtMins
            # elif aggressInfusTimeDay1Cis[0] + drug15Mins <= mcs < aggressInfusTimeDay1Cis[0] + cIVFirstPoint: # plateau for 5.7m
            #     IVtMins = 5.59 # constant for ~6 mins mins; highest and first data point after infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            #     print 'plateau cis IVtMins=',IVtMins
            # elif aggressInfusTimeDay1Cis[0] + cIVFirstPoint <= mcs < aggressInfusTimeDay1Cis[1]:
            # # aggressInfusTimeDay1Cis[0] + cEndDataSet:        # prior to end of IV data set
            #     tMins=((mcs - aggressInfusTimeDay1Cis[0])/CisGem1Min) - (5.742+15) # take away added infusion and plateau time so fit is correct; use floats
            #     IVtMins = -1.154e-06*tMins**3 + 0.0005737*tMins**2 - 0.09922*tMins + 5.973 # Casper, 1984
            #     if IVtMins < 0:
            #         IVtMins=0 # in case time frame goes past where fit becomes negative
            #     print 'decay cis tMins=',tMins
            #     print 'decay cis IVtMins= ',IVtMins

            # update IV conc
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtHrs             # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

            def finish(self):
                # Finish Function gets called after the last MCS
                pass

    #IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min #(Sugarbaker)
    # IVtMins=2.1996*tMins # linear fit for first 5 min (from [C]=0.0 to [C]=~11muM, t=min)#(Sugarbaker)
    #             IVtMins = 0.9731*tMins # linear fit for first 5 min (Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class DiffusionSolverFESteeringGemcitabineIV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
        print 'starting DiffusionSolverFESteeringGemcitabineIV'
    def start(self):
        pass
    def step(self,mcs):
        #print 'drug30Mins = ',drug30Mins
        #print 'aggressInfusTimesGem[2]=',aggressInfusTimesGem[2]
        if (aggressInfusTimesGem[0] < mcs < aggressInfusTimesGem[1]) or (aggressInfusTimesGem[2] < mcs < aggressInfusTimesGem[3]) or (aggressInfusTimesGem[4] < mcs < aggressInfusTimesGem[5]):  # FLOATS; USE CONDITIONALS WITHOUT "="

            # first infusion
            if aggressInfusTimesGem[0] <= mcs < drug30Mins + aggressInfusTimesGem[0]:
                tMins = (mcs - aggressInfusTimesGem[0])/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3 # linear fit infusion period of 30 mins (Fan et al., 2010)
                # print 'infusion 1 gem IVtMins=',IVtMins
            elif drug30Mins + aggressInfusTimesGem[0] <= mcs < aggressInfusTimesGem[1]: # end of infusion to end of decay
                tMins = (mcs - aggressInfusTimesGem[0])/CisGem1Min # DO NOT SUBTRACT INFUSION TIME
                IVtMins = 101.3452 * math.exp(- 0.0676 * tMins) # Fan, 2010
                # print 'decay 1 gem IVtMins=',IVtMins
                
            # second infusion
            elif aggressInfusTimesGem[2] <= mcs < drug30Mins + aggressInfusTimesGem[2]:
                tMins = (mcs - aggressInfusTimesGem[2])/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3
            elif aggressInfusTimesGem[2] + drug30Mins <= mcs < aggressInfusTimesGem[3]:
                tMins = ((mcs - aggressInfusTimesGem[2])/CisGem1Min)
                IVtMins = 101.3452 * math.exp(- 0.0676 * tMins)

            # third infusion
            elif aggressInfusTimesGem[4] <= mcs < drug30Mins + aggressInfusTimesGem[4]:
                tMins = (mcs - aggressInfusTimesGem[4])/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3
            elif aggressInfusTimesGem[4] + drug30Mins <= mcs < aggressInfusTimesGem[5]:
                tMins = ((mcs - aggressInfusTimesGem[4])/CisGem1Min)
                IVtMins = 101.3452 * math.exp(- 0.0676 * tMins)

            # update IV conc
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins             # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

    def finish(self):
        # Finish Function gets called after the last MCS
        pass




# *****************************
# CELLS ACCUMULATE CISPLATIN AND REMOVE IT FROM SURFACE OF CELL
# limit of cellular cisplatin absorption not known
class SecretionSteppableCisplatin(SecretionBasePy,SteppableBasePy):

# MEDIUM is a null pointer, not a true cell type, and cannot secrete (conversation w/CC3D developer Maciek Swat, July 2014); must add another medium-like cell type to secrete
    def __init__(self,_simulator,_frequency=1):
        SecretionBasePy.__init__(self,_simulator, _frequency)

    def start(self):
        print "This function (SecretionSteppableCisplatin) is called at every MCS"

    def step(self,mcs):
        ## START PROFILER
        # profile = cProfile.Profile()
        # profile.enable()
        if aggressInfusTimeDay1Cis[0] < mcs < aggressInfusTimeDay1Cis[1]: # at correct time in regimen
        # if mcs > aggressInfusTimeDay1Cis[0] and mcs < aggressInfusTimeDay1Cis[1]:
            for cell in self.cellList:
                if cell.type!=1 and cell.type!=2 and cell.type!=3:  # Vessel(no accum), LungNormal(below), Dead (below)
                    comPt=CompuCell.Point3D()
                    field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual)
                    comPt.x=int(cell.xCM)
                    comPt.y=int(cell.yCM)
                    comPt.z=int(cell.zCM)
                    cisplatin=field.get(comPt) # get concentration at center of mass
                    attrSecretor=self.getFieldSecretor("Cisplatin")
                    if cisplatin > 0:
                        dictionaryAttrib = CompuCell.getPyAttrib(cell)
                        # ADD EMPIRICALLY-DETERMINED FRACTION OF CURRENT CONCENTRATION AT CELL COM TO ACCUMULATED CONCENTRATION IN CELL (DICTIONARY)
                        accumC=(cisplatin * cell.dict["accumRtCis"]) #  siteConcCis * microM/MCS = (-0.8242 * IC50 + 67.2261) * siteConcCis/50 * 1/1.5E6 * 1/10^9 * 1/$B$6 * $B$9 * 10^6    =	microM cis accumulation / MCS * frac50uMCis
                        cell.dict["cisAccum"]+=accumC
                        # REMOVE ACCUMULATED DRUG FROM EXTERNAL CONCENTRATION
                        attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
                # lung, dead(filled with active phagocytes), IC50cis, and IC50gem, equal to middle-of-the-road-least-sensitive line SCRG_SW780
                if cell.type==2 or cell.type==3:
                    comPt=CompuCell.Point3D()
                    field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    comPt.x=int(cell.xCM)
                    comPt.y=int(cell.yCM)
                    comPt.z=int(cell.zCM)
                    cisplatin=field.get(comPt)
                    attrSecretor=self.getFieldSecretor("Cisplatin")
                    if cisplatin > 0:
                        dictionaryAttrib = CompuCell.getPyAttrib(cell)
                        accumC=(cisplatin * cispAccumFrac_SCRG_SW780) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                        cell.dict["cisAccum"]+=accumC
                        attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary

        # print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["cisAccum"],'microM cisplatin.'

                 #     # REMOVE EMPIRICALLY-DETERMINED EFFLUXED DRUG FROM ACCUMULATED DRUG IN CELL
                #     # IF IMPLEMENTING EFFLUX AFTER EXTERNAL DRUG = 0, ADD EFFLUX FROM TOTAL DRUG ACCUMULATED PER CELL BACK INTO EXTERNAL DRUG CONCENTRATION. Net accumulation empirical measurements include efflux.
                # #EFFLUX (REMOVE FROM DICTIONARY, ADD TO FIELD)
                # dictionaryAttrib = CompuCell.getPyAttrib(cell)
                 # # efflux=3.58E-06*dictionaryAttrib[3] # for normal tissue
                # efflux=1.38E-06*dictionaryAttrib[3] # for tumor tissue
                # if efflux>0:
                #     dictionaryAttrib[3]-=efflux
                #     attrSecretor.secreteInsideCellAtBoundary(cell,efflux) # uM secretion from pixels at outer boundary of cell




# *****************************
# CELLS ACCUMULATE GEMCITABINE AND REMOVE IT FROM SURFACE OF CELL
# limit of cellular gemcitabine absorption not known
class SecretionSteppableGemcitabine(SecretionBasePy,SteppableBasePy):
# MEDIUM is a null pointer, not a true cell type, and cannot secrete (conversation w/CC3D developer Maciek Swat, July 2014); must add another medium-like cell type to secrete
    def __init__(self,_simulator,_frequency=1):
        SecretionBasePy.__init__(self,_simulator, _frequency)

    def start(self):
        print "This function (SecretionSteppableGemcitabine) is called at every MCS"

    def step(self,mcs):
        if (aggressInfusTimesGem[0] < mcs < aggressInfusTimesGem[1]) or (aggressInfusTimesGem[2] < mcs < aggressInfusTimesGem[3]) or (aggressInfusTimesGem[4] < mcs < aggressInfusTimesGem[5]):
             for cell in self.cellList:
                 if cell.type!=1 and cell.type!=2 and cell.type!=3:  # Vessel(no accum), LungNormal(below), Dead (below)
                     comPt=CompuCell.Point3D()
                     field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                     # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual)
                     comPt.x=int(cell.xCM)
                     comPt.y=int(cell.yCM)
                     comPt.z=int(cell.zCM)
                     gemcitabine=field.get(comPt) # get concentration at center of mass
                     attrSecretor=self.getFieldSecretor("Gemcitabine")
                     if gemcitabine > 0:
                         dictionaryAttrib = CompuCell.getPyAttrib(cell)
                         # ADD EMPIRICALLY-DETERMINED FRACTION OF CURRENT CONCENTRATION AT CELL COM TO ACCUMULATED CONCENTRATION IN CELL (DICTIONARY)
                         accumG=(gemcitabine * cell.dict["accumRtGem"]) 
                         cell.dict["gemAccum"]+=accumG
                         # REMOVE ACCUMULATED DRUG FROM EXTERNAL CONCENTRATION
                         attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
                     # lung, dead(filled with active phagocytes): equal to middle-of-the-road-least-sensitive line SCRG_SW780
                     if cell.type==2 or cell.type==3:
                         comPt=CompuCell.Point3D()
                         field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                         comPt.x=int(cell.xCM)
                         comPt.y=int(cell.yCM)
                         comPt.z=int(cell.zCM)
                         gemcitabine=field.get(comPt)
                         attrSecretor=self.getFieldSecretor("Gemcitabine")
                         if gemcitabine > 0:
                             dictionaryAttrib = CompuCell.getPyAttrib(cell)
                             accumG=(gemcitabine * gemAccumFrac_SCRG_SW780) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                             cell.dict["gemAccum"]+=accumG
                             attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary

            # print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["gemAccum"],'microM gemcitabine.'




######################################### CELL TYPES CHANGE AT GEMCITABINE IC50 THRESHOLD
class ChangeAtGemIC50Steppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
          # 2 lines following are mentioned in the manual, but are not in Listing 17,
          #    or the demo (cellsort_2D_field_modules.py)
          # self.dim=self.simulator.getPotts().getCellFieldG().getDim()
          # self.fieldName="FGF"
    def start(self):
        pass
    def step(self,mcs):
        for cell in self.cellList:
            if cell.type!=0 and cell.type!=1 and cell.type!=2 and cell.type!=3: # all cell types accumulate cisplating except for Vessel, LungNormal, Dead, respectively
                # print 'inside gemAcuum: celltype=',cell.type,', cell.dict=',cell.dict
                if cell.dict["gemAccum"] > cell.dict["IC50Gem"]:
                    # if cell.type == 14, cell.type = 16, else
                    cell.type=13




######################################### CELL TYPES CHANGE AT CISPLATIN IC50 THRESHOLD
class ChangeAtCisIC50Steppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
    def start(self):
        pass
    def step(self,mcs):
        for cell in self.cellList:
            if cell.type!=0 and cell.type!=1 and cell.type!=2 and cell.type!=3: # all cell types accumulate cisplating except for Vessel, LungNormal, Dead, respectively
                # print 'inside cisAcuum: print 'celltype=',cell.type,', cell.dict=',cell.dict
                if cell.dict["cisAccum"] > cell.dict["IC50Cis"]:
                    cell.type=12




# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& DATA RECORDING FUNCTIONS
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




class CispAccumVisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField=self.createScalarFieldCellLevelPy("AccumulatedCisplFractionIC50")
    def step(self,mcs):
        self.scalarCLField.clear()
        for cell in self.cellList:
            if cell.type!=1 and cell.type!=2: # all cell types except Vessel, LungNormal, respectively
                # print cell.type, cell.id
                # print cell.dict
                self.scalarCLField[cell]=cell.dict["cisAccum"]/cell.dict["IC50Cis"]




class GemAccumVisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField=self.createScalarFieldCellLevelPy("AccumulatedGemFractionIC50")
    def step(self,mcs):
        self.scalarCLField.clear()
        for cell in self.cellList:
            # print cell.dict
            if cell.type!=1 and cell.type!=2: # all cell types except Vessel, LungNormal, respectively
                self.scalarCLField[cell]=cell.dict["gemAccum"]/cell.dict["IC50Gem"]




class PlotCellPops(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=(10)):
    def __init__(self,_simulator,_frequency=(3941)): # = MCS per min * 60 min
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        # hex color codes from http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
        self.pW=self.addNewPlotWindow(_title='Cell Populations',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Variables', _xScaleType='linear',_yScaleType='linear')
#        self.pW.addPlot('DATA_SERIES_1',_style='Dots',_color='red',_size=5)
#        self.pW.addPlot('DATA_SERIES_2',_style='Steps',_size=1)
        self.pW.addPlot("SCSG_BFTC_905_pop",_style='Dots',_color='Aquamarine',_size=5)
        self.pW.addPlot("SCSG_J82_pop",_style='Dots',_color='Green Yellow',_size=5)
        self.pW.addPlot("RCRG_RT4_pop",_style='Dots',_color='Medium Blue',_size=5)
        self.pW.addPlot("RCRG_HT_1197_pop",_style='Dots',_color='Dark Slate Blue',_size=5)
        self.pW.addPlot("SCRG_SW780_pop",_style='Dots',_color='Light Sea Green',_size=5)
        self.pW.addPlot("SCRG_KU_19_19_pop",_style='Dots',_color='Forest Green',_size=5)
        self.pW.addPlot("RCSG_LB831_BLC_pop",_style='Dots',_color='Lawn Green',_size=5)
        self.pW.addPlot("RCSG_DSH1_pop",_style='Dots',_color='Chartreuse',_size=5)
        self.pW.addPlot("IC50Cis_pop",_style='Dots',_color='Blue',_size=5)
        self.pW.addPlot("IC50Gem_pop",_style='Dots',_color='Blue Violet',_size=5)
        self.pW.addPlot("CisResist_pop",_style='Dots',_color='Light Sky Blue',_size=5)
        self.pW.addPlot("GemResist_pop",_style='Dots',_color='Orchid',_size=5)
        self.pW.addPlot("DualResist_pop",_style='Dots',_color='Lavender',_size=5)
        self.pW.addPlot("LungNormal_pop",_style='Dots',_color='Deep Pink',_size=5)
        self.pW.addPlot("Vessel_pop",_style='Dots',_color='red',_size=5)
        self.pW.addPlot("Dead_pop",_style='Dots',_color='Light Slate Gray',_size=5)

    def step(self, mcs):

        # TypeId="4" TypeName="SCSG_BFTC_905"
        # TypeId="5" TypeName="SCSG_J82"
        # TypeId="6" TypeName="RCRG_RT4"
        # TypeId="7" TypeName="RCRG_HT_1197"
        # TypeId="8" TypeName="SCRG_SW780"
        # TypeId="9" TypeName="SCRG_KU_19_19"
        # TypeId="10" TypeName="RCSG_LB831_BLC"
        # TypeId="11" TypeName="RCSG_DSH1"
        # TypeId="12" TypeName="IC50Cis"
        # TypeId="13" TypeName="IC50Gem"
        # # initialize in case cell type isn't yet present
        # SCSG_BFTC_905_pop = 0
        # SCSG_J82_pop = 0
        # RCRG_RT4_pop = 0
        # RCRG_HT_1197_pop = 0
        # SCRG_SW780_pop = 0
        # SCRG_KU_19_19_pop = 0
        # RCSG_LB831_BLC_pop = 0
        # RCSG_DSH1_pop = 0
        # IC50Cis_pop = 0
        # IC50Gem_pop = 0
        # LungNormal_pop = 0
        # Vessel_pop = 0
        # Dead_pop = 0

        SCSG_BFTC_905_pop = float(len(self.cellListByType(self.SCSG_BFTC_905)))
        SCSG_J82_pop = float(len(self.cellListByType(self.SCSG_J82)))
        RCRG_RT4_pop = float(len(self.cellListByType(self.RCRG_RT4)))
        RCRG_HT_1197_pop = float(len(self.cellListByType(self.RCRG_HT_1197)))
        SCRG_SW780_pop = float(len(self.cellListByType(self.SCRG_SW780)))
        SCRG_KU_19_19_pop = float(len(self.cellListByType(self.SCRG_KU_19_19)))
        RCSG_LB831_BLC_pop = float(len(self.cellListByType(self.RCSG_LB831_BLC)))
        RCSG_DSH1_pop = float(len(self.cellListByType(self.RCSG_DSH1)))
        IC50Cis_pop = float(len(self.cellListByType(self.IC50CIS)))
        IC50Gem_pop = float(len(self.cellListByType(self.IC50GEM)))
        CisResist_pop = float(len(self.cellListByType(self.CISRESIST)))
        GemResist_pop = float(len(self.cellListByType(self.GEMRESIST)))
        DualResist_pop = float(len(self.cellListByType(self.DUALRESIST)))
        LungNormal_pop = float(len(self.cellListByType(self.LUNGNORMAL)))
        Vessel_pop = float(len(self.cellListByType(self.VESSEL)))
        Dead_pop = float(len(self.cellListByType(self.DEAD)))

        hrs=mcs/CisGem1Min/60.0
        days=hrs/24.0
        self.pW.addDataPoint("SCSG_BFTC_905_pop",days,SCSG_BFTC_905_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCSG_J82_pop",days,SCSG_J82_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCRG_RT4_pop",days,RCRG_RT4_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCRG_HT_1197_pop",days,RCRG_HT_1197_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCRG_SW780_pop",days,SCRG_SW780_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCRG_KU_19_19_pop",days,SCRG_KU_19_19_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCSG_LB831_BLC_pop",days,RCSG_LB831_BLC_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCSG_DSH1_pop",days,RCSG_DSH1_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("IC50Cis_pop",days,IC50Cis_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("IC50Gem_pop",days,IC50Gem_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("CisResist_pop",days,CisResist_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("GemResist_pop",days,GemResist_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("DualResist_pop",days,DualResist_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("LungNormal_pop",days,LungNormal_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("Vessel_pop",days,Vessel_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("Dead_pop",days,Dead_pop) # arguments are (name of the data series, x, y)

        # self.pW.eraseAllData()

        self.pW.savePlotAsPNG('CellPops.png',1000,1000) # here we specify size of the image saved (1000x1000) - default is 400 x 400
        self.pW.savePlotAsData('CellPops.txt')

    def finish(self):
        pass




class PlotDrugs(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=(10)):
    def __init__(self,_simulator,_frequency=(3941)): # = MCS per min * 60 min
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        # hex color codes from http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
        self.pW=self.addNewPlotWindow(_title='Cell Populations',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Variables', _xScaleType='linear',_yScaleType='linear')
#        self.pW.addPlot('DATA_SERIES_1',_style='Dots',_color='red',_size=5)
#        self.pW.addPlot('DATA_SERIES_2',_style='Steps',_size=1)
        self.pW.addPlot("SCSG_BFTC_905_pop",_style='Dots',_color='Aquamarine',_size=5)
        self.pW.addPlot("SCSG_J82_pop",_style='Dots',_color='Green Yellow',_size=5)
        self.pW.addPlot("RCRG_RT4_pop",_style='Dots',_color='Medium Blue',_size=5)
        self.pW.addPlot("RCRG_HT_1197_pop",_style='Dots',_color='Dark Slate Blue',_size=5)
        self.pW.addPlot("SCRG_SW780_pop",_style='Dots',_color='Light Sea Green',_size=5)
        self.pW.addPlot("SCRG_KU_19_19_pop",_style='Dots',_color='Forest Green',_size=5)
        self.pW.addPlot("RCSG_LB831_BLC_pop",_style='Dots',_color='Lawn Green',_size=5)
        self.pW.addPlot("RCSG_DSH1_pop",_style='Dots',_color='Chartreuse',_size=5)
        self.pW.addPlot("IC50Cis_pop",_style='Dots',_color='Blue',_size=5)
        self.pW.addPlot("IC50Gem_pop",_style='Dots',_color='Blue Violet',_size=5)
        self.pW.addPlot("LungNormal_pop",_style='Dots',_color='Deep Pink',_size=5)
        self.pW.addPlot("Vessel_pop",_style='Dots',_color='red',_size=5)
        self.pW.addPlot("Dead_pop",_style='Dots',_color='Light Slate Gray',_size=5)

    def step(self, mcs):

        # TypeId="4" TypeName="SCSG_BFTC_905"
        # TypeId="5" TypeName="SCSG_J82"
        # TypeId="6" TypeName="RCRG_RT4"
        # TypeId="7" TypeName="RCRG_HT_1197"
        # TypeId="8" TypeName="SCRG_SW780"
        # TypeId="9" TypeName="SCRG_KU_19_19"
        # TypeId="10" TypeName="RCSG_LB831_BLC"
        # TypeId="11" TypeName="RCSG_DSH1"
        # TypeId="12" TypeName="IC50Cis"
        # TypeId="13" TypeName="IC50Gem"
        # # initialize in case cell type isn't yet present
        # SCSG_BFTC_905_pop = 0
        # SCSG_J82_pop = 0
        # RCRG_RT4_pop = 0
        # RCRG_HT_1197_pop = 0
        # SCRG_SW780_pop = 0
        # SCRG_KU_19_19_pop = 0
        # RCSG_LB831_BLC_pop = 0
        # RCSG_DSH1_pop = 0
        # IC50Cis_pop = 0
        # IC50Gem_pop = 0
        # LungNormal_pop = 0
        # Vessel_pop = 0
        # Dead_pop = 0

        SCSG_BFTC_905_pop = float(len(self.cellListByType(self.SCSG_BFTC_905)))
        SCSG_J82_pop = float(len(self.cellListByType(self.SCSG_J82)))
        RCRG_RT4_pop = float(len(self.cellListByType(self.RCRG_RT4)))
        RCRG_HT_1197_pop = float(len(self.cellListByType(self.RCRG_HT_1197)))
        SCRG_SW780_pop = float(len(self.cellListByType(self.SCRG_SW780)))
        SCRG_KU_19_19_pop = float(len(self.cellListByType(self.SCRG_KU_19_19)))
        RCSG_LB831_BLC_pop = float(len(self.cellListByType(self.RCSG_LB831_BLC)))
        RCSG_DSH1_pop = float(len(self.cellListByType(self.RCSG_DSH1)))
        IC50Cis_pop = float(len(self.cellListByType(self.IC50CIS)))
        IC50Gem_pop = float(len(self.cellListByType(self.IC50GEM)))
        LungNormal_pop = float(len(self.cellListByType(self.LUNGNORMAL)))
        Vessel_pop = float(len(self.cellListByType(self.VESSEL)))
        Dead_pop = float(len(self.cellListByType(self.DEAD)))

        ##### CIS DRUG CONCENTRATIONS AFTER IV DELIVERY:
        if aggressInfusTimeDay1Cis[0] < mcs < aggressInfusTimeDay1Cis[1]: # at correct time in regimen

            # infusion
            if aggressInfusTimeDay1Cis[0] < mcs < aggressInfusTimeDay1Cis[0] + drug15Mins:
                tMins= (mcs - aggressInfusTimeDay1Cis[0]) / CisGem1Min # time since injection
                IVtMins = 0.3725*tMins # linear fit for 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
                print 'infusion cis IVtMins=',IVtMins
            elif aggressInfusTimeDay1Cis[0] + drug15Mins <= mcs < aggressInfusTimeDay1Cis[0] + cIVFirstPoint: # plateau for 5.7m
                IVtMins = 5.59 # constant for ~6 mins mins; highest and first data point after infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
                print 'plateau cis IVtMins=',IVtMins
            elif aggressInfusTimeDay1Cis[0] + cIVFirstPoint <= mcs < aggressInfusTimeDay1Cis[1]:
            # aggressInfusTimeDay1Cis[0] + cEndDataSet:        # prior to end of IV data set
                tMins=((mcs - aggressInfusTimeDay1Cis[0])/CisGem1Min) - (5.742+15) # take away added infusion and plateau time so fit is correct; use floats
                IVtMins = -1.154e-06*tMins**3 + 0.0005737*tMins**2 - 0.09922*tMins + 5.973 # Casper, 1984
                if IVtMins < 0:
                    IVtMins=0 # in case time frame goes past where fit becomes negative
                print 'decay cis tMins=',tMins
                print 'decay cis IVtMins= ',IVtMins


            # GEM first infusion
            if mcs < drug30Mins:
                tMins = mcs/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3 # linear fit infusion period of 30 mins (Fan et al., 2010)
                # print 'infusion 1 gem IVtMins=',IVtMins
            elif drug30Mins <= mcs < aggressInfusTimesGem[1]: # end of infusion to end of decay
                tMins=(mcs/CisGem1Min) - 30.0 # take away infusion time so fit is correct, starting at t = 0; use floats
                IVtMins =101.3452 * math.exp(- 0.0676 * tMins) # Fan, 2010
                # print 'decay 1 gem IVtMins=',IVtMins
            # second infusion
            elif aggressInfusTimesGem[2] <= mcs < drug30Mins + aggressInfusTimesGem[2]:
                tMins = (mcs - aggressInfusTimesGem[2])/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3
            elif aggressInfusTimesGem[2] + drug30Mins <= mcs < aggressInfusTimesGem[3]:
                tMins=((mcs - aggressInfusTimesGem[2])/CisGem1Min) - 30.0
                IVtMins =101.3452 * math.exp(- 0.0676 * tMins)
            # third infusion
            elif aggressInfusTimesGem[4] <= mcs < drug30Mins + aggressInfusTimesGem[4]:
                tMins = (mcs - aggressInfusTimesGem[2])/CisGem1Min
                IVtMins = 6.8*(tMins/15 - 1) + 7.3
            elif aggressInfusTimesGem[4] + drug30Mins <= mcs < aggressInfusTimesGem[5]:
                tMins=((mcs - aggressInfusTimesGem[2])/CisGem1Min) - 30.0
                IVtMins =101.3452 * math.exp(- 0.0676 * tMins)

        hrs=mcs/CisGem1Min/60.0
        days=hrs/24.0
        self.pW.addDataPoint("SCSG_BFTC_905_pop",days,SCSG_BFTC_905_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCSG_J82_pop",days,SCSG_J82_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCRG_RT4_pop",days,RCRG_RT4_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCRG_HT_1197_pop",days,RCRG_HT_1197_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCRG_SW780_pop",days,SCRG_SW780_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("SCRG_KU_19_19_pop",days,SCRG_KU_19_19_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCSG_LB831_BLC_pop",days,RCSG_LB831_BLC_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("RCSG_DSH1_pop",days,RCSG_DSH1_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("IC50Cis_pop",days,IC50Cis_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("IC50Gem_pop",days,IC50Gem_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("LungNormal_pop",days,LungNormal_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("Vessel_pop",days,Vessel_pop) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("Dead_pop",days,Dead_pop) # arguments are (name of the data series, x, y)

        # self.pW.eraseAllData()

        self.pW.savePlotAsPNG('CellPops.png',1000,1000) # here we specify size of the image saved (1000x1000) - default is 400 x 400
        self.pW.savePlotAsData('CellPops.txt')

    def finish(self):
        pass



    
class PrintCellData(SteppableBasePy):
    def __init__(self,_simulator,_frequency=(3941)):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        #        self.simulator=_simulator
        #        self.inventory=self.simulator.getPotts().getCellInventory()
        #        self.cellList=CellList(self.inventory)
        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=('CellData_' + datestring)
        #        self.openFileInSimulationOutputDirectory('%s.txt' %fileName,'w')
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        # self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        self.file.write('MCS CellType CellID AgeHrs cycleHrs generation numDivisions HrsSinceDeath cisAccum gemAccum cisResistance gemResistance IC50Cis IC50Gem accumRtCis accumRtGem cell.lambdaVolume cell.targetVolume\n')
    def step(self, mcs):
        for cell in self.cellList:
            print cell.type,            cell.id,            cell.dict["AgeHrs"],            cell.dict["cycleHrs"],            cell.dict["generation"],            cell.dict["numDivisions"],            cell.dict["HrsSinceDeath"],            cell.dict["cisAccum"],            cell.dict["gemAccum"],            cell.dict["cisResistance"],            cell.dict["gemResistance"],            cell.dict["IC50Cis"],            cell.dict["IC50Gem"],            cell.dict["accumRtCis"],            cell.dict["accumRtGem"],            cell.lambdaVolume,            cell.targetVolume
            self.file.write("%d %d %d %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f \n" %(mcs,
                                                                       cell.type,
                                                                       cell.id,
                                                                       cell.dict["AgeHrs"],
                                                                       cell.dict["cycleHrs"],
                                                                       cell.dict["generation"],
                                                                       cell.dict["numDivisions"],
                                                                       cell.dict["HrsSinceDeath"],
                                                                       cell.dict["cisAccum"],
                                                                       cell.dict["gemAccum"],
                                                                       cell.dict["cisResistance"],
                                                                       cell.dict["gemResistance"],
                                                                       cell.dict["IC50Cis"],
                                                                       cell.dict["IC50Gem"],
                                                                       cell.dict["accumRtCis"],
                                                                       cell.dict["accumRtGem"],
                                                                       cell.lambdaVolume,
                                                                       cell.targetVolume))
                                                                                      #            self.file.write("%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n" %(mcs,cell.type,cell.id
                                                                                      #,,,,,,))

    def finish(self):
        self.file.close() # close the file
