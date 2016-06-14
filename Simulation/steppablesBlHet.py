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

""" "TO GET ALL CELL ATTRIBUTES" (to see what cell attributes can be accessed/changed in Python):
print dir(cell)
output 7-2-2012:
['__class__', '__del__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattr__', '__getattribute__', '__hash__', '__init__', '__module__',
'__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__swig_destroy__', '__swig_getmetho
ds__', '__swig_setmethods__', '__weakref__', 'angle', 'averageConcentration', 'clusterId', 'clusterSurface', 'ecc', 'extraAttribPtr', 'flag', 'fluctAm
pl', 'iXX', 'iXY', 'iXZ', 'iYY', 'iYZ', 'iZZ', 'id', 'lX', 'lY', 'lZ', 'lambdaClusterSurface', 'lambdaSurface', 'lambdaVecX', 'lambdaVecY', 'lambdaVec
Z', 'lambdaVolume', 'pyAttrib', 'subtype', 'surface', 'targetClusterSurface', 'targetSurface', 'targetVolume', 'this', 'type', 'volume', 'xCM', 'xCOM'
, 'xCOMPrev', 'yCM', 'yCOM', 'yCOMPrev', 'zCM', 'zCOM', 'zCOMPrev']
"""


#STRINGS FOR OUTPUT FILES
drug='cisplatin'
#size='lgSphere'
#vess='451v_0.99pctA'
#vess='151v_2.17pctA'
# size='smSphereIV'
# size='smSphereIP'
# size='smSphereIPIV'
# size='smSphereIV'
# size='smSphereIVinfus'
# size='smSphereIPinfus'
# size='R140CubeEdgeIVinfus'
# size='R140CubeEdgeIPinfus'
size='R140CubeCenterIVinfus'
# size='R140CubeCenterIPinfus'
# size='smSphereIPIV'
# vess='18v_2.26pctA'
# vess='17v_2.26pctA'
# vess='16v_7.64pctA'
# vess='11v_10.89pctA'
# vess='11v_6.65pctA'
# vess='11v_4.53pctA'
# vess='10v_2.26pctA'
# vess='8v_4.38pctA'
# vess='6v_3.82pctA'
# vess='3v_10.04pctA'
# vess='2v_1.27pctA'

# vess='1v_3pctvol_2pctA'
# vess='7v_5pctvol_4pctA'
# vess='6v_8pctvol_6pctA'
# vess='15v_10pctvol_8pctA'
# vess='19v_12pctvol_10pctA'
# vess='0v_0.0pctA'
vess='_10pctA' # for cubes



#INITIALIZE COUNTING VARIABLES
numPCancer=0
numIC50Cancer=0


#INITIALIZE CHEMICAL THRESHOLDS

# FINAL CISPLATIN IC50 = MISTRY 1992
cisplatinIC50=52.71 # FINAL USED: muM (equiv to (equitoxic) 2h-IC-50 Mistry, 1992
# cisplatinIC50=38.3 # (MATCHES 2H DRUG TIME COURSE): muM (SD=12.6, in SKOV-3, 2h exposure; Table 1, Mistry, 1992)
# cisplatinIC50=0.0001 #test
# cisplatinIC50=126.63 # muM (+/- 12.06 micromols/L, in SKOV3ip1, 48h exposure; Fig. 3, Xu, 2008)
# cisplatinIC50  =3.33 # muM (Nakaro 1997, minimum effective concentration for DNA damage)

# DIFFUSION COEFFICIENTS
CisGem1Min = 65.678 # 65.678 mcs = 1 min of diffusion time for cells of diameter T24 bladder cancer cell line, for drugs with diffusion coeff. of sodium fluorescein

#IV CISPLATIN
Infusion15Mins=15*CisGem1Min # 18107.745 MCS; subtract from runtime once begin using concentration time course
cIVFirstPoint=(5.742+15)*CisGem1Min # 6931.79 MCS; five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
cEndDataSet=(170.862*CisGem1Min)+Infusion15Mins
#### first five minutes (no data, linear fit from t=0, [iv]=0 to first time point
# cFirst5Mins=5*465.189 # five minutes worth of MCSes in normal tissue (Swabb 1974?)
cFirst5Mins=5*CisGem1Min # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
####
cFirst20Mins=20*CisGem1Min # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
# cIVZero=math.exp(-16.094/-3.338)*1207.183 ##mcs mins, MCS scaled for diffusion in tumor tissue;Excel check, yintercept = 124.1449668
# x-intercept of cisplatin IP-post-IV concentration function
# cIPpostIV_zero=(0.386/0.5431)*465.... #mcs hrs, MCS scaled for diffusion in normal tissue
# cIPpostIVZero= (0.386/0.0091)*1207.183 #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after this time point (42.42mins = 51205.8mcs)#  cIPpostIVZero= 896.305 #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after t = 896.305 (Wolfram Alpha cubic polynomial solver).  Will not use because becomes zero after IV goes to zero; then soln. is zero
#  cIVpostPZero=~ 7hrs #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after t = 7 hrs

#ip CISPLATIN
# cIPFirstPt=0.2126*60*1207.183 #hrs, first time point
# cIPFirstPt=7.386*1207.183 #7.386 min, first time point
cIPFirstPtPlus15=(7.386+15.0)*CisGem1Min # = mcs to get to 7.386 min, first time point in IV infusion, plus 15 mins infusion time; initialize with floats
#cIVFirstPtPlus15=(5.742+15.0)*1207.183 # = 25039.389786 MCS = 5.742 min, first time point, plus 15 mins infusion time

# IV GEMCITABINE
gem30Mins=30*CisGem1Min

## CELL IDs
# TypeId="4" TypeName="SCSG_BFTC_# 905"
# TypeId="5" TypeName="SCSG_J82"
# TypeId="6" TypeName="RCRG_RT4"
# TypeId="7" TypeName="RCRG_HT_1197"
# TypeId="8" TypeName="SCRG_SW780"
# TypeId="9" TypeName="SCRG_KU_19_19"
# TypeId="10" TypeName="RCSG_LB831_BLC"
# TypeId="11" TypeName="RCSG_DSH1"
# TypeId="12" TypeName="IC50Cis"
# TypeId="13" TypeName="IC50Gem"

## cisplatin, platinum accumulation per cell per time step, based on IC50 of bladder cancer cell line;
## also = concentration removed from voxel; microM/MCS * siteConcCis(microM)
cispAccumFrac_SCSG_BFTC_905 = 0.3147446166  # (sens cis and gem)	2.575477619	IC50 microM	cisplatin
cispAccumFrac_SCSG_J82 = 0.3033715255	    # (sens cis and gem)	5.42972235	IC50 microM	cisplatin				
cispAccumFrac_RCRG_RT4 = 0.2154448963	    # (resist cis and gem)	27.49620513	IC50 microM	cisplatin				
cispAccumFrac_RCRG_HT_1197 = 0.1498013738   # (resist cis and gem)	43.97041406	IC50 microM	cisplatin				
cispAccumFrac_SCRG_SW780 = 0.2691142767     # (sens cis resist gem)	14.02708355	IC50 microM	cisplatin					
cispAccumFrac_SCRG_KU_19_19 = 0.2847440491  # (sens cis resist gem)	10.10456195	IC50 microM	cisplatin
cispAccumFrac_RCSG_LB831_BLC = 0.02925373453# (resist cis sens gem)	225.1619678	IC50 microM	cisplatin
cispAccumFrac_RCSG_DSH1 = 0.04556319402     # (resist cis sens gem)	144.5646771	IC50 microM	cisplatin				

## gemcitabine (possibly dFdCtP) accumulation per cell per time step, based on IC50 of bladder cancer cell line;
## also = concentration removed from voxel; microM/MCS * siteConcGem(microM)
gemAccumFrac_SCSG_BFTC_905 = 4.41575E-03  # (sensitive cis and gem)	5.15E-06	IC50 microM	gemcitabine			
gemAccumFrac_SCSG_J82 = 4.41565E-03       # (sensitive cis and gem)
gemIC50_SCSG_J82 = 0.007409799	# IC50 microM	gemcitabine				
gemAccumFrac_RCRG_RT4 = 4.22858E-03       # (resistant cis and gem)	13.84278281	IC50 microM	gemcitabine				
gemAccumFrac_RCRG_HT_1197 = 4.39190E-03   # (resistant cis and gem)	1.764248816	IC50 microM	gemcitabine				
gemAccumFrac_SCRG_SW780 = 2.68443E-03     # (sens cis resist gem)	128.0487435	IC50 microM	gemcitabine				
gemAccumFrac_SCRG_KU_19_19 = 4.18843E-03  # (sens cis resist gem)	16.81235445	IC50 microM	gemcitabine				
gemAccumFrac_RCSG_LB831_BLC = 4.41518E-03 # (resist cis sens gem)	0.041854289	IC50 microM	gemcitabine				
gemAccumFrac_RCSG_DSH1 = 4.41445E-03      # (resist cis sens gem)	0.096498675	IC50 microM	gemcitabine				


## CELL PARAMETERS
T24BCCellVol = 1 # bladder cancer cell volume (units = voxels)
MCSFractionOfHour = 0.0002537615293 # hours per MCS, based on diffusion time for one T24 cell diameter of sodium fluorescein, proxy for cisplatin and gemcitabine

#divisionCycleTimeHrs = 30 # average time to division / replication from several cancer cell lines in vitro
divisionCycleTimeHrs = 0.005 # TEST average time to division / replication from several cancer cell lines in vitro
phagocytosisEndTime = 24 # dead cells removed at 24 hours


# PRINT SIMULATION START TIME
now = datetime.now()
print "SIMULATION START TIME =",now




class CispIC50VisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField=self.createScalarFieldCellLevelPy("NetAccumulatedCispl")
    def step(self,mcs):
        self.scalarCLField.clear()
        for cell in self.cellList:
            if cell.type==4 or cell.type==5 or cell.type==6 or cell.type==7 or cell.type==8 or cell.type==9 or cell.type==10 or cell.type==11:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                self.scalarCLField[cell]=dictionaryAttrib[3]




# SETCELLCLOCKS
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

            # SAME CLOCK IN ALL CELLS
            # ALSO CHANGE IN MITOSISSTEPPABLE CLASS (~LINE 361?)
            # x = gauss(30,1)
            # y = uniform(0,30) # age of cells initialized into simulation
            # cell.dict["AgeHrs"]=y
            # cell.dict["HrsSinceDeath"]=0

            # # test non-dividing and dead cells
            cell.dict["AgeHrs"]=29.99
            cell.dict["HrsSinceDeath"]=23.99
            cell.dict["cisAccum"]=0
            cell.dict["gemAccum"]=0

            # for cell in self.cellList:
            print 'cell.id=',cell.id,' dict=',cell.dict




# *****************************
# INCREMENT AGE AND TIME SINCE DEATH
class IncrementClocks(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (IncrementClocks) is called at every MCS"

    def step(self,mcs):
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:
            cell.dict["AgeHrs"]+= MCSFractionOfHour
            if cell.type==3:
                cell.dict["HrsSinceDeath"]+= MCSFractionOfHour
        for cell in self.cellList:
            # if cell.id > 125:
                print 'cell.id=',cell.id,'cell.type=',cell.type,' dict=',cell.dict, 'vol=',cell.targetVolume,'volLambda=',cell.lambdaVolume




# ***************************** ALSO UNCOMMENT CELL DICTIONARY KEYS 1 AND 2 FOR CELL AGE ABOVE
# VOLUMEPARAMSTEPPABLE
# SETS CELL TARGET VOLUMES OR COMPRESSIBILITIES
# VOLUME CONSTRAINTS
class VolumeParamSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (VolumeParamSteppable) is called once before simulation"
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:
            # if cell.type==4 or cell.type==5 or cell.type==6 or cell.type==7 or cell.type==8 or cell.type==9 or cell.type==10 or cell.type==11:
                cell.targetVolume=T24BCCellVol
                cell.lambdaVolume=100.0
                print 'cell.type=',cell.type,'cell.id=',cell.id,'cell.volume=',cell.targetVolume,'cell.lambdaVolume=',cell.lambdaVolume

    # def step(self,mcs):
    #     for cell in self.cellList:
    #         print "MCS",mcs,'cell.type=',cell.type,'cell.id=',cell.id,'cell.volume=',cell.targetVolume,'cell.lambdaVolume=',cell.lambdaVolume




# *****************************
# GROW CELLS WHEN READY TO DIVIDE
class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            if cell.dict["AgeHrs"]>divisionCycleTimeHrs:
                cell.targetVolume=2*T24BCCellVol
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
        print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:

            # print 'cell.id=',cell.id,' dict=',cell.dict

            if cell.dict["AgeHrs"]>divisionCycleTimeHrs:
                if cell.volume==2*T24BCCellVol:
                    cells_to_divide.append(cell)
                    print 'cell is dividing at AgeHrs',cell.dict["AgeHrs"]
                
        for cell in cells_to_divide:
            if cell.type==12 or cell.type==13:  # if cells are IC50Cis or IC50Gem
                deathChance = uniform(0,1)
                print 'deathChance=',deathChance
                if deathChance<=0.5:
                    cell.type=3 # cell dies with 50% chance
            elif cell.type!=3 and cell.type!=1 and cell.type!=2: # all cell types divide except for Vessel, LungNormal, Dead, IC50Cis, and IC50Gem
                # to change mitosis mode leave one of the below lines uncommented
                self.divideCellRandomOrientation(cell)
                # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
                # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
                # self.divideCellAlongMinorAxis(cell)                               # this is a valid option

    def updateAttributes(self):
        self.parentCell.targetVolume /= 2.0 # reduce parent target volume by increasing; = ratio to parent vol
        self.cloneParent2Child()
        # set parent lambda volume post-division
        self.parentCell.lambdaVolume = 2000 # make sure parent stays in place
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )

        self.childCell.type=self.parentCell.type
        self.childCell.dict["AgeHrs"]=0
        self.childCell.dict["HrsSinceDeath"]=0
        self.childCell.dict["cisAccum"]=0
        self.childCell.dict["gemAccum"]=0
        ## for cell in self.cellList:
        print 'childCell.id=',self.childCell.id,' dict=',self.childCell.dict,'childCell.targetVolume=', self.childCell.targetVolume,'childCell.lambdaVolume=', self.childCell.lambdaVolume

        if self.parentCell.type==5:
            self.childCell.type=6
        elif self.parentCell.type==6:
            self.childCell.type=7
            # else:
        #     self.childCell.type=1


        
        
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
                print 'I am dead cell.id', cell.id, 'and I died',cell.dict["HrsSinceDeath"],'hrs ago.'
                if cell.dict["HrsSinceDeath"]>=phagocytosisEndTime:
                    print 'removing dead cell.id', cell.id,'with cell volume',cell.volume,'cell.targetVolume',cell.targetVolume,'and cell.lambdaVolume',cell.lambdaVolume
                    cell.targetVolume=0
                    cell.lambdaVolume=1000
                    



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

        for cell in self.cellList:
            if cell.type==4:  # SCSG_BFTC_905
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
                    accumC=(cisplatin * cispAccumFrac_SCSG_BFTC_905) #  microM/MCS * siteConcCis = (-0.8242 * IC50 + 67.2261) * siteConcCis/50 * 1/1.5E6 * 1/10^9 * 1/$B$6 * $B$9 * 10^6    =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC

                    print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["cisAccum"],'microM cisplatin.'

                    # REMOVE ACCUMULATED DRUG FROM EXTERNAL CONCENTRATION
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==5: # SCSG_J82
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_SCSG_J82) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==6: # RCRG_RT4
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_RCRG_RT4) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==7: # RCRG_HT_1197
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_RCRG_HT_1197) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==8: # SCRG_SW780
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
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==9: # SCRG_KU_19_19
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_SCRG_KU_19_19) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==10: # RCSG_LB831_BLC
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_RCSG_LB831_BLC) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==11: # RCSG_DSH1
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Cisplatin")
                if cisplatin > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumC=(cisplatin * cispAccumFrac_RCSG_DSH1) #  microM/MCS * siteConcCis  =	microM cis accumulation / MCS * frac50uMCis
                    cell.dict["cisAccum"]+=accumC
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell
                    
            # lung, dead(filled with active phagocytes), IC50cis, and IC50gem, equal to middle-of-the-road-least-sensitive line SCRG_SW780
            if cell.type==2 or cell.type==3 or cell.type==12 or cell.type==13: 
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

                    print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["cisAccum"],'microM cisplatin.'




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
        ## START PROFILER
        # profile = cProfile.Profile()
        # profile.enable()

        for cell in self.cellList:
            if cell.type==4:  # SCSG_BFTC_905
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
                    accumG=(gemcitabine * gemAccumFrac_SCSG_BFTC_905) #  microM/MCS * siteConcGem = (-0.8242 * IC50 + 67.2261) * siteConcGem/50 * 1/1.5E6 * 1/10^9 * 1/$B$6 * $B$9 * 10^6    =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG

                    print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["gemAccum"],'microM gemcitabine.'

                    # REMOVE ACCUMULATED DRUG FROM EXTERNAL CONCENTRATION
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==5: # SCSG_J82
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_SCSG_J82) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==6: # RCRG_RT4
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_RCRG_RT4) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==7: # RCRG_HT_1197
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_RCRG_HT_1197) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==8: # SCRG_SW780
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
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==9: # SCRG_KU_19_19
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_SCRG_KU_19_19) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==10: # RCSG_LB831_BLC
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_RCSG_LB831_BLC) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
            if cell.type==11: # RCSG_DSH1
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Gemcitabine")
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                gemcitabine=field.get(comPt)
                attrSecretor=self.getFieldSecretor("Gemcitabine")
                if gemcitabine > 0:
                    dictionaryAttrib = CompuCell.getPyAttrib(cell)
                    accumG=(gemcitabine * gemAccumFrac_RCSG_DSH1) #  microM/MCS * siteConcGem  =	microM gem accumulation / MCS * frac50uMGem
                    cell.dict["gemAccum"]+=accumG
                    attrSecretor.uptakeInsideCellAtCOM(cell,accumG,1.0) # uM secretion from pixels at outer boundary of cell
                    
            # lung, dead(filled with active phagocytes), IC50gem, and IC50gem, equal to middle-of-the-road-least-sensitive line SCRG_SW780
            if cell.type==2 or cell.type==3 or cell.type==12 or cell.type==13: 
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

                    print "I am cell.id",cell.id,'cell.type',cell.type,'and I have accumulated',cell.dict["gemAccum"],'microM gemcitabine.'




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class DiffusionSolverFESteeringCisplatinIV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
    def start(self):
        pass
    def step(self,mcs):

        ##### DRUG CONCENTRATIONS AFTER IV DELIVERY:
        # INTRAVENOUS DRUG CONCENTRATION
        # tMins=mcs/465.189 # diffusion time for one cell diameter in normal tissue
#         if 0<=mcs<cFirst5Mins:
        if 0<=mcs<Infusion15Mins:
            tMins=mcs/CisGem1Min
            IVtMins = 0.3725*tMins # linear fit for 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

        elif Infusion15Mins<=mcs<cIVFirstPoint:
            IVtMins = 5.59 # highest and first data point after infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

        elif cIVFirstPoint<=mcs<cEndDataSet:        # prior to end of IV data se
            tMins=(mcs/CisGem1Min) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats
            IVtMins = -1.154e-06*tMins**3 + 0.0005737*tMins**2 - 0.09922*tMins + 5.973 # Casper, 1984
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins             # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
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
    def start(self):
        pass
    def step(self,mcs):

        ##### DRUG CONCENTRATIONS AFTER IV DELIVERY:
        # INTRAVENOUS DRUG CONCENTRATION
        # tMins=mcs/465.189 # diffusion time for one cell diameter in normal tissue
#         if 0<=mcs<cFirst5Mins:
        if 0<=mcs<gem30Mins:
            tMins=mcs/CisGem1Min
            IVtMins = 6.8*(tMins/15 - 1) + 7.3 # linear fit infusion period of 30 mins (Fan et al., 2010)
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

        elif gem30Mins<=mcs:
            tMins=(mcs/CisGem1Min) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats
            IVtMins =101.3452 * math.exp(- 0.0676 * tMins) # Fan, 2010
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel']))
            IVxml=IVtMins             # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Gemcitabine'],['SecretionData'],['ConstantConcentration','Type','Vessel'])
            self.updateXML()

    def finish(self):
        # Finish Function gets called after the last MCS
        pass




#########################################CELL TYPES CHANGE AT DRUG IC50 THRESHOLD
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
            if cell.type==5:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                if (cell.dict["gemAccum"]>=gemIC50_SCSG_J82):
                    cell.type=13 # IC50Gem




# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& DATA RECORDING FUNCTIONS
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PrintAllCells(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        pass
    def step(self,mcs):
        #type here the code that will run every _frequency MCS
        for cell in self.cellList:
            print "cell.type=", cell.type, ";","cell.id=",cell.id
    def finish(self):
        # Finish Function gets called after the last MCS
        pass



# INFORPRINTERSTEPPABLE
class InfoPrinterSteppable(SteppablePy):
    def __init__(self,_simulator,_frequency=60):
        SteppablePy.__init__(self,_frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (InfoPrinterSteppable) is called once before simulation"
        # numPCancer=0
        # numIC50Cancer=0
        # totalVolumeECM=0
        # totalVolumeMesothelial=0
        # totalVolumeEndothelial=0
    def step(self,mcs):
#        print "This function is called every 60 MCS"
        for cell in self.cellList:
#            dictionaryAttrib = CompuCell.getPyAttrib(cell)
            if cell.type==2:
                numPCancer+=1
#                dictionaryAttrib = CompuCell.getPyAttrib(cell)
#                print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            if cell.type==8:
                numIC50Cancer+=1
#                dictionaryAttrib = CompuCell.getPyAttrib(cell)
#                print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            # if cell.type==1:
            #     totalVolumeMesothelial+=cell.volume
            # if cell.type==11:
                # totalVolumeEndothelial+=cell.volume
 #               dictionaryAttrib = CompuCell.getPyAttrib(cell)
 #               print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            # if cell.type==9:
            #     totalVolumeECM+=cell.volume
        # print "TOTAL # PCancer cells=",numPCancer
        # print "total ECM vol = ",totalVolumeECM
        # print "total Mesothelium vol = ",totalVolumeMesothelial
        # print "total new capillary (Endothelial) vol = ",totalVolumeEndothelial

            # file.write("CELL ID=%d CELL TYPE=%d volume=%d\n" %(cell.id,cell.type,cell.volume))




class CisplatinToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=1207): # one minutes of diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=24240): # twenty minutes of diffusion, 1 cell diameter/MCS in tumor tissue
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring+'_'+ drug+'_' +size+'_' +vess+'_' +'Concentrations_')
#        self.openFileInSimulationOutputDirectory('%s.txt' %fileName,'w')
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        # self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        # self.file.write('CellID CellType CellVolume\n')
        self.file.write('MCS CellType CellID COMx COMy COMz Cisplatin(uM)\n')
    def step(self,mcs):
        # print "cisplatin-concentration-to-file function is called every 4650 MCS (10 real-time minutes for Cisplatin in normal tissue)"
        print "cisplatin-concentration-to-file function is called every 12070 MCS (10 real-time minutes for Cisplatin in tumor tissue)"
        for cell in self.cellList:
            comPt=CompuCell.Point3D()
            field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
            # WORKS WHEN cell vol = 1 voxel; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual)
            # comPt.x=int(cell.xCM)
            # comPt.y=int(cell.yCM)
            # comPt.z=int(cell.zCM)
            comPt.x=int(round(cell.xCM/float(cell.volume))) # "divide by float = 0 error, 4-19-2013"
            comPt.y=int(round(cell.yCM/float(cell.volume)))
            comPt.z=int(round(cell.zCM/float(cell.volume)))
            cisplatin=field.get(comPt) # get concentration at center of mass
            # print 'cisplatin =',cisplatin
            # self.file.write('CELL ID=%d CELL TYPE=%d volume=%d\n' %(cell.id,cell.type,cell.volume))
            self.file.write('%d %d %d %d %d %d %f \n' %(mcs,cell.type,cell.id,comPt.x,comPt.y,comPt.z,cisplatin))
    def finish(self):
        # pass
        self.file.close() # close the file




class PlotTreatedCellsSteppable(SteppableBasePy):
# supported SVG color names at https://pythonhosted.org/ete2/reference/reference_svgcolors.html

    # def __init__(self,_simulator,_frequency=4650): # once per 10 minute approx. (1/465.189 min for Cisplatin diffusion 1 cell diam in normal tissue)
    def __init__(self,_simulator,_frequency=1207): # once per minute approx. (1/465.189 min for Cisplatin diffusion 1 cell diam in tumor tissue)
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
            # self.pW=self.addNewPlotWindow(_title='Number of cells exposed to cumulative therapeutic drug quantity',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Tumor cells')
            # self.pW.addPlot('TreatedCells',_style='Dots',_color='Red',_size=5)
            # self.pW=self.addNewPlotWindow(_title='Untreated Cells',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Tumor cells')
            # self.pW.addPlot('UntreatedCells',_style='Dots',_color='Blue',_size=5)
            # self.pW=self.addNewPlotWindow(_title='Intravenous (IV) Cisplatin',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 1/465.189 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            self.pW=self.addNewPlotWindow(_title='Intravenous (IV) Cisplatin',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 1/1207.183 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            self.pW.addPlot('Cisplatin',_style='Dots',_color='MediumSpringGreen',_size=4)
            # self.pW=self.addNewPlotWindow(_title='SKOV3 Cisplatin Accumulation',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            # self.pW.addPlot('SKOV3accumPerConc',_style='Dots',_color='MediumSpringGreen',_size=5)
            # self.pW.addPlot('TreatedCells',_style='Dots',_color='red',_size=5)
            # self.pW.addPlot('Cell1Vol',_style='Steps',_color='black',_size=5)
            # self.pW.setYAxisLogScale()
            # self.pW.set_xlim([5, 15])

    def step(self,mcs):
        # diffusion time for one cell diameter
        # tMins=mcs/465.189 # in normal tissue
        # tMins=mcs/1207.183 # in tumor tissue
        tHrs=mcs/1207.183/60 # in tumor tissue
        if 0<=mcs<cIPFirstPt2Hrs: # x-intercept (drug conc. = 0, becoming negative) > mcs > first IV data point for Cisplatin
            cisplatinIPtHrs=1856.1*tHrs
            self.pW.addDataPoint("Cisplatin",mcs,cisplatinIPtHrs) # name of the data series, x, y
            self.pW.showAllPlots()
            fileName="CisplatinIPplot"+str(mcs)+".txt"
            self.pW.savePlotAsPNG(fileName,1000,1000) # here we specify size of the image
        #if cIVZero>mcs>cFirst5Mins: # x-intercept (drug conc. = 0, becoming negative) > mcs > first IV data point for Cisplatin
        else:
            cisplatinIPtHrs=451.51*math.exp(-0.551*tHrs)
            # field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
            # cisplatinIV=max(field[])
            #Skov3accumperconc = -1.2e-07*cisplatin**3 + 1.3e-05*cisplatin**2 + 0.00058*cisplatin + 0.076
            # self.pW.addDataPoint("Treated Cells",mcs,treatedCells) # name of the data series, x, y
            # self.pW.addDataPoint("UntreatedCells",mcs,untreatedCells) # name of the data series, x, y
            self.pW.addDataPoint("Cisplatin",mcs,cisplatinIPtHrs) # name of the data series, x, y
            # self.pW.addDataPoint("Cisplatin",mcs,SKOV3accumPerConc) # name of the data series, x, y
            # self.pW.addDataPoint("Cell1Vol",mcs,cell1.volume) #name of the data series, x, y
            self.pW.showAllPlots()
            # fileName="NumCellsIC50Cisplatin_Plot"+str(mcs)+".png"
            fileName="CisplatinIPplot"+str(mcs)+".png"
            self.pW.savePlotAsPNG(fileName,1000,1000) # here we specify size of the image
            # fileName="CisplatinMax_txt"+str(mcs)+".txt"
            # self.pW.savePlotAsData(fileName)
            # fileName="CisplatinIVConc_txt"+str(mcs)+".txt"
            # self.pW.savePlotAsData(fileName)



class CellAccumToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of Cisplatin diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=121): # 0.1 minute of Cisplatin diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=113): # 0.1 hours of IgG diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=1): # for testing
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring + '_' + 'Accumulated' + drug + '_' + size + '_' + vess)
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        self.file.write('MCS cell.type CisplatinConcentrationMcrmolar CisplatinAccumulatedMcrmolar\n')
    def step(self,mcs):
        print "cell-list-to-file function is called every 113 MCS (1 real-time minutes for IgG in tumor)"
        for cell in self.cellList:
            if cell.type==2:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt) # get concentration at center of mass
                self.file.write('%d %d %f %f \n' %(mcs,cell.type,cisplatin,dictionaryAttrib[3])) # WRITE A LINE FOR EACH CELL 
            # if cell.type==3:
            #     dictionaryAttrib = CompuCell.getPyAttrib(cell)
            #     self.file.write('%d %d %d \n' %(mcs,cell.type,dictionaryAttrib[4]))
            if cell.type==8:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt) # get concentration at center of mass
                self.file.write('%d %d %f %f \n' %(mcs,cell.type,cisplatin,dictionaryAttrib[3])) # WRITE A LINE FOR EACH CELL 
    def finish(self):
        # pass
        self.file.close() # close the file




# PYTHONOUTPUTTOFILESTEPPABLE
class CellListToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=121): # 0.1 minute of diffusion, 1 cell diameter/MCS in tumor tissue
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
        # 2 lines following are mentioned in the manual
        # file=open('CellTotals.txt','a')
        ## FILE WILL BE IN ROOT OR C DRIVE:
        ## FILENAME "CellTotals_year-month-day_hour_minute_second.txt"

        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring+'_'+ drug+'_' +size+'_' +vess)
        # self.file=open('%s.txt' %fileName,'w')
        # self.openFileInSimulationOutputDirectory('%s.txt' %fileName,'w')
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        #pass
        self.file.write('MCS #UntreatedCancerCells #TreatedCancerCells\n')
    def step(self,mcs):
        print "cell-list-to-file function is called every 4650 MCS (10 real-time minutes for Cisplatin)"
        numPCancer=0
        numIC50Cancer=0
        for cell in self.cellList:
            if cell.type==2:
                numPCancer+=1
            if cell.type==8:
                numIC50Cancer+=1
            # self.file.write('CELL ID=%d CELL TYPE=%d volume=%d\n' %(cell.id,cell.type,cell.volume))
        self.file.write('%d %d %d \n' %(mcs,numPCancer,numIC50Cancer)) # WRITE AFTER ALL CELLS COUNTED
    def finish(self):
        # pass
        self.file.close() # close the file



