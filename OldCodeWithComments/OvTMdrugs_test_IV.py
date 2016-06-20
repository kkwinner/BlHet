"""PARAMETER SOURCES
FINALIZE BEFORE FINAL RUNS

<!-- CISPLATIN -->
<!-- (Diffusion coefficient of Creatinine (used as Proxy for Cisplatin in Morrison, 1986) = 1.9 x 10^- 6 cm^2/s = 0.000114 cm^2/min = 11400 um^2/min = ) -->
<!-- D_{tissue} 
     = 1.778 x 10^{-4} (300.5)^{-0.75} = 0.000002466246792 cm^2/s = 246.624679177701 um^2/s (Swabb 1974 via Thurber 2011) 
     = 4652 (SKOV3.ip1 cell diameter)^2 / (10 min) -->
<!-- Cisplatin (Pt) accumulation at 5 $\mu$M in SKOV3  (Fig.3, Mistry, 1992)
     = 9.4E-15 $\mu$M/cell/min = 1.57E-16 muM/cell/sec 
     = 9.4E-14muM/cell/10 min -->
<!-- Cisplatin: accumulation at 100 micromolar for 50min in SKOV3 
     = 8 +/- 0.8 nmol/mg protein/50 mins: 8 +/- 0.8 nM * 1e-3uM/nM /mg protein * 0.21mg protein/uL * 1uL/10^9um^3 * 179.1 um^3/ SKOV3 cell /50min 
     = 6.02E-12 uM Pt/SKOV3 cell/min = 1.00E-13 uM Pt/SKOV3 cell/s 
     = 6.02E-11muM Pt/SKOV3 cell/10 min-->
<!-- Peak ave cisplatin intraperitoneal 
     (90mg/m^2 IP for 4 hour dwell = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; Howell, 1982): 
     40.5mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol = 134.977503749375 muM = 8098.65022496251 (concentration * 60 sec/MCS)-->
<!-- Cisplatin peak values IV
     (90mg/m^2IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; 1.6mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol =  5.3324445925679 uM,  
     = 319.946675554074 (concentration * 60 sec/MCS) -->
<!-- Minimum effective concentration of cisplatin 
     = 1 mug/mL; 1ug/ml*1000ml/L*1g/10^6ug* 1/300.5g/mol*10^6umol/mol = 3.33 muM-->
 <!-- secretion = micromolar = umol/
      (90mg/m^2IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL	1.6	1.6mug/mL*1g/10^6mug*1000mL/L/Xg/mol cisplatin*10^6umol/mol	5.3324445925679	uM	319.946675554074	(concentration * 60 sec/MCS)
      IC50 for SKOV3.ip1 cells: Xu W, et al.: ATP7B Antisense Oligodeoxynucleotides Increase the Cisplatin Sensitivity of Human Ovarian Cancer Cell Line SKOV3ip1. Int. J. of Gyn. Canc.-->
<!-- Go, Ata, RXMed:Platinol
     decay rate = ln(2)/half life = 0.23 / 10 min-->

<!-- TRASTUZUMAB -->
<!--http://www.druglib.com/druginfo/herceptin/description_pharmacology: In a study of women receiving adjuvant therapy for breast cancer, a mean half-life of trastuzumab of 16 days (range: 11-23 days) was observed after an initial dose of 8 mg/kg followed by a dose of 6 mg/kg every three weeks.  Between weeks 6 and 37, trastuzumabn serum concentrations reached a steady state with mean trough and peak concentrations of 63 mcg/mL and 216 mcg/mL, respectively.

mean of trough and peak plasma conc = (63+216)/2 = 139.5 mcg/ml
half-life = 1/(16days * 24h/day * 60min/h) = 0.0000434/min -->
<!-- A Systems Approach for Tumor Pharmacokinetics; Greg Michael Thurber, Ralph Weissleder; 2011
     antibody diffusion 
     = 10um^2/s 
     = 600 um^2/min 
     = 188.622 (SKOV3.ip1 cell diameter)^2 / (10 min) -->
<!-- Leveque, 2008
     half life = 0.000248 / 10 min
     decay rate = ln(2)/half life = 0.00017 / 10 min-->
<!-- absorption unknown - guessing 1% of min effective conc-->
<!-- peritoneal concentration unknown - guessing 10% less than plasma -->
"""





""" FOR PYTHON-ONLY SIMULATION
    INSERT AFTER COMPUCELL3DELMNT
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")

    # Basic properties of CPM (GGH) algorithm
    PottsElmnt.ElementCC3D("Dimensions",{"x":"100","y":"100","z":"1"})
    PottsElmnt.ElementCC3D("Steps",{},"1000")
    PottsElmnt.ElementCC3D("Temperature",{},"10.0")
    PottsElmnt.ElementCC3D("NeighborOrder",{},"1")
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})

    # Listing all cell types in the simulation
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"1","TypeName":"Tumor"})

    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
   # Module tracking center of mass of each cell

    #     PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"PixelTracker"})
    #     # Module tracking pixels of each cell

    #     PluginElmnt_3=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"BoundaryPixelTracker"})
    #     # Module tracking boundary pixels of each cell

    #     PluginElmnt_3.ElementCC3D("NeighborOrder",{},"2")

    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    # Specification of adhesion energies
    PluginElmnt_2.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"},"10.0")
    PluginElmnt_2.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Tumor"},"10.0")
    PluginElmnt_2.ElementCC3D("Energy",{"Type1":"Tumor","Type2":"Tumor"},"10.0")
    PluginElmnt_2.ElementCC3D("NeighborOrder",{},"1")

    PluginElmnt_3=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Chemotaxis"})

    # You may repeat ChemicalField element for each chemical field declared in the PDE solvers
    # Specification of chemotaxis properties of select cell types.

    PluginElmnt_4=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Secretion"})

    # Specification of secretion properties of select cell types.
    # You may repeat Field element for each chemical field declared in the PDE solvers
    # Specification of secretion properties of individual cells can be done in Python
    FieldElmnt=PluginElmnt_4.ElementCC3D("Field",{"Name":"Cisplatin"})
    FieldElmnt.ElementCC3D("ConstantConcentration",{"Type":"Tumor"},"10.0")
    FieldElmnt_1=PluginElmnt_4.ElementCC3D("Field",{"Name":"Trastuzumab"})
    FieldElmnt_1.ElementCC3D("ConstantConcentration",{"Type":"Tumor"},"10.0")

    SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"DiffusionSolverFE"})


#CISPLATIN
    # Specification of PDE solvers
    DiffusionFieldElmnt=SteppableElmnt.ElementCC3D("DiffusionField")
    DiffusionDataElmnt=DiffusionFieldElmnt.ElementCC3D("DiffusionData")
    DiffusionDataElmnt.ElementCC3D("FieldName",{},"Cisplatin")
    DiffusionDataElmnt.ElementCC3D("GlobalDiffusionConstant",{},"0.1")
    DiffusionDataElmnt.ElementCC3D("GlobalDecayConstant",{},"1e-05")
    # Additional options are:
    # DiffusionDataElmnt.ElementCC3D("InitialConcentrationExpression",{},"x*y")
    # DiffusionDataElmnt.ElementCC3D("ConcentrationFileName",{},"INITIAL CONCENTRATION FIELD - typically a file with path Simulation/NAME_OF_THE_FILE.txt")
    DiffusionDataElmnt.ElementCC3D("DiffusionCoefficient",{"CellType":"PCancerGFP"},"0.1")
    DiffusionDataElmnt.ElementCC3D("DecayCoefficient",{"CellType":"PCancerGFP"},"0.0001")
    SecretionDataElmnt=DiffusionFieldElmnt.ElementCC3D("SecretionData")
    # When secretion is defined inside DissufionSolverFEall secretio nconstants are scaled automaticly to account for extra calls of the solver when handling large diffusion constants
    
    # Uniform secretion Definition
    SecretionDataElmnt.ElementCC3D("Secretion",{"Type":"Tumor"},"0.1")
    # SecretionDataElmnt.ElementCC3D("SecretionOnContact",{"SecreteOnContactWith":"Tumor","Type":"Tumor"},"0.2")
    # SecretionDataElmnt.ElementCC3D("ConstantConcentration",{"Type":"Tumor"},"0.1")
    BoundaryConditionsElmnt=DiffusionFieldElmnt.ElementCC3D("BoundaryConditions")
    PlaneElmnt=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"X"})
    PlaneElmnt.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    PlaneElmnt_1=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"Y"})
    PlaneElmnt_1.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_1.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    PlaneElmnt_2=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"Z"})
    PlaneElmnt_2.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_2.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})    # Other options are (examples):
    # PlaneElmnt_1.ElementCC3D("Periodic")
    # PlaneElmnt_1.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":"10.0"})


#TRASTUZUMAB
    DiffusionFieldElmnt_1=SteppableElmnt.ElementCC3D("DiffusionField")
    DiffusionDataElmnt_1=DiffusionFieldElmnt_1.ElementCC3D("DiffusionData")
    DiffusionDataElmnt_1.ElementCC3D("FieldName",{},"Trastuzumab")
    DiffusionDataElmnt_1.ElementCC3D("GlobalDiffusionConstant",{},"0.1")
    DiffusionDataElmnt_1.ElementCC3D("GlobalDecayConstant",{},"1e-05")
    # Additional options are:
    # DiffusionDataElmnt_1.ElementCC3D("InitialConcentrationExpression",{},"x*y")
    # DiffusionDataElmnt_1.ElementCC3D("ConcentrationFileName",{},"INITIAL CONCENTRATION FIELD - typically a file with path Simulation/NAME_OF_THE_FILE.txt")
    DiffusionDataElmnt_1.ElementCC3D("DiffusionCoefficient",{"CellType":"Tumor"},"0.1")
    DiffusionDataElmnt_1.ElementCC3D("DecayCoefficient",{"CellType":"Tumor"},"0.0001")
    SecretionDataElmnt_1=DiffusionFieldElmnt_1.ElementCC3D("SecretionData")
    # When secretion is defined inside DissufionSolverFEall secretio nconstants are scaled automaticly to account for extra calls of the solver when handling large diffusion constants
    
    # Uniform secretion Definition
    # SecretionDataElmnt_1.ElementCC3D("Secretion",{"Type":"Tumor"},"0.1")
    # SecretionDataElmnt_1.ElementCC3D("SecretionOnContact",{"SecreteOnContactWith":"Tumor","Type":"Tumor"},"0.2")
    SecretionDataElmnt_1.ElementCC3D("ConstantConcentration",{"Type":"Tumor"},"0.1")
    BoundaryConditionsElmnt_1=DiffusionFieldElmnt_1.ElementCC3D("BoundaryConditions")
    PlaneElmnt_2=BoundaryConditionsElmnt_1.ElementCC3D("Plane",{"Axis":"X"})
    PlaneElmnt_2.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":"10.0"})
    PlaneElmnt_2.ElementCC3D("ConstantValue",{"PlanePosition":"Max","Value":"5.0"})
    # Other options are (examples):
    # PlaneElmnt_2.ElementCC3D("Periodic")
    # PlaneElmnt_2.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"10.0"})
    PlaneElmnt_3=BoundaryConditionsElmnt_1.ElementCC3D("Plane",{"Axis":"Y"})
    PlaneElmnt_3.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"10.0"})
    PlaneElmnt_3.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"5.0"})
    # Other options are (examples):
    # PlaneElmnt_3.ElementCC3D("Periodic")
    # PlaneElmnt_3.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":"10.0"})
    # CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    
    # CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    

    # SteppableElmnt_1=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"BlobInitializer"})
    SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"BlobInitializer"})
    
    # Initial layout of cells in the form of spherical (circular in 2D) blob
    # RegionElmnt=SteppableElmnt_1.ElementCC3D("Region")
    RegionElmnt=SteppableElmnt.ElementCC3D("Region")
    RegionElmnt.ElementCC3D("Center",{"x":"5","y":"5","z":"5"})
    RegionElmnt.ElementCC3D("Radius",{},"2")
    RegionElmnt.ElementCC3D("Gap",{},"0")
    RegionElmnt.ElementCC3D("Width",{},"5")
    RegionElmnt.ElementCC3D("Types",{},"VesselWall")
"""    

def configureSimulation(sim):
    import CompuCellSetup
    from XMLUtils import ElementCC3D    

    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"version":"3.7.2"})
    # CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    
    


# INCORPORATE XML DATA SPECIFICATION FILE
import sys
from os import environ

from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])



import CompuCellSetup

# path for Win
CompuCellSetup.setSimulationXMLFileName("Simulation/OvTMdrugs_test.xml")

# # path for Mac
# CompuCellSetup.setSimulationXMLFileName("/Users/squirrel/Desktop/Dropbox/research/CurrentSimulationsDropbox/OvTM/OvTMdrugs_test/Simulation/OvTMdrugs_test.xml")

#add additional

sim,simthread = CompuCellSetup.getCoreSimulationObjects()

#add additional attributes
# add extra attributes here
        
pyAttributeDictionaryAdder,dictAdder=CompuCellSetup.attachDictionaryToCells(sim)

CompuCellSetup.initializeSimulationObjects(sim,simthread)

#notice importing CompuCell to main script has to be done after
# call to getCoreSimulationObjects()
import CompuCell



#add lattice monitors here
#### changeWatcherREgistry is also initialized at start of Python steppables
#### CC3D crashes when changeWatcherRegistry in one of these lines is not commented out;
changeWatcherRegistry=CompuCellSetup.getChangeWatcherRegistry(sim)
stepperRegistry=CompuCellSetup.getStepperRegistry(sim)


#Add Python steppables here
from PySteppablesExamples import SteppableRegistry
#### changeWatcherREgistry is also initialized at start of lattice monitors
#### CC3D crashes when changeWatcherRegistry in one of these lines is not commented out;
# changeWatcherRegistry=CompuCellSetup.getChangeWatcherRegistry(sim)
steppableRegistry=SteppableRegistry()



from steppables_drugspheroid import CispIC50VisualizationSteppable
cispIC50VisualizationSteppable=CispIC50VisualizationSteppable (sim)
steppableRegistry.registerSteppable(cispIC50VisualizationSteppable)

# from steppables_drugspheroid import VolumeParamSteppable
# volumeParamSteppable=VolumeParamSteppable(sim)
# steppableRegistry.registerSteppable(volumeParamSteppable)

from steppables_drugspheroid import SetCellDictionaries
setCellDictionaries=SetCellDictionaries(sim)
steppableRegistry.registerSteppable(setCellDictionaries)

from steppables_drugspheroid import DiffusionSolverFESteeringCisplatinIV
diffusionSolverFESteeringCisplatinIV=DiffusionSolverFESteeringCisplatinIV(sim)
steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIV)

# from steppables_drugspheroid import DiffusionSolverFESteeringCisplatinIP
# diffusionSolverFESteeringCisplatinIP=DiffusionSolverFESteeringCisplatinIP(sim)
# steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIP)

from steppables_drugspheroid import SecretionSteppableCisplatin
secretionSteppableCisplatin=SecretionSteppableCisplatin(sim)
steppableRegistry.registerSteppable(secretionSteppableCisplatin)

# from steppables_drugspheroid import DiffusionSolverFESteeringCisplatinIPplusIV
# diffusionSolverFESteeringCisplatinIPplusIV=DiffusionSolverFESteeringCisplatinIPplusIV(sim)
# steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIPplusIV)

# from steppables_drugspheroid import SetCellConcentrations
# setCellConcentrations=SetCellConcentrations(sim)
# steppableRegistry.registerSteppable(setCellConcentrations)

# from steppables_drugspheroid import VolumeConstraintSteppable
# volumeConstraint=VolumeConstraintSteppable(sim)
# steppableRegistry.registerSteppable(volumeConstraint)

from steppables_drugspheroid import ChangeWithCisplatinSteppable
changeWithCisplatinSteppable=ChangeWithCisplatinSteppable(sim)
steppableRegistry.registerSteppable(changeWithCisplatinSteppable)

# from steppables_drugspheroid import MitosisSteppable
# mitosisSteppable=MitosisSteppable(sim)
# steppableRegistry.registerSteppable(mitosisSteppable)

# from steppables_drugspheroid import PrintAllCells
# printAllCells=PrintAllCells(sim)
# steppableRegistry.registerSteppable(printAllCells)

from steppables_drugspheroid import CellListToFileSteppable
cellListToFileSteppable=CellListToFileSteppable(sim)
steppableRegistry.registerSteppable(cellListToFileSteppable)

# from steppables_drugspheroid import CisplatinToFileSteppable
# cisplatinToFileSteppable=CisplatinToFileSteppable(sim)
# steppableRegistry.registerSteppable(cisplatinToFileSteppable)

from steppables_drugspheroid import CellAccumToFileSteppable
cellAccumToFileSteppable=CellAccumToFileSteppable(sim)
steppableRegistry.registerSteppable(cellAccumToFileSteppable)

# from steppables_drugspheroid import PlotTreatedCellsSteppable
# plotTreatedCellsSteppable=PlotTreatedCellsSteppable(sim)
# steppableRegistry.registerSteppable(plotTreatedCellsSteppable)

# from steppables_drugspheroid_chgWdrug_3_17_14 import InfoPrinterSteppable
# infoPrinterSteppable=InfoPrinterSteppable(sim)
# steppableRegistry.registerSteppable(infoPrinterSteppable)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)

