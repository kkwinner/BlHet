def configureSimulation(sim):
    import CompuCellSetup
    from XMLUtils import ElementCC3D

    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"version":"3.7.5"})
    # CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)
    

# INCORPORATE XML DATA SPECIFICATION FILE
import sys
from os import environ

from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

# path for Win (path for Mac may need full path)
CompuCellSetup.setSimulationXMLFileName("Simulation/BlHet.xml")

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


############################################### STEPPABLES


################################# DRUG DYNAMICS

# from steppablesBlHet import CispIC50VisualizationSteppable
# cispIC50VisualizationSteppable=CispIC50VisualizationSteppable (sim)
# steppableRegistry.registerSteppable(cispIC50VisualizationSteppable)

from steppablesBlHet import SetCellDictionaries
setCellDictionaries=SetCellDictionaries(sim)
steppableRegistry.registerSteppable(setCellDictionaries)

from steppablesBlHet import VolumeParamSteppable
volumeParamSteppable=VolumeParamSteppable(sim)
steppableRegistry.registerSteppable(volumeParamSteppable)

"""
from steppablesBlHet import DiffusionSolverFESteeringCisplatinIV
diffusionSolverFESteeringCisplatinIV=DiffusionSolverFESteeringCisplatinIV(sim)
steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIV)

# from steppablesBlHet import DiffusionSolverFESteeringCisplatinIP
# diffusionSolverFESteeringCisplatinIP=DiffusionSolverFESteeringCisplatinIP(sim)
# steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIP)

from steppablesBlHet import SecretionSteppableCisplatin
secretionSteppableCisplatin=SecretionSteppableCisplatin(sim)
steppableRegistry.registerSteppable(secretionSteppableCisplatin)

# from steppablesBlHet import DiffusionSolverFESteeringCisplatinIPplusIV
e# diffusionSolverFESteeringCisplatinIPplusIV=DiffusionSolverFESteeringCisplatinIPplusIV(sim)
# steppableRegistry.registerSteppable(diffusionSolverFESteeringCisplatinIPplusIV)

# from steppablesBlHet import SetCellConcentrations
# setCellConcentrations=SetCellConcentrations(sim)
# steppableRegistry.registerSteppable(setCellConcentrations)

from steppablesBlHet import ChangeWithCisplatinSteppable
changeWithCisplatinSteppable=ChangeWithCisplatinSteppable(sim)
steppableRegistry.registerSteppable(changeWithCisplatinSteppable)


################################# CELL DYNAMICS

# from steppablesBlHet import VolumeConstraintSteppable
# volumeConstraint=VolumeConstraintSteppable(sim)
# steppableRegistry.registerSteppable(volumeConstraint)

# from steppablesBlHet import MitosisSteppable
# mitosisSteppable=MitosisSteppable(sim)
# steppableRegistry.registerSteppable(mitosisSteppable)


################################# OUTPUTS

# from steppablesBlHet import PrintAllCells
# printAllCells=PrintAllCells(sim)
# steppableRegistry.registerSteppable(printAllCells)

from steppablesBlHet import CellListToFileSteppable
cellListToFileSteppable=CellListToFileSteppable(sim)
steppableRegistry.registerSteppable(cellListToFileSteppable)

# from steppablesBlHet import CisplatinToFileSteppable
# cisplatinToFileSteppable=CisplatinToFileSteppable(sim)
# steppableRegistry.registerSteppable(cisplatinToFileSteppable)

from steppablesBlHet import CellAccumToFileSteppable
cellAccumToFileSteppable=CellAccumToFileSteppable(sim)
steppableRegistry.registerSteppable(cellAccumToFileSteppable)

# from steppablesBlHet import PlotTreatedCellsSteppable
# plotTreatedCellsSteppable=PlotTreatedCellsSteppable(sim)
# steppableRegistry.registerSteppable(plotTreatedCellsSteppable)

# from steppablesBlHet_chgWdrug_3_17_14 import InfoPrinterSteppable
# infoPrinterSteppable=InfoPrinterSteppable(sim)
# steppableRegistry.registerSteppable(infoPrinterSteppable)
"""




CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)





"""PARAMETER SOURCES
FINALIZE BEFORE FINAL RUNS
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
