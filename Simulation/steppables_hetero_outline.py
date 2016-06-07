##### A model of bladder tumor metastasis heterogeneity under the most aggressive gemcitabine-cisplatin regimen.
#####
##### Based on NCCN May 2016 regimen guidelines for metastatic bladder cancer:
##### 21-day cycle for 6 cycles (metastasis), Gem 1250mg/m^2 over 30 min on days 1,8, followed by CISplatin 70 mg/m^2 IV over 60 min on day 1 or day 2
#####
##### components:
##### concentration IV (cisplatin or gemcitabine)
##### accumulation (cisplatin or gemcitabine, cell-line dependent)
##### removal from local concentration
#####
##### IC50 cell type switch -- 50% probability that I will not replicate (I will continue to collect drug however)
##### fork to IC50 non-replicating cells: die with certain probability
##### fork to dead cells: time frame to zero volume (instantaneous -- assume that the phagocytes (probably dedicated phagocytes such as macrophages, immature dendritic cells, and neutrophils () will be in the space, and then leave fluid)
##### 
##### outputs to file:
##### concentration at each cell
##### accumulation at each cell
#####     Although cellular accumulation of gemcitabine plateaus between 15-20 micromolar,
#####     the in vivo concentration during infusion of 1250 mg/m^2 does not go over
#####     ~ 15 micromolar, so we do not have to model a plateau.
##### cell counts of each cell type
##### 
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
#####
