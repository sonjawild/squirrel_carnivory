# Competitive ability predicts the spread of novel carnivory in California ground squirrels
Code and data for analysing prevalence of carnivory on mammalian prey in California ground squirrels and the spread of carnivory on voles through the population

## Folder: Raw data:
- Combined behavioral and DNA data on vole consumption_4-20-26: contains detailed records of each individual with and without evidence of carnivory on mammal prey from 100 fecal samples analysed. Additionally, we have added records of behavioral observations of vole hunting or consumption observed by members of our study population. Together, these were used to construct input data for analyses (see below)

## Folder: Input data
### edge_list.csv: contains social network data of focal squirrels observed on at least 3 observation days
- focal: identity of squirrel 1
- other: identity of squirrel 2 in the dyad
- assoc: simple ratio association index calculated from affiliate dyadic interactions
- trial: 1 for all (needed for NBDA)
- spatial: spatial overlap in home ranges (not used for this analysis)

### event.data.csv: contains information on when each individual was first seen hunting or consuming voles (observation and metagenomics)
- id: identity of squirrel
- time: time of first evidence of vole hunting or consumption (in days). 52 for right-censored individuals that never learned (end-time +1)
- trial: again 1 for all

### genome_summary.csv: table with evidence of mammalian prey in fecal samples (each row is a sample)
- age: A=adult, P=juvenile
- sex: F for female, M for male
- colony: Crow or Paradise (which colony an individual was trapped in)
- date: date fecal sample was collected
- mark: identity of squirrel
- Microtus californicus: >0 if evidence for vole DNA in fecal sample
- Thomomys bottae: >0 if evidence for presence of DNA
- Reithrodontomys megalotis: >0 if evidence for presence of DNA
- Mus musculus: >0 if evidence for presence of DNA
- Peromyscus truei: >0 if evidence for presence of DNA

### ILV_tv.csv: time-varying individual level variable (cumulative number of days present in study area + # fecal samples) - gets updated at every time step
- trial: again 1 for everyone
- id: identity of the squirrel
- time: time step in days
- presence: standardized presence within study (consisting of number of days observed on observation days + number of fecal samples collected up to each time step)

## NBDA vole eating: R code for replicating all analyses


