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

### rowwise_event.csv: contains information on each contact (observation, trapping day with or without sample) of each individual
- id: identity of squirrel
- date: date of contact
- vole_pos: whether we have evidence for vole consumption (1) or not (0) from observation
- mammal_pos: whether there is metagenomic evidence for carnivory of any mammalian prey (1/0)
- carnivory_pos: metagenomic and observational evidence combined - any evidence of prey hunting and or consumption (1/0)
- sex: individual sex (F or M)
- individual age: A for adult, P for juvenile
- age_sex: age sex cateogry (combined)
- hum_reactivity: behavioral type reactivity to humans in the trap
- agon_PCA: PC1 score representing competitive ability (composite measure of win rate and initiation rate)
- n_vole_burrows: number of vole burrows detected within 15m of landmarks that the squirrel used (weighed by relative space use), divided by 10


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

### ILVs: Individual-level variables for each individual (static)
- id: identity of squirrel
- sex: F for female, M for male
- age: A for adult, P for juvenile
- age_sex: age sex category with abbreviations as above
- n_obs_positive: number of observations with vole hunting or consumption
- n_microtus_positive: number of fecal samples positive for vole DNA
- total_pos_voles: sum of n_obs_positive and n_microtus_positive
- hum_reactivity: proportion of trapping events during which individual showed fear responses (chatter, call, struggle)
- total_interactions: how many affiliative interactions individual was involved in
- initiation_af: initiation rate of affiliative interactions
- initiation_ag: initiation rate of agonistic interactions
- win_rat: win rate of agonistic interactions
- agon_PCA: score of PC1 resulting from initiation rate (agonistic) and win rate
- n_vole_burrows: number of vole burrows detected within 15m of landmarks that the squirrel used (weighed by relative space use), divided by 10
- n_trap_days: number of distinct days it was trapped
- n_obs_days: number of distinct observation days present in the study area
- n_feces: number of fecal samples collected
- n_mammals_DNA: number of fecal samples with mammalian DNA present
- n_vole_obs: same as n_obs_positive
- total_carnivory: total number of evidence for carnivory (from positive fecal samples and observations of hunting/consuming)

## NBDA vole eating: R code for replicating all analyses




