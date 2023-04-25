## GS_modules
Modules of GWS using JWAS for prediction. Warning :warning: : These modules have been tailored to work with the Cucumber breeding program data some modifications might be needed so it can work for other crops.

**Domino environment** : PMC and Duplicate of R 3.5 and Python 3.7 -- EVA and Asreml - Revision #6

To run the master table using the python and origin selection code (diallel_mod) the domino environment is Duplicate of R 3.5 and Python 3.7 -- EVA and Asreml
check on the $Environment variables VAULT_APP_ROLE_ID and SECRET to make sure you have those credentials

**Phenotype module** :pencil:

![alt text](https://github.platforms.engineering/ELBFA/GS_modules/blob/main/pheno_module.png?raw=true)

**Genotype module** :dna:
step 1 use inbred geno format convert HD format of the DH lines into 0,1,2
step 2 use the origin_git code, this code will:

- create a table of the hybrids with GermID parent 1 and GermID parent 2 to be used as a reference to create the synth hybrids
- generate synth hybrids in HD/Beagle format
- generate synth hybrids converted into 0,1,2 format

This code uses sql to retrieve pedigree name and GermID parameters for several generations of input data
Also the code generates a file which contains the name of missing parents that could not be found in the genotype data, this can be shares with the breeding team/genomics team to follow up with remediation, usually the breeding team collects more than one sammple of the parents so they could be re-genotyped 

![alt text](https://github.platforms.engineering/ELBFA/GS_modules/blob/main/geno_module.png?raw=true)

**Statistical module**

Implementation of statistical models in JWAS
- GBLUP

- ABLUP 

- SSGBLUP

**Formatf90 module**
Code to format input genos and phenos to be used as input data in the Single step model using blupF90 (JingYang)

**Diallel module**
Example code that were run to provide results of diallel table for different scenarios
- PCM3 x PCM3
- PCM1 x new testers same season
- PCM1 x testers of a different season 

Master table and Origin prediction : example files of how the master table code from Xing and the Origin prediction code from Yujing can be run within R to obtain a diallel table with additional information on marker data and if available het groups

**Accuracy module**
Using experimental stage subsets for validation of models

- GBLUP

- ABLUP 

- SSGBLUP

