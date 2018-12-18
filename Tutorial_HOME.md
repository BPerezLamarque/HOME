
**HOME Tutorial** -	**HO**st-**M**icrobiota **E**volution
====




This document indicates how to use our model of **HO**st-**M**icrobiota **E**volution. The first part details how to simulate mock microbiota. The second part details how to perform an empirical application; it is illustrated with the great apes microbiota data. The last part details how to interpret the HTML outputs of HOME.


**Citation:** Perez-Lamarque B., Morlon H. Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota (in prep.).


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


# Contents:
**Installation;**\
**Run Simulations;**\
**Run Empirical application;**\
**Interpretation of Results:**\
     *       Example 1: Results from a simulation with horizontal transmission;\
      *      Example 2: Results from a simulation with strict vertical transmission;\
       *     Example 3: Results from a simulation with environmental acquisition;


# Installation:
  
Our model is part on the R package RPANDA (Morlon et al., 2016) availbale on the CRAN or from gitHub.
  
  
```r
library(devtools)
install_github("hmorlon/PANDA",ref="Benoit")

```


# Simulations:


You can *provide a host tree* (e.g. an empirical tree) and simulate the evolution of a mock microbiota on it. Your tree must be binary, rooted and ultrametric. You can either directly provide a host tree in the function *sim_microbiota*, or you can have the host tree saved in your working directory (thus, its filename must be well-formated **host_tree_"name".tre** in a Newick format).


If you don't provide a host tree, the function *sim_microbiota* will randomly simulate a host tree (with a pure-birth process) according to the number of tips (i.e. number of host species) you have provided. 





## Parameters for simulations

```r

setwd("/my_working_directory/") # working directory where all the files and folders will be created. 


# name of your simulation
name <- "simulation_tree_1" 

# name of the different microbial OTU of your simulations (hereafter you will model the evolution of 6 OTUs on your host tree)
name_index <-  c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6")  

# simulated scenarii for each OTU: 
simul <- c(0,1,3,5,"indep","indep")
# each simulated OTU has its own evolutionary history:
# i) "0" simulates strict vertical transmission (i.e. O host-switch)
# ii) any positive integer simulates horizontal transmissions (with this given number of host-switches)
# iii) "indep" simulates environmental acquisition (i.e. independent evolutions) 



# Load the host tree to simulate the evolution of the OTUs on it 
host_tree <- read.tree("my_tree.tre")

# if you don't provide a host tree, you must add the number of tips of the simulated host tree
# n <- 20  



# simulated substitution rate 
simulated_mu <- 1  # simulated_mu=1 corresponds to on average 1 mutation per nucleotide

# total number of nucleotides in the alignment 
N <- 300  

# proportion of variable nucleotides in the alignment
proportion_variant <- 0.1 


# number of cores to run the analyses
nb_cores <- 1 # if you don't run in a multi-cores machine, the default value is 1

# seed used for simulations 
seed <- 1 



```

## Simulation of a microbiota

```r

sim_microbiota(name=name, name_index=name_index, simul=simul, 
provided_tree=host_tree, mu=simulated_mu, N=N, proportion_variant=proportion_variant, 
seed=seed, nb_cores=nb_cores)



```



Then, you can proceed to the **parameters estimation** by running HOME on the simulated microbiota.


```r



#  possible numbers of host-switches to test during the inference
lambda <- c(1:25) 

# number of trees (for Monte-Carlo estimation of the number of switches)
nb_tree <- 10000 

# number of randomizations in the model selection testing independent evolutions (R parameter)
nb_random <- 10

# rarefactions on the number of trees
raref <- FALSE # if TRUE rarefactions on the number of trees will be performed

```

```r

HOME_model(name=name, name_index=name_index, nb_tree=nb_tree, lambda=lambda, empirical=FALSE, raref=raref, 
nb_random=nb_random, seed=seed, nb_cores=nb_cores)

```


# Run Empirical applications: 

## Creating alignments for each OTU


In order to run HOME, you need first to create the microbial alignments for each OTU of the empirical microbiota. The first step consists in the usual processing the raw data into OTUs (here we propose to combine [QIIME](http://qiime.org) and [UPARSE](https://www.drive5.com/uparse/) thanks to the [BMP](http://www.brmicrobiome.org/16sillumina) pipelines). Thus, the second step makes the OTU alignments for running HOME on them (using our own bash script: for each core OTU, it picks the most abundant sequence for every host and align them). An example of scripts to use for preparing alignments before running HOME is available in a bash script [here](https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/make_clusters_OTU.sh).


All OTU alignments must be stored in a specific folder (i.e. the working directory where you will run HOME) with the filenames formatted as **alignment_"name"_"OTUXXX".fas**, where "name" will be the name our your HOME run ("name"), and "OTUXXX" is the name of the specific OTU ("name_index"). 



## Example of empirical applications - great apes microbiota: 

Let's run HOME for 3 OTUs from the great apes microbiota. First, your working directory must contain the host tree saved with a filename **host_tree_"name".tre** (Newick format), and all the OTU alignments with the filenames **alignment_"name"_"OTUXXX".fas** (FASTA format). You can also directly provide the host tree in the function *HOME_model*. 


For instance, in this empirical application, you must provide:

- */my_working_directory/host_tree_great_apes.tre*
- */my_working_directory/alignment_great_apes_OTU0001.fas*
- */my_working_directory/alignment_great_apes_OTU0002.fas*
- */my_working_directory/alignment_great_apes_OTU0003.fas*




You can directly download this example from RPANDA:  

```r

setwd("/my_working_directory/") # working directory where all the files and folders will be created. 


example_great_apes_microbiota(name="great_apes")

```


## Parameters of the empirical application 

```r

# name of the empirical application 
name <- "great_apes" 

# name of the different OTUs
name_OTU <-  c("OTU0001","OTU0002","OTU0003")


# possible numbers of host-switches to test during the inference
lambda <- c(1:25) 

# number of trees (for Monte-Carlo estimation of the number of switches)
nb_tree <- 10000 

# number of randomizations in the model selection testing independent evolutions 
nb_random <- 10

# rarefactions on the number of trees
raref <- FALSE # if TRUE rarefactions on the number of trees are performed
  
# number of cores to run the analyses
nb_cores <- 1 # if you don't run in a multi-cores machine, the default value is 1

# seed used for simulations 
seed <- 1 

```

## Run HOME on the empirical application 

```r


HOME_model(name=name, name_index=name_OTU, nb_tree=nb_tree, lambda=lambda, empirical=TRUE, raref=raref, 
nb_random=nb_random,  seed=seed, nb_cores=nb_cores)



```


# Interpret Results: 


The results of each HOME run are available in the folder "/my_working_directory/figures/". For each OTU, a HTMH file summaries all the results with tables and figures.


# Example 1: Results from a simulated OTU with horizontal transmission

## Description of the data - Simulation

<center>

Parameters  | Values
------------- | -------------
Number of symbiont-host association:|	20
Simulated substitution rate:|	1
Number of simulated switches:|	3
Seed for simulations:|	30
Sequences length:	|300
Probability variant sites:|	0.1
Number of invariant sites:|	281
Number of strictly variant sites:|	19

</center>

## Summary of the most likely scenario:

Parameters  | Values
------------- | -------------
Inferred substitution model:|	K80
Transition/Transversion ratio:|	0.6558
Estimated substitution rate:|	0.8519
Estimated number of switches:|	3
Associated likelihood:|	137.151
p-value: Strict vertical transmission:|	6e-04
p-value: Independent evolutions (nb switches):|	0
p-value: Independent evolutions (subs. rate):|	0


## Host-switches inference:
Most likely scenario estimated by the host-switches estimation: 

ksi | mu | -log(Likelihood)
------------- | ------------- | -------------
3 | 0.8519 | 137.151


Minus log likelihood as a function of the number of switches:

<p align="center">
    <img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_simul_C1_S65.png" width="500">
</p>


Estimated substitution rate as a function of the number of switches:

<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_mu_simul_C1_S65.png" width="500">
</p>




## Strict vertical transmission model:

Likelihood ratio test testing the model of strict vertical transmission (ksi=0). Strict vertical transmission is rejected if p-value < 0.05.


<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/pure_vertical_transmission_simul_C1_S65.png" width="500">
</p>


Results of the likelihood ratio test. The grey curve correponds to the Chi2 distribution with df=1. The dark blue line (resp. light) stands for the 0.05 (resp. 0.01) p-value threshold and the dashed orange line is the observed LRT ratio.


## Host-symbiont independent evolutions:

Model selection on independent evolutions.

Test  | p-values
------------- | -------------
Empirical ranking (ksi distribution)|	0
Empirical ranking (mu distribution)|	0



<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/independent_evolution_bootstrap_simul_C1_S65.png" width="500">
</p>

Representation of the estimated numbers of switches and the estimated substitution rates for the randomized alignments (in blue) and the empirical alignment (in orange). Independent evolutions can be rejected if the orange dot stands alone in the bottom left corner (i.e. rejected if p-values < 0.05).



## Estimated substitution model

*Estimated rate matrix:*

Rates  | A  |C |G |T
------------- | -------------| -------------| -------------| -------------
A |	-1  |	0.17  |	0.66 |	0.17 |	
C |	0.17  |	-1  |	0.17 |0.66	 |	
G |	0.66  |	 0.17 |	-1 |	0.17 |	
T |	0.17  |	0.66  |	0.17 |	-1 |	


*Nucleotide frequencies:*

A  |C |G |T
------------- | -------------| -------------| -------------
0.25  |	0.25  |	0.25 |	0.25 

## Simulated switches:

Simulated switch(es):  | 1 |2  |3
------------- | -------------| -------------| -------------
Branch origin  | 1 |36| 35
Branch arrival |35| 2 |28
Absolute position| 0.022 |0.089 |0.13



Host tree: 

<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/host_tree_simul_C1_S65.png" width="500">
</p>




# Example 2: Results from a simulated OTU with strict vertical transmission

## Description of the data - Simulation

Parameters  | Values
------------- | -------------
Number of symbiont-host association:|	20
Simulated substitution rate:|	1
Number of simulated switches:|	0
Seed for simulations:|	30
Sequences length:	|300
Probability variant sites:|	0.1
Number of invariant sites:|	279
Number of strictly variant sites:|	21

## Summary of the most likely scenario:

Parameters  | Values
------------- | -------------
Inferred substitution model:|	K80
Transition/Transversion ratio:|	0.70
Estimated substitution rate:|	0.70
Estimated number of switches:|	0
Associated likelihood:|	149.8
p-value: Strict vertical transmission:|	1
p-value: Independent evolutions (nb switches):|	0
p-value: Independent evolutions (subs. rate):|	0

*CONCLUSION: STRICT VERTICAL TRANSMISSION.*

## Host-switches inference:
Most likely scenario estimated by the host-switches estimation: 

ksi | mu | -log(Likelihood)
------------- | ------------- | -------------
3 | 0.70 | 149.8




Minus log likelihood as a function of the number of switches:
<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_simul_C1_S5.png" width="500">
</p>

Estimated substitution rate as a function of the number of switches:
<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_mu_simul_C1_S5.png" width="500">
</p>



## Strict vertical transmission model:

Likelihood ratio test testing the model of strict vertical transmission (ksi=0). Strict vertical transmission is rejected if p-value < 0.05.


<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/pure_vertical_transmission_simul_C1_S5.png" width="500">
</p>


Results of the likelihood ratio test. The grey curve correponds to the Chi2 distribution with df=1. The dark blue line (resp. light) stands for the 0.05 (resp. 0.01) p-value threshold and the dashed orange line is the observed LRT ratio.



## Host-symbiont independent evolutions:

Model selection on independent evolutions.

Test  | p-values
------------- | -------------
Empirical ranking (ksi distribution)|	0
Empirical ranking (mu distribution)|	0


<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/independent_evolution_bootstrap_simul_C1_S5.png" width="500">
</p>

Representation of the estimated numbers of switches and the estimated substitution rates for the randomized alignments (in blue) and the empirical alignment (in orange). Independent evolutions can be rejected if the orange dot stands alone in the bottom left corner (i.e. rejected if p-values < 0.05).




## Estimated substitution model

*Estimated rate matrix:*

Rates  | A  |C |G |T
------------- | -------------| -------------| -------------| -------------
A |	-1  |	0.15  |	0.7 |	0.15 |	
C |	0.15  |	-1  |	0.15 |0.7	 |	
G |	0.7  |	 0.15 |	-1 |	0.15 |	
T |	0.15  |	0.7  |	0.15 |	-1 |	


*Nucleotide frequencies:*

A  |C |G |T
------------- | -------------| -------------| -------------
0.25  |	0.25  |	0.25 |	0.25 




# Example 3: Results from a simulated OTU with environnemental acquisition (independent evolution)

## Description of the data - Simulation

Parameters  | Values
------------- | -------------
Number of symbiont-host association:|	20
Simulated substitution rate:|	1
Number of simulated switches:|	independent
Seed for simulations:|	30
Sequences length:	|300
Probability variant sites:|	0.1
Number of invariant sites:|	278
Number of strictly variant sites:|	22

## Summary of the most likely scenario:

Parameters  | Values
------------- | -------------
Inferred substitution model:|	K80
Transition/Transversion ratio:|	0.54
Estimated substitution rate:|	3.3793
Estimated number of switches:|	25
Associated likelihood:|	270.3
p-value: Strict vertical transmission:|	0
p-value: Independent evolutions (nb switches):|	0.56
p-value: Independent evolutions (subs. rate):|	0.89


*CONCLUSION:INDEPENDENT EVOLUTIONS (BASED ON THE NUMBER OF SWITCHES AND THE SUBSTITUTION RATE).*

## Host-switches inference:
Most likely scenario estimated by the host-switches estimation: 

ksi | mu | -log(Likelihood)
------------- | ------------- | -------------
25 | 3.38 | 270.3

Minus log likelihood as a function of the number of switches:

<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_simul_C1_S124.png" width="500">
</p>


Estimated substitution rate as a function of the number of switches:
<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/profil_model_switches_mu_simul_C1_S124.png" width="500">
</p>

## Strict vertical transmission model:

Likelihood ratio test testing the model of strict vertical transmission (ksi=0). Strict vertical transmission is rejected if p-value < 0.05.

<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/pure_vertical_transmission_simul_C1_S124.png" width="500">
</p>


Results of the likelihood ratio test. The grey curve correponds to the Chi2 distribution with df=1. The dark blue line (resp. light) stands for the 0.05 (resp. 0.01) p-value threshold and the dashed orange line is the observed LRT ratio.




## Host-symbiont independent evolutions:

Model selection on independent evolutions.

Test  | p-values
------------- | -------------
Empirical ranking (ksi distribution)|	0.56
Empirical ranking (mu distribution)|	0.88

<p align="center">
<img src="https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/independent_evolution_bootstrap_simul_C1_S124.png" width="500">
</p>

Representation of the estimated numbers of switches and the estimated substitution rates for the randomized alignments (in blue) and the empirical alignment (in orange). Independent evolutions can be rejected if the orange dot stands alone in the bottom left corner (i.e. rejected if p-values < 0.05).



## Estimated substitution model

*Estimated rate matrix:*

Rates  | A  |C |G |T
------------- | -------------| -------------| -------------| -------------
A |	-1  |	0.23  |	0.54 |	0.23 |	
C |	0.23  |	-1  |	0.23 |0.54	 |	
G |	0.54  |	 0.23 |	-1 |	0.23 |	
T |	0.23  |	0.54  |	0.23 |	-1 |	


*Nucleotide frequencies:*

A  |C |G |T
------------- | -------------| -------------| -------------
0.25  |	0.25  |	0.25 |	0.25 


