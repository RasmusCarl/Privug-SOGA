# Privug-SOGA
Privug backend featuring an exact Bayesian inference engine based on Gaussian mixtures.

## Contents

This repository contains the code to repoduce the experiments in the paper "Quantifying Privacy Risks with Gaussian Mixtures".

To reproduce all the plots in the paper run the notebook `SOSYM_Experiments`.

All programs used in the experiments can be found in the folder `programs/SOGA/Privug/rasmus`. 

<!---## Reproducibilty

- Create and start a new docker container based on the following steps:

docker container create -i -t --name SOGA bistrulli/soga:0.1
docker start SOGA

- then enter the container 

docker attach SOGA

- for reproduing Table 3 issue the following command

cd /root/SOGA/experimemts
python3 reproduce.py

- after executing this command the Table will be saved in /root/SOGA/experimemts/results/Table3.csv with the same structure of the Table 3 of the paper   

## Comparison with PSI

We removed the dependency with proprietary third-party tools to have a self-contained package. To this end, we replicated the PSI experiments via Mathics (https://mathics.org/) instead of Mathematica. However, if one would still like to use Mathematica for computing the PSI formula, the tool and a trial license can be requested on the Mathematica website (https://www.wolfram.com/mathematica/trial/). Once the license has been obtained, it is possible to run Mathematica to compute the symbolic formulas produced by PSI. For replicability SOGA, we save each PSI  formula in the folder "/root/SOGA/experiments/results/psi_formula" so that these can then be executed once a license for the tool has been obtained. -->


