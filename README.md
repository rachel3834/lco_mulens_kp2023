# lco_mulens_kp2023
## Simulations for the LCO Microlensing Key Project Proposal 2023

This repository contains code used to simulate the target selection and observing strategy for a Key Project proposal to the Las Cumbres Observatory Telescope Network. 

The kp_simulator_* Python notebooks contain simulations of the proposed observing strategy as implemented for targets discovered by several different surveys. 

Each simulator generates a large sample of alerts of microlensing events in progress, using realistic distributions of event parameters and the survey's spatial footprint. 
Simulated survey lightcurves for each event are then generated, using the realistic distribution of survey time sampling to different regions of the sky, and photometric properties. 
The simulator then generates the lightcurves that would be produced by the LCO network, implementing the proposed follow-up strategy, and taking into account 
the targets visibility from all sites in the LCO telescope network.  
Lastly, a microlensing model is fitted to both simulated lightcurves for each event, and the resulting parameters evaluated to determine how well the event parameters would be measured. 

The estimated_number_*_events* notebooks use the simulation results and the catalogs of (real) past events discovered by each survey to estimate the number of events that each survey 
would be likely to discover during the Key Project.  These notebooks also estimate the amount of time on LCO facilities required for the 
follow-up strategy, based on the simulated datasets. 

The kp_observing_time* notebooks combine the time required to follow-up targets from each of the different surveys, and calculate the overall time required per semester. 

Notebooks with numbers in suffix to the name indicate implementations of alternative strategies that were explored. 

The *.fits and *.hdf5 files contain the simulated datasets and event tables produced by the most recent iteration of follow-up strategy.
