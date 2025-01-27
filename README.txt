This repository contains files related to the publication:

H. O. Caldag, S. Yesilyurt, Acoustic radiation forces on magnetically actuated helical swimmers, 
Physics of Fluids , 32, 092012, 2020. https://doi.org/10.1063/5.0020930

The paper introduces a methodology called chain-of-spheres. The method is an analytical approach to compute the acoustic radiation forces on intricate geometeries. Here we use the approach to compute the radiation forces on helical geometries and combine the model with a resistive force theory-based model for swimming to evaluate the trajectories of magnetized helices under acousto-magnetic actuation.

There are three folders:

- Force Calculation Only: This is the code used for the verification part in the paper. The code in this folder computes
the acoustic radiation force for a helix in the standing wave field. It does not simulate the trajectories.

- Standing Wave Simulation: The code in this folder simulates a magnetically actuated helical swimmer in a standing wave field over time in bulk swimming conditions. Main code of interest is traj_helix_standingwave.m

- Travelling Wave Simulation: The code in this folder simulates a magnetically actuated helical swimmer in a travelling wave field over time in bulk swimming conditions. Main code of interest is traj_helix_travellingwave.m

The codes themselves have many comments, check them out. Also check our paper covering this study (info below).

If you have any questions, write to hakanosmanc@gmail.com

Also if you use this methodology, don't forget to cite our study.
