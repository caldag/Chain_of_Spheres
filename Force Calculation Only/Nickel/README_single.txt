- This folder contains the codes for the verification part in the paper (nickel).

- helix_comp_main.m computes the radiation force using the chain-of-spheres approach in a standing wave field.
  - The code is set up to try different minor radii and put forth as an example.
  - Some small modifications needed to test other geometric parameter sweeps.
- To generate Fig. 3 in the article (containing all comparisons), run arfcom_comp_all.m
- The .txt files contain the forces evaluated from the finite-element-based model.
  - The first column in the data files correspond to the geometric parameter swept (e.g. minor diameter)
  - The second column in the data files correspond to the acoustic radiation force.
- Check the codes, they have some explanatory comments! Also check our article if you are confused.
- And you can always write to hakan.caldag@sabanciuniv.edu. I will gladly try to help.

---------------------------------------------------------------------------

Explanations for other bits of code in the folder:

- sphbes1.m and sphbes2.m compute the spherical Bessel functions of the first and second kind.
- Fcalc.m computes the 

If you use this method in your studies, please cite our article:

H. O. Caldag, S. Yesilyurt, Acoustic radiation forces on magnetically actuated helical swimmers, 
Physics of Fluids , 32, 092012, 2020. https://doi.org/10.1063/5.0020930