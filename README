=====================================
              DIFFUSE 
=====================================

1D explicit/implicit/time-independent 2nd-order diffusion solver

(c) Peter Woitke, 2017  (pw31@st-and.ac.uk)

------------------------------------------------------------------------
To checkout the git repository and compile the code, type

> git clone https://github.com/pw31/Diffusion
> cd Diffusion/src
> cp makefile.prodimo makefile
> make
> cd ..

The makefile.prodimo is for ifort compiler, adjust your own makefile 
if you want to compile e.g. with gfortran.

------------------------------------------------------------------------
To run a few test problems, type

> ./diffuse test1.in
> ./diffuse test2.in
> ./diffuse test3.in
> ./diffuse test4.in

------------------------------------------------------------------------
To plot the results after each run, type

> python Plot.py

test3 and test4 have analytic solutions, they will be automatically
overplotted by Plot.py

------------------------------------------------------------------------
To run a diffusion test problem in a planetary atmosphere, type

> ./diffuse atmos1.in
> ./diffuse atmos2.in
> ./diffuse atmos3.in

followed by

> python Plot.py
