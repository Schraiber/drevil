# Inferring demography from large samples using DR EVIL
## Installation
All required functions for inference are in `rare_variant_model.r`. The libraries `parallel`, `cubature`, `nloptr`, and `tidyverse` are required, and can be installed using `install.pcakges`.

All required functions for simulation are in `sim_wf.r`. 

A number of additional useful functions are found in ``rare_alleles_figures_revision.ipynb``, mostly for plotting but some of them are useful for handling the output of inference. 

## Population size file format

Input to simulation and output of inference result in a population size history in three vectors, `a`, `r`, and `t`. 
The vector `t` denotes the start of each population size epoch, starting at `t[1] = 0` (the furthest in the past) and ending at `t[length(t)`], which indicates the present time.
Note that, depending on the use case, times may be in diffusion time units or in generations.
The vector `a`, the same length as `t` denotes the size of the populatoin at the start of each epoch. 
Note that in many cases, the last entry in `a` will be a duplicate of the second to last entry; this indicates the poulation size at the present is the same as the population size in the epoch prior to the present.
Note that, depending on the use case, sizes may be in diffusion units (i.e. relative to a fixed `N[0]`, or in raw population sizes).
The vector `r`, the same length as `t`, denotes the exponential growth rate during each epoch. During epoch i, the population grows as `a[i]*exp(r[i]*(t-t[i]))`. 
Note that, depending on the use case, rates may be in diffusion units, or raw growth rates.
Note that for now, built-in inference functions assume piecewise constant population histories (i.e. `r` is identically zero for all entries).
However, in principle site frequency spectra can be computed using non-zero growth rates in each epoch. 
Moreover, simulation with non-zero growth rates is also functional.

## Simulating data

Data can be simulated using the simple Wright-Fisher simulator in `sim_wf.r`. 
Any arbitrary population history can be simulated as long as an appropriate `N(t)` is specified, with `N(0)` being the most ancient point. 
Population sizes should be specified in raw population sizes and times in generations. 
To generate an automatic population size function given population size vectors `a`, `t`, and `r`, you can use the function `pop_size`,
i.e. `N = function(s) {pop_size(s,a,r,t)}`.

The function `sim_alleles(p0,N,s,h,mu1,mu2,tmax,ss=NULL)` can be used to simulate independent sites. 
`p0` is a vector of initial allele frequencies; the lenght of `p0` determines the number of sites to be simulated
`N` is the population size function described above
`s` is the selection coefficient and `h` is the dominance coefficient. Fitness is parameterized such that the fitness of the aa homozyogte is 1, Aa heterozygote is 1+hs, and the AA homozygote is 1+s. 
`mu1` is the mutation rate from a to A
`mu2` is the mutation rate from A to a.
`tmax` is the length of time to run the simulations.
`ss` is the sample size of a hypergeometric sample to be taken from the population. If `ss` is `NULL` (the default), no sample will be taken. 
