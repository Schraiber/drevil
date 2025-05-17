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

## Inferrence 

### Data format

To use `DR EVIL` for inference, site frequency spectra must bein a specific file format. An example, derived from the gnomAD v4.1 dataset, can be found in `real_data/SFS_to_analyze.tsv.gz`. 
The first two columns,`ref_context` and `alt_context`, give the (possibly trinucleotide) context of the mutation. Note that these are essentially just labels, and can be any way of separating mutations into discrete classes.
`methylation` is the methylation level of the allele. In the published analysis, it is `NA` for sites that are not CpG transversions, and ranges in values from 1 to 16 for CpG transversions. Note that this column is again simply a label, and not used quantitatively at any step. 
Note that DR EVIL assume every unique combination of `ref_context`, `alt_context`, and `methylation_level` is a unique mutational context for which to assume a unique mutation rate.
`type` is used to indicate if a site is synonymous, missense, etc. Again, this is just a label and can be used to indicate whatever is desired. However, it is *NOT* used internally by DREVIL to stratify sites, unlike the previous labels.
`AC` is the allele count, ranging from 0 through the maximum allele frequency considered. NOTE: this column is *required* to have every value between 0 and the allele frequency cutoff for every context, and will cause mysterious errors if that is not true.
`n` is the number of sites observed for that context with `AC` alternative alleles. 

### Inference of mutation rates

If demographic history is fixed, mutation rates can be inferred using a suite of functions, `mutation_rate_MOM`, `mutation_rate_FO`, and `mutation_rate_MLE`.

`mutation_rate_MOM` uses the method of moments estimator of the mutation rate, and takes several arguments.
`counts` counts is a tibble of the format explained above
`K` is the allele frequency cutoff
`AN` is the sample size 
`Ne`, `R`, and `T` are the population history vectors described above (called `a`, `r`, and `t`, above)
`N0` can be given to convert results from diffusion units into real units
`nc` is the number of cores to use

`mutation_rate_FO` uses the first order estimator, i.e. estimates the mutation rate as proportional to the fraction of segregating sites.
All parameters are the same as for `mutation_rate_MOM`, except the parameters `K` and `nc` are not required.

`mutation_rate_MLE` computes the full maximum likelihood mutation rate estimate. 
All parameters are the same as `mutation_rate_MOM`, with the addition of the parameter `start` which says which allele frequency bin to use as the lowest frequency considered. `start = 0` includes monomorphic sites, while `start = 1` excludes monomorphic sites but uses all segregating sites.

### Inference of a distribution of mutation rates

Using the output of the previous functions, which estimate a single mutation rate for each context, it is possible to infer a distribution of mutation rates using the function `mutation_rate_dist_MLE`. 
First, you must join your input data with the mutation rate estimates.  If your input data is `SFS` and your mutation rates are `mr`, then you can simply use `SFS_with_rates = inner_join(SFS,mr)`.
Then, all parameters are the same as `mutation_rate_MLE`, except the first argument, `counts_with_theta`, is this new object with joined SFS and rates.
*NOTE* that if the population history parameters are in diffusion units, then mutation rates must also be in diffusion units. This is what the default will be if you use inferred demographic parameters to estimate mutation rates and set `N0 = NULL` when estimating rates.

### Inference of sequencing error

Similarly to inferring mutation rates, sequencing error rates can be inferred using the output of a function like `mutation_rate_MLE`. 
Note that it is not possible to simultaneously infer a distribution of mutation rates and sequencing error rates.
The function `sequencing_error_fit` takes all the same arguments as `mutation_rate_MLE`; unlike `mutation_rate_dist_MLE`, it takes the input data and mutation rate tables as separate input arguments, `SFS` and `mutation_rates`, respectively.

### Inference of population size history
We recommend a two step procedure for inferring population size history, first using `optimize_with_MOM_theta` to get near the MLE and following it up with `optimize_full_likelihood_const`.
In all cases, only piecewise constant population histories are considered.

For inference, an input, initial population history must be supplied. This essentially specifies the time bins for inference as well as any bits of population history that are not to be inferred. 
In the manuscript, we advise using a population history that is at equilibrium (aka "flat" past a certain time in the past).
Such a history can be initialized using the following code. Note that population sizes must be initialized in diffusion time units.
```
num_bin = 10
start_time = 0.0004
end_time = 0.025
burnin = 10

T = rev(c(0,exp(seq(log(start_time),log(end_time),len=num_bin)),end_time+burnin))

T = T[1]-T #flip the T

Ne = rep(1,length(T))

R = rep(0,length(Ne))
```
`num_bin` specifies the number of pieces of population history to infer. 
`start_time` defines the "present" in inference. Because bins are broken up on the log scale, this cannot be 0. We recommend a small number---here, we choose 0.0004 because it is corresponds to approximately 8 generations assuming an effective population size of 10,000. 
`end_time` corresponds to the most ancient time at which population sizes will be inferred. More anciently than this time, the population is assumed to be at equilibrium. 
`burnin` indicates the amount of time to "burn in" to equilibrum more anciently than `end_time`. We note that for rare alleles, a value as small as 1 or 2 is typically fine, but we recommend 10 as a standard setting. Typically, making this number larger does not increase the run time, as there is no simulation going on, it is simply accounted for with a numerical integral. 

With a base population size in hand, first use the function `optimize_with_MOM_theta`, which has several arguments.


## Running analyses from the manuscript
