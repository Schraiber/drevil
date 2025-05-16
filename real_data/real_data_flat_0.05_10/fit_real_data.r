source("/scratch1/schraibe/rare_alleles/sims/rare_variant_model.r")
library(vroom)


#get args
args = commandArgs(trailingOnly=TRUE)

#get outfile
outfile = args[1]
num_bin = as.numeric(args[2])
start_time = as.numeric(args[3])
end_time = as.numeric(args[4])
burnin = as.numeric(args[5])

T = rev(c(0,exp(seq(log(start_time),log(end_time),len=num_bin)),end_time+burnin))

T = T[1]-T #flip the T

Ne = rep(1,length(T))

R = rep(0,length(Ne))

#read in data
SFS_to_analyze = vroom("/project/edgem_352/schraibe/gnomad_v4/SFS_to_analyze.tsv.gz")

#set up some parameters
K = 1000
max_AC = K
AN = 1000000
num_replace=num_bin
cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]
perturb_start = 1
start = 0

#filter data
SFS_to_analyze = SFS_to_analyze %>% 
    ungroup() %>% 
    filter(AC >= start,AC<=K)

#get first round
print("Running optimization with MOM theta")
Ne_opt = optimize_with_MOM_theta(SFS_to_analyze,num_replace,AN,0,Ne,R,T,nc=20,K=1000,perturb_start=1)

#set up second round
par = exp(Ne_opt$par)
new_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])

#get full optimization
print("Do the full optimization")
full_opt = optimize_full_likelihood_const(SFS_to_analyze, num_replace, AN, 0, new_Ne, R, T, nc = 20,rel_reduce = sqrt(.Machine$double.eps),max_iter=2,perturb_start = 0,K = K)

#write output
num_par = length(full_opt$Ne_opt$par)

return_tibble = tibble(par = 1:num_par, est = exp(full_opt$Ne_opt$par), lnL = full_opt$Ne_opt$value)

write_tsv(return_tibble,outfile)
