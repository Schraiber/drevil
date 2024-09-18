source("../sims/rare_variant_model.r")
library(vroom)


#get args
args = commandArgs(trailingOnly=TRUE)

#get outfile
outfile = args[1]

#starting model
Ne = c(14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,
      19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,
      29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,
      5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,
      14310,13292,14522,613285,5000000,5000000,5000000,5000000)
T = c(70000,55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,
     21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,
     5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,
     1752,1579,1413,1252,1096,945,798,656,517,383,252,124,50,25,12,0)

T = (T[1]-T)/(2*Ne[1])

Ne = Ne/Ne[1]

R = rep(0,length(Ne))

#read in data
SFS_to_analyze = vroom("SFS_to_analyze.tsv.gz")

#set up some parameters
K = 1000
max_AC = K
AN = 1000000
num_replace=5
cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]
perturb_start = 1
start = 0

#filter data
SFS_to_analyze = SFS_to_analyze %>% 
    ungroup() %>% 
    filter(AC >= start,AC<=K)

#get first round
Ne_opt = optimize_with_MOM_theta(SFS_to_analyze,num_replace,AN,0,Ne,R,T,nc=20,K=1000,perturb_start=1)

#set up second round
par = exp(Ne_opt$par)
new_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])

#get full optimization
full_opt = optimize_full_likelihood_const(SFS_to_analyze, num_replace, AN, 0, new_Ne, R, T, nc = 20,rel_reduce = sqrt(.Machine$double.eps),max_iter=2,perturb_start = 0,K = K)

#write output
num_par = length(full_opt$Ne_opt$par)

return_tibble = tibble(par = 1:num_par, est = exp(full_opt$Ne_opt$par), lnL = full_opt$Ne_opt$value)

write_tsv(return_tibble,outfile)
