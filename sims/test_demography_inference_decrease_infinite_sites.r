source("../sim_wf.r")
source("../rare_variant_model.r")

#get args
args = commandArgs(trailingOnly=TRUE)

#get sample size
AN = as.numeric(args[1])
#get K
K = as.numeric(args[2])
#get target size
TS = as.numeric(args[3])
#get outfile
outfile = args[4]

num_replace = 5

#define model with decrease
Ne = c(14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,
      19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,
      29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,
      5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,
      14310,13292,13292,17803713,5362375)
T = c(70000,55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,
     21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,
     5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,
     1752,1579,1413,1252,1096,945,798,656,517,383,252,12,0)

T = (T[1]-T)#/(2*Ne[1])

# Ne = Ne/Ne[1]

R = rep(0,length(Ne))
R[length(R)-1] = -.1
R[length(R)-2] = .03

#simulate data
print("Simulate SFS1")
res1 = unlist(mclapply(1:10,function(i){sim_alleles(numeric(TS/10),function(s){pop_size(s,Ne,R,T)},0,.5,1e-9,0,max(T),AN)},mc.cores=10))
SFS1 = hist(res1,breaks=0:(AN+1),plot=FALSE,right=FALSE)

print("Simulate SFS2")
res2 = unlist(mclapply(1:10,function(i){sim_alleles(numeric(TS/10),function(s){pop_size(s,Ne,R,T)},0,.5,1e-8,0,max(T),AN)},mc.cores=10))
SFS2 = hist(res2,breaks=0:(AN+1),plot=FALSE,right=FALSE)

#convert to tibble for inference
SFS1_tibble = tibble(ref_context="AAA",alt_context="ATA",methylation_level=NA,AC=0:AN,n = SFS1$counts)
SFS2_tibble = tibble(ref_context="AAA",alt_context="ACA",methylation_level=NA,AC=0:AN,n = SFS2$counts)

#combine tibbles
SFS_tibble = bind_rows(SFS1_tibble,SFS2_tibble)

#define inference params
num_bin = 10
start_time = 0.0004
end_time = 0.025
burnin = 10

#make starting inference demography
T = rev(c(0,exp(seq(log(start_time),log(end_time),len=num_bin)),end_time+burnin))

T = T[1]-T #flip the T

Ne = rep(1,length(T))

R = rep(0,length(Ne))

#get estimates
#first, do a bunch of loops to get the best estimate
num_replace=num_bin
best_lnL = Inf
for (i in 1:40) {
	print(paste0("OPTIMIZATION TRY ", i))
	cur_opt = Ne_opt <- optimize_with_infinite_sites(SFS_tibble,num_replace,AN,0,Ne,R,T,nc=20,K=K,perturb_start=1)
	if (cur_opt$value < best_lnL) {
		best_lnL = cur_opt$value
		best_opt = cur_opt
	}
}

print(best_opt)

num_par = length(best_opt$par)

return_tibble = tibble(par = 1:num_par, est = exp(best_opt$par), lnL = best_opt$value, AN = AN)

print("Writing pop size")
write_tsv(return_tibble,outfile)

print("Get first order thetas")
par = exp(best_opt$par)
cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]
new_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
theta_opt = mutation_rate_FO(SFS_tibble,AN,new_Ne,R,T)

print("Writing theta")
write_tsv(mutate(theta_opt,AN=AN),paste0(outfile,".theta"))
