source("sim_wf.r")
source("rare_variant_model.r")

#get args
args = commandArgs(trailingOnly=TRUE)

#get mean mutation rate
mean_mu = as.numeric(args[1])
#get sample size
AN = as.numeric(args[2])
#get K
K = as.numeric(args[3])
#get target size
TS = as.numeric(args[4])
#get CV
CV = as.numeric(args[5])
#get outfile
outfile = args[6]


#define Schiffels-Durbin-Agarwal model

Ne = c(14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,
      19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,
      29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,
      5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,
      14310,13292,14522,613285,5000000,5000000)
T = c(70000,55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,
     21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,
     5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,
     1752,1579,1413,1252,1096,945,798,656,517,383,252,124,50,0)


#Ne = c(10000, 100000, 1000000, 10000000, 10000000)

#T = c(4000, 3000, 2000, 1000, 0)

#Reverse time so it's in the right direction
T = (T[1]-T)

#no growth in any window
R = rep(0,length(Ne))

#get dist of mutation rates
sd_mu = mean_mu*CV
alpha = mean_mu^2/sd_mu^2
beta = mean_mu/sd_mu^2
mu = matrix(rgamma(TS,alpha,beta),nrow=10)

#simulate data
sim_counts <- unlist(mclapply(1:10, function(i){sim_alleles(numeric(TS/10),function(s){pop_size(s,Ne,R,T)},0,.5,mu[i,],0,max(T),AN)},mc.cores=10))

SFS = hist(sim_counts,breaks=0:AN,plot=FALSE,right=FALSE)$counts

#convert to tibble for inference
SFS_tibble = tibble(ref_context="AAA",alt_context="ATA",methylation_level=NA,AC=0:K,n=SFS[1:(K+1)])

#change demography to appropriate units
t = T/(2*Ne[1])

rho = Ne/Ne[1]

r = R = rep(0,length(Ne))

#get estimates

print("Get MLE")
MLE = mutation_rate_MLE(SFS_tibble,K,AN,rho,r,t,N0=Ne[1],start=0,nc=1)


print("Make data for dist")
#attach theta to original
SFS_for_dist = SFS_tibble %>% 
    inner_join(MLE)


print("Get dist")
MLE_dist = mutation_rate_dist_MLE(SFS_for_dist,K,AN,rho,r,t,nc=1)


print("Write output")

res = inner_join(MLE,MLE_dist,by=c("ref_context","alt_context","methylation_level"),suffix=c("_single","_dist"))

res %>% write_tsv(outfile)
 

