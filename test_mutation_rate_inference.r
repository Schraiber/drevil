source("sim_wf.r")
source("rare_variant_model.r")

#get args
args = commandArgs(trailingOnly=TRUE)

#get mutation rate
mu = as.numeric(args[1])
#get sample size
AN = as.numeric(args[2])
#get K
K = as.numeric(args[3])
#get target size
TS = as.numeric(args[4])
#get outfile
outfile = args[5]


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

#simulate data
sim_counts <- unlist(mclapply(1:10,function(i){sim_alleles(numeric(TS/10),function(s){pop_size(s,Ne,R,T)},0,.5,mu,0,max(T),AN)},mc.cores=10))

SFS = hist(sim_counts,breaks=0:AN,plot=FALSE,right=FALSE)$counts

#convert to tibble for inference
SFS_tibble = tibble(ref_context="AAA",alt_context="ATA",methylation_level=NA,AC=0:K,n=SFS[1:(K+1)])

#change demography to appropriate units
t = T/(2*Ne[1])

rho = Ne/Ne[1]

r = R = rep(0,length(Ne))

#get estimates

print("Get MOM")
MOM = mutation_rate_MOM(SFS_tibble,AN,rho,r,t,N0=Ne[1]) %>% mutate(method="MOM")

print("Get FO")
FO = mutation_rate_FO(SFS_tibble,AN,rho,r,t,N0=Ne[1]) %>% mutate(method="FO")

print("Get MLE")
MLE = mutation_rate_MLE(SFS_tibble,K,AN,rho,r,t,N0=Ne[1],start=0,nc=1) %>% mutate(method="MLE_target")

print("Get MLE_no_target")
MLE_no_target = mutation_rate_MLE(SFS_tibble,K,AN,rho,r,t,N0=Ne[1],start=1,nc=1) %>% mutate(method="MLE_no_target")

print("Write output")

res = bind_rows(MOM,FO,MLE,MLE_no_target) %>% mutate(mu=mu,AN=AN,TS=TS,K=K)

res %>% write_tsv(outfile)
 

