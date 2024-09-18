library(parallel)
library(cubature)
library(nloptr)
library(tidyverse)
library(parallel)

normalize_vec = function(vec) { vec/sum(vec,na.rm=TRUE)}

pop_size = function(s,a,r,t) {
    num_epoch = length(a)
    if (length(r) != num_epoch | length(t) != num_epoch) {
        stop('a, r, and t should be the same length')
    }
    
    max_epoch = findInterval(s, t)
            
    pop_size = a[max_epoch]*exp(r[max_epoch]*(s-t[max_epoch])) #TODO: This is bad on the boundaries
    
    return(pop_size)
    
    
    
}

A = function(s,tmax,n,g,a,r,t,c=700) {
    num_epoch = length(a)
    if (length(r) != num_epoch | length(t) != num_epoch) {
        print(a)
        print(r)
        print(t)
        stop('a, r, and t should be the same length')
    }
    
    # I = sapply(s,function(x){max(which(t<=x))})
    I = findInterval(s, t)
            
    i = 1:(num_epoch-1)

    #only works for piecewise constant
    if (g == 0) {
        integral_bits = .5*(t[i+1]-t[i])/a[i]
        last_integral = .5*(s-t[I])/a[I]
    } else {
        integral_bits = (exp(-g/2*(t[i]-tmax))-exp(-g/2*(t[i+1]-t[i])-g/2*(t[i]-tmax)))/(a[i]*g)
        last_integral = (exp(-g/2*(t[I]-tmax))-exp(-g/2*(s-t[I])-g/2*(t[I]-tmax)))/(a[I]*g)
    }

    #works for piecewise exponential
    # integral_bits = ifelse(
    #     r[i] == 0 & g == 0,
    #     .5*(t[i+1]-t[i])/a[i],
    #     # exp(log(exp(-g/2*(t[i]-tmax))*(exp(-c)-exp(-g/2*(t[i+1]-t[i])-r[i]/2*(t[i+1]-t[i])-c))/(a[i]*(g+r[i])))+c)
    #     (exp(-g/2*(t[i]-tmax))-exp(-g/2*(t[i+1]-t[i])-r[i]/2*(t[i+1]-t[i])-g/2*(t[i]-tmax)))/(a[i]*(g+r[i]))
    #     # exp( -g/2*(t[i]-tmax) + log(1-exp(-g/2*(t[i+1]-t[i])-r[i]/2*(t[i+1]-t[i]))) )/(a[i]*(g+r[i]))
        
    # )        
    
    # last_integral = ifelse(
    #     r[I] == 0 & g == 0,
    #     .5*(s-t[I])/a[I],
    #     (exp(-g/2*(t[I]-tmax))-exp(-g/2*(s-t[I])-r[I]/2*(s-t[I])-g/2*(t[I]-tmax)))/(a[I]*(g+r[I]))
    #     # exp( -g/2*(t[I]-tmax) + log(1-exp(-g/2*(s-t[I])-r[I]/2*(s-t[I]))) )/(a[I]*(g+r[I]))
    # )

    cumulative_sum = cumsum(integral_bits)

    # Ensure that indexing is correct by appending a 0 at the start for the case when I = 1
    cumulative_sum = c(0, cumulative_sum)
    
    # Use the I vector to index into the cumulative_sum and add last_integral
    int = cumulative_sum[I] + last_integral
    
    # int = sapply(I, function(I) {sum(integral_bits[0:(I-1)])}) + last_integral
    
    return(int) #NB: I factor out the n
    
}

d_integrand = function(u,s,i,n,g,a,r,t) {
    Adiff = A(s,s,n,g,a,r,t)-A(u,s,n,g,a,r,t)
    res = g/2*(s-u) + (i-1)*log(Adiff) - (i+1)*log(1/n+(Adiff))
    # res = ifelse(is.na(res),0,res)
    res[is.na(res)]=0
    return(exp(res))
}

d_tilde = function(i,s,theta,n,g,a,r,t) {
    theta/(2*n)*hcubature(function(u) {d_integrand(u,s,i,n,g,a,r,t)}, lowerLimit = 0, upperLimit = s, 
        tol = .Machine$double.eps^0.25,vectorInterface=TRUE,maxEval=0)$integral
}

bell_polynomial_tilde_DP = function(n,d) {
    if (length(d)!=n){
        stop("Need same number of variables as the order you want to compute to")
    }
    
    b = numeric(length=n+1)
    
    b[1] = 1
        
    for (j in 2:(n+1)) {
        i = 1:(j-1)
        b[j] = sum(i/(j-1)*b[j-i]*d[i])
    }
    
    return(b)
    
}

bell_gamma_int_DP = function(n,c,alpha,beta,betap) {
    #the c are just the d but evaluated at theta = 1
    if (length(c)!=n) {
           stop("Need same number of variables as the order you want to compute!")
    }

    F = numeric(length=n)

    F[1] = exp(alpha*(log(beta)-log(betap)))

    for (j in 2:n) {
        i = 1:(j-1)
        F[j] = 1/betap*sum(((alpha-1)*i+j-1)/(j-1)*c[i]*F[j-i])
    }

    return(F)
}

exp_integrand = function(u,s,theta,n,g,a,r,t) {
    #NB: the factor of n is pulled out of A
    (exp(g/2*(s-u)))/(1/n+(A(s,s,n,g,a,r,t)-A(u,s,n,g,a,r,t)))
}

exp_integral = function(s,theta,n,g,a,r,t) {
    theta/2 * hcubature(function(u) {exp_integrand(u,s,theta,n,g,a,r,t)}, lowerLimit = 0, upperLimit = s, 
        tol = .Machine$double.eps^0.25,vectorInterface=TRUE)$integral
}

#compute p with fixed d and exp_integral
p_d = function(d,ei,theta,log=FALSE) {
    k = length(d)
    B_i = bell_polynomial_tilde_DP(k,theta*d)

    res = -theta*ei + log(B_i)
    if (log) {
        return(res)
    } else {
        return(exp(res))
    }
}

p_d_gamma = function(d,ei,alpha,beta,log=FALSE) {
    k = length(d)-1
    F_i = bell_gamma_int_DP(k+1,d,alpha,beta,beta+ei)

    res = log(F_i)
    if (log) {
        return(res)
    } else {
        return(exp(res))
    }
}

#compute p for a vector of thetas, they are all assumed to have the same d
#the idea here is different mutation categories
p_sep = function(k,s,theta,n,g,a,r,t,log=FALSE,nc=1) {
    num_tries = 10
    for (i in 1:num_tries) {
        tryCatch(
            {
                d_i = as.numeric(mclapply(1:k,d_tilde, s=s,theta=1,n=n,g=g,a=a,r=r,t=t,mc.cores=nc))
                break
            },
            error = function(cond) { message("failed"); NA }
        )
    }
    ei = exp_integral(s,1,n,g,a,r,t)

    matrix(sapply(theta,function(theta){p_d(d_i,ei,theta)}),nrow=k+1)
}

#Multicore version
p_mc = function(k,s,theta,n,g,a,r,t,nc=1,log=FALSE) {
    #compute integral, independent of k
    # exp_int = exp_integral(s,theta,n,g,a,r,t)
    exp_int = exp_integral(s,theta,n,g,a,r,t)
    #get all the ds
    d_i = as.numeric(mclapply(1:k,d_tilde, s=s,theta=theta,n=n,g=g,a=a,r=r,t=t,mc.cores=nc))
    #get all the bell polynomials
    B_i = bell_polynomial_tilde_DP(k,d_i)
    
    # return(list(exp=exp(-exp_int),bell=B_i,factorial=1/factorial(0:k)))
    
    #return the boy
    res = -exp_int + log(B_i)
    if (log) {
        return(res)
    } else {
        return(exp(res))
    }
}

p_inf_sites = function(k,s,theta,n,g,a,r,t,log=FALSE,nc=1) {
    #ONLY RETURNS THE VARIABLE SITES
    d_i = as.numeric(mclapply(1:k,d_tilde, s=s,theta=1,n=n,g=g,a=a,r=r,t=t,mc.cores=nc,mc.preschedule=TRUE))
    
    res = log(theta) + log(d_i)
        
    if (log) {
        return(res)
    } else {
        return(exp(res))
    }
    
}

p_gamma_mc = function(k,s,n,g,a,r,t,alpha,beta,log=FALSE,nc=1) {
    #compute integral, independent of k, setting theta = 1
    exp_int = exp_integral(s,theta=1,n,g,a,r,t)
    #get all the cs, which are ds with theta = 1
    #I DON'T FULLY UNDERSTAND WHY I NEED TO INDEX UP TO K+1
    c_i = as.numeric(mclapply(1:(k+1),d_tilde, s=s,theta=1,n=n,g=g,a=a,r=r,t=t,mc.cores=nc))
    #get the gamma integrals
    #I DON'T FULLY UNDERSTAND WHY I NEED TO INDEX UP TO K+1
    F_i = bell_gamma_int_DP(k+1,c_i,alpha,beta,beta+exp_int)
    
    #return
    # res = alpha*log(beta)-lgamma(alpha) + log(F_i)
    res = log(F_i)
    if (log) {
        return(res)
    } else {
        return(exp(res))
    }
    
    
}

mutation_rate_MOM = function(counts,AN,Ne,R,T,N0 = NULL) {
    ei = exp_integral(max(T),1,AN,0,Ne,R,T)

    res = counts %>%
        group_by(ref_context,alt_context,methylation_level) %>%
        summarise(p0 = n[AC==0]/sum(n)) %>%
        mutate(theta = -log(p0)/ei)
    

    if (is.null(N0)) {
        return(res)
    } else {
        return(res %>% mutate(mu = theta/(4*N0)))
    }
    
}

mutation_rate_FO = function(counts,AN,Ne,R,T,N0 = NULL) {
    ei = exp_integral(max(T),1,AN,0,Ne,R,T)

    res = counts %>%
        group_by(ref_context,alt_context,methylation_level) %>%
        summarise(p0 = n[AC==0]/sum(n)) %>%
        mutate(theta = (1-p0)/ei)
    

    if (is.null(N0)) {
        return(res)
    } else {
        return(res %>% mutate(mu = theta/(4*N0)))
    }
}

mutation_rate_MLE = function(counts,K,AN,Ne,R,T,N0 = NULL,start=0,nc=1) {

    d_i = as.numeric(mclapply(1:K,d_tilde, s=max(T),theta=1,n=AN,g=0,a=Ne,r=R,t=T,mc.cores=nc))
    ei = exp_integral(max(T),1,AN,0,Ne,R,T)
    
    
    res = counts %>%
        filter(AC>=start,AC<=K) %>% 
        group_by(ref_context,alt_context,methylation_level) %>%
        do({
            opt = optim(
                par = runif(1),
                fn = function(par){
                    theory = p_d(d_i,ei,par[1])
                    theory_norm = normalize_vec(theory[(start+1):length(theory)])
                    -sum(.$n*log(theory_norm))
                },
                method="Brent",
                lower=1e-5,
                upper=1
            )
            tibble(theta = opt$par, LL = opt$value, num_sites = sum(.$n))
        }) 

        if (is.null(N0)) {
            return(res)
        } else {
            return(res %>% mutate(mu = theta/(4*N0)))
        }
}

mutation_rate_dist_MLE = function(counts_with_theta,K,AN,Ne,R,T,N0 = NULL,start=0,nc=1) {

    message("Precompute")
    
    d_i = as.numeric(mclapply(1:(K+1),d_tilde, s=max(T),theta=1,n=AN,g=0,a=Ne,r=R,t=T,mc.cores=nc))
    ei = exp_integral(max(T),1,AN,0,Ne,R,T)
    
    message("Optimize")
    
    res = counts_with_theta %>%
        filter(AC>=start,AC<=K) %>% 
        group_by(ref_context,alt_context,methylation_level) %>%
        do({
            opt = optim(
                par = log(c(.$theta[1],.$theta[1])),
                fn = function(par){
                    par = exp(par)
                    alpha = par[1]^2/par[2]^2
                    beta = par[1]/par[2]^2
                    theory = p_d_gamma(d_i,ei,alpha,beta)
                    theory_norm = normalize_vec(theory[(start+1):length(theory)])
                    -sum(.$n*log(theory_norm))
                },
                method="L-BFGS-B",
                lower=log(c(1e-7,1e-7)),
                upper = log(rep(.1,2)),
                control = list(factr=10,ndeps=rep(1e-8,2))
            )
            tibble(mean_theta = exp(opt$par[1]), sd_theta = exp(opt$par[2]), LL = opt$value, n = sum(.$n))
        }) 

    message("Done!")

    if (is.null(N0)) {
        return(res)
    } else {
        return(res %>% mutate(mu = mean_theta/(4*N0)))
    }
}

optimize_with_infinite_sites = function(SFS, num_replace, n, g, Ne, R, T, nc = 1, ftol_rel = 1e-10, perturb_start = 0.1, K = 0.005 * n, start = 1) {
    if( start < 1) {
        print("Error: cannot optimize infinite sites with target size!")
        return(1)
    }
    
    counts = SFS %>% 
        filter(AC >= start,AC<=K) %>%  
        group_by(AC) %>%
        summarise(n = sum(n)) %>% #collapse all the mutational types together
        pull(n) 

    message(start)

    cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]

    # x0 = rnorm(num_replace,log(Ne[(length(Ne)-num_replace):(length(Ne)-1)]),perturb_start)
    # x0 = pmin(x0,rep(log(10000),num_replace))
    # x0 = pmax(x0,rep(log(1e-2),num_replace))

    if (perturb_start!=0) {
        best_LL = -Inf
	for (i in 1:20) {
            par = runif(num_replace,log(10),log(10000))
            par = exp(par)
            test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
            theory = p_inf_sites(K,max(T),1,n,g,test_Ne,R,T,nc=nc)
            theory_norm = normalize_vec(theory[start:length(theory)])
            LL = sum(counts*log(theory_norm),na.rm=TRUE)
            if (LL > best_LL) {
                best_LL = LL
                x0 = log(par)
            }
        }
    } else {
        x0 = log(Ne[(length(Ne)-num_replace):(length(Ne)-1)])
    }

    Ne_opt = optim(
        par = x0,
        fn = function(par) {
            par = exp(par)
            message(paste(par,collapse=" "))
            test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
            theory = p_inf_sites(K,max(T),1,n,g,test_Ne,R,T,nc=nc)
            theory_norm = normalize_vec(theory[start:length(theory)])
            LL = sum(counts*log(theory_norm),na.rm=TRUE)
            message(LL)
            -LL
        },
        method="L-BFGS-B",
        upper=rep(log(10000),num_replace),
        lower = rep(log(1e-2),num_replace),
        control=list(ndeps=rep(1e-8,num_replace))
    )

    return(Ne_opt)
}

optimize_with_MOM_theta = function(SFS, num_replace, n, g, Ne, R, T, nc = 1, ftol_rel = 1e-10, perturb_start = 0.1, K = 0.005 * n, start = 0) {
    counts = SFS %>% 
        filter(AC >= start,AC<=K) %>%  
        pull(n) %>% 
        matrix(nrow=K-start+1)

    cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]

    # x0 = rnorm(num_replace,log(Ne[(length(Ne)-num_replace):(length(Ne)-1)]),perturb_start)
    # x0 = pmin(x0,rep(log(10000),num_replace))
    # x0 = pmax(x0,rep(log(1e-2),num_replace))

    if (perturb_start!=0) {
        best_LL = -Inf
        #try to choose good starting parameters
	for (i in 1:20) {
            par = runif(num_replace,log(10),log(10000))
            par = exp(par)
            test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
            ei = exp_integral(max(T),1,n,0,test_Ne,R,T)
            p0 = counts[1,]/colSums(counts)
            theta = -log(p0)/ei
            theory = p_sep(K,max(T),theta,n,g,test_Ne,R,T,nc=nc)
            theory_norm = apply(as.matrix(theory[(start+1):nrow(theory),]),2,normalize_vec)
            LL = sum(counts*log(theory_norm),na.rm=TRUE)
            if (LL > best_LL) {
                best_LL = LL
                x0 = log(par)
            }
        }
    } else {
        x0 = log(Ne[(length(Ne)-num_replace):(length(Ne)-1)])
    }

    Ne_opt = optim(
        par = x0,
        fn = function(par) {
            par = exp(par)
            message(paste(par,collapse=" "))
            test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
            ei = exp_integral(max(T),1,n,0,test_Ne,R,T)
            p0 = counts[1,]/colSums(counts)
            theta = -log(p0)/ei
            theory = p_sep(K,max(T),theta,n,g,test_Ne,R,T,nc=nc)
            theory_norm = apply(as.matrix(theory[(start+1):nrow(theory),]),2,normalize_vec)
            LL = sum(counts*log(theory_norm),na.rm=TRUE)
            message(LL)
            -LL
        },
        method="L-BFGS-B",
        upper=rep(log(10000),num_replace),
        lower = rep(log(1e-2),num_replace),
        control=list(ndeps=rep(1e-8,num_replace))
    )

    return(Ne_opt)
}

optimize_full_likelihood_const = function(SFS, num_replace, n, g, Ne, R, T, nc = 1,rel_reduce = sqrt(.Machine$double.eps),max_iter=20,perturb_start = .1,K = 0.005*n,start=0) {
    
    #Get counts matrix
    counts = SFS %>% 
        arrange(ref_context,alt_context) %>% 
        filter(AC >= start,AC<=K) %>%  
        pull(n) %>% 
        matrix(nrow=K-start+1)

    message("Initialize")

    theta_opt = mutation_rate_MLE(SFS,K=K,AN=n,Ne=Ne,R=R,T=T,start=start,nc=nc)

    message(sum(theta_opt$LL))
    
    cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]    
    
    Ne_opt = bobyqa(
        x0 = pmin(rnorm(num_replace,log(Ne[(length(Ne)-num_replace):(length(Ne)-1)]),perturb_start),rep(log(10000),num_replace)),
        fn = function(par) {
            par = exp(par)
            message(paste(par,collapse=" "))
            test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
            theory = p_sep(K,max(T),theta_opt$theta,n,0,test_Ne,R,T,nc=nc)
            theory_norm = apply(as.matrix(theory[(start+1):nrow(theory),]),2,normalize_vec)
            LL = sum(counts*log(theory_norm))
            message(LL)
            -LL
        },
        upper=rep(log(10000),num_replace),
        control=list(ftol_rel = 1e-10)#sqrt(.Machine$double.eps))
    )

    old_LL = Ne_opt$value

    for (i in 2:max_iter) {
        message(paste0("Iteration ", i))
        Ne = c(Ne[1:(length(Ne)-(num_replace+1))],exp(Ne_opt$par[1:num_replace]),exp(Ne_opt$par[num_replace]))

        theta_opt = mutation_rate_MLE(SFS,K=K,AN=n,Ne=Ne,R=R,T=T,start=start,nc=nc)

        message(sum(theta_opt$LL))
    
        cur_Ne = Ne[1:(length(Ne)-(num_replace+1))]    

        old_Ne_opt = Ne_opt
        
        Ne_opt = bobyqa(
            x0 = pmin(rnorm(num_replace,log(Ne[(length(Ne)-num_replace):(length(Ne)-1)]),perturb_start/sqrt(i)),rep(log(10000),num_replace)),
            fn = function(par) {
                par = exp(par)
                message(paste(par,collapse=" "))
                test_Ne = c(cur_Ne,par[1:num_replace],par[num_replace])
                theory = p_sep(K,max(T),theta_opt$theta,n,0,test_Ne,R,T,nc=nc)
                theory_norm = apply(as.matrix(theory[(start+1):nrow(theory),]),2,normalize_vec)
                LL = sum(counts*log(theory_norm))
                message(LL)
                -LL
            },
            upper=rep(log(10000),num_replace),
            control=list(ftol_rel = 1e-10)#sqrt(.Machine$double.eps))
        )

        new_LL = Ne_opt$value

        rel_change = (new_LL-old_LL)/old_LL

        #If we did worse, just try again!
        if (rel_change > 0) {
            Ne_opt = old_Ne_opt
            next
        }

        message(rel_change)

        if (abs(rel_change) < rel_reduce) { break } 

        old_LL = new_LL
    
    }

    return(list(theta_opt=theta_opt,Ne_opt=Ne_opt))
}
