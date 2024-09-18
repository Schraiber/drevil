pop_size = function(s,a,r,t) {
    num_epoch = length(a)
    if (length(r) != num_epoch | length(t) != num_epoch) {
        stop('a, r, and t should be the same length')
    }
    
    max_epoch = sapply(s,function(x){max(which(t<=x))})
            
    pop_size = a[max_epoch]*exp(r[max_epoch]*(s-t[max_epoch])) #TODO: This is bad on the boundaries
    
    return(pop_size)
    
    
    
}

after_mutation = function(mu1,mu2,p) {
    p*(1-mu2)+(1-p)*mu1
}

after_selection = function(s,h,p) {
    p*(1 + s*h*(1 - p) + s*p)/(1 + s*p*(2*h*(1 - p) + p))
}

sim_alleles = function(p0,N,s,h,mu1,mu2,tmax,ss=NULL) {
    p = p0
    n = length(p0)
    message("generation 1")
    update_freq = floor(tmax/10)
    cur_update = update_freq
    for (t in 1:tmax) {
        if (t == cur_update) {
            message(paste0("generation ",t))
            cur_update = cur_update + update_freq
        }
        curN = floor(N(t))
        # message(paste(p,collapse=" "))
        p = after_mutation(mu1,mu2,p)
        # message(paste(p,collapse=" "))
        p = after_selection(s,h,p)
        # message(paste(p,collapse=" "))
        x = rbinom(n,2*curN,p)
        if (any(is.na(x))) {
            message(paste(p,collapse=" "))
            message(paste(x, collapse=" "))
            message(curN)
            message(t)
            invisible(readline(prompt="Press [enter] to continue"))
        }
        p = x/(2*curN)
        # message(paste(p,collapse=" "))
        # invisible(readline(prompt="Press [enter] to continue"))
    }
    if (is.null(ss)) {
        x
    } else {
        print(ss)
        print(curN)
        rhyper(n,x,2*curN-x,ss)
    }
}
