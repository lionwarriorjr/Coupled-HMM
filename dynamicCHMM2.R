require(MASS)
require(ggplot2)
############################################################################
#                       Forward/Backward Recursions                        #
############################################################################
compute_obs_prob <- function(y_t,states){
  num_states <- length(states)
  mu <- 0
  for(m in 1:M){
    start_index <- starts[m]
    end_index <- ends[m]
    W_m <- W[,start_index:end_index]
    one_hot <- rep(0,dim[m])
    one_hot[states[m]] <- 1
    one_hot <- matrix(one_hot)
    mu <- mu+W_m%*%one_hot
  }
  result <- (det(C)^-0.5)%*%((2*pi)^(-D/2))%*%exp(-0.5%*%t(y_t-mu)%*%solve(C)%*%(y_t-mu))
  result <- result[1]
  result
}
forward <- function(){
  alphas <- array(0,c(dim,T_))
  for(t in 1:T_){
    alphas <- forward_unwrap(alphas,t,NULL)
  }
  alphas
}
forward_unwrap <- function(alphas,t,chains){
  if(length(chains) == M){
    obs.prob <- compute_obs_prob(Y[t],chains)
    if(t == 1){
      prior.prod <- 1.0
      sapply(1:M,function(m){
        prior.prod <<- prior.prod*PI[[m]][chains[m]]
      })
      indexes <- matrix(c(chains,t),1)
      alphas[indexes] <- obs.prob*prior.prod
    }else{
      recur <- forward_recursion(alphas,chains,NULL,t)
      indexes <- matrix(c(chains,t),1)
      alphas[indexes] <- obs.prob*recur
    }
  }else{
    index <- length(chains)+1
    for(k in 1:dim[index]){
      states <- c(chains,k)
      alphas <- forward_unwrap(alphas,t,states)
    }
  }
  alphas
}
forward_recursion <- function(alphas,curr,prev,t){
  if(length(prev) == M){
    indexes <- matrix(c(prev,t-1),1)
    result <- alphas[indexes]
  }else{
    index <- length(prev)+1
    next_state <- curr[index]
    transmat <- P[[index]]
    K <- dim[index]
    result <- 0
    for(kminus in 1:K){
      trans.prob <- transmat[kminus,next_state]
      prv <- c(prev,kminus)
      result <- result+trans.prob*forward_recursion(alphas,curr,prv,t)
    }
  }
  result
}
backward <- function(){
  betas <- array(0,c(dim,T_))
  for(t in T_:1){
    betas <- backward_unwrap(betas,t,NULL)
  }
  betas
}
backward_unwrap <- function(betas,t,chains){
  if(length(chains) == M){
    indexes <- matrix(c(chains,t),1)
    betas[indexes] <- backward_recursion(betas,NULL,chains,t)
  }else{
    index <- length(chains)+1
    for(kminus in 1:dim[index]){
      states <- c(chains,kminus)
      betas <- backward_unwrap(betas,t,states)
    }
  }
  betas
}
backward_recursion <- function(betas,curr,prev,t){
  if(length(curr) == M){
    result <- compute_obs_prob(Y[t],curr)
    if(t < T_){
      indexes <- matrix(c(curr,t+1),1)
      result <- result*betas[indexes]
    }
  }else{
    index <- length(curr)+1
    prev_state <- prev[index]
    transmat <- P[[index]]
    K <- dim[index]
    result <- 0
    for(k in 1:K){
      trans.prob <- transmat[prev_state,k]
      nxt <- c(curr,k)
      result <- result+trans.prob*backward_recursion(betas,nxt,prev,t)
    }
  }
  result
}
normalize_gamma <- function(x,n,part.func){
  d <- dim(x)
  d.new <- d[-length(d)]
  block.size <- prod(d.new)
  x[(block.size*(n-1)+1):(block.size*n)] <-
    x[(block.size*(n-1)+1):(block.size*n)]/part.func
  array(x,dim=d)
}
compute_gamma <- function(alphas,betas){
  gammas <<- array(0,c(dim,T_))
  for(t in 1:T_){
    result <- gamma_unwrap(gammas,alphas,betas,t,NULL)
    gammas <<- normalize_gamma(gammas,t,result)
  }
  gammas
}
gamma_unwrap <- function(gammas,alphas,betas,t,chains){
  if(length(chains) == M){
    indexes <- matrix(c(chains,t),1)
    gammas[indexes] <<- alphas[indexes]*betas[indexes]
    result <- alphas[indexes]*betas[indexes]
  }else{
    index <- length(chains)+1
    result <- 0
    for(k in 1:dim[index]){
      states <- c(chains,k)
      result <- result+gamma_unwrap(gammas,alphas,betas,t,states)
    }
  }
  result
}
############################################################################
#                              E step updates                              #
############################################################################
compute_E_St <- function(t,m,chains){
  if(length(chains) == M){
    one_hot <- rep(0,dim[m])
    one_hot[chains[m]] <- 1
    indexes <- matrix(c(chains,t),1)
    result <- one_hot*gamma.mat[indexes]
  }else{
    index <- length(chains)+1
    result <- rep(0,dim[m])
    for(k in 1:dim[index]){
      states <- c(chains,k)
      result <- result+compute_E_St(t,m,states)
    }
  }
  result
}
get_e_St <- function(){
  e_St <- list()
  for(t in 1:T_){
    e_St[[t]] <- list()
    for(m in 1:M){
      e_St[[t]][[m]] <- compute_E_St(t,m,NULL)
    }
  }
  e_St
}
compute_E_St_E_St <- function(t,m,n,chains){
  if(length(chains) == M){
    o1 <- rep(0,dim[m])
    o2 <- rep(0,dim[n])
    o1[chains[m]] <- 1
    o2[chains[n]] <- 1
    indexes <- matrix(c(chains,t),1)
    result <- matrix(o1)%*%t(matrix(o2))*gamma.mat[indexes]
  }else{
    index <- length(chains)+1
    result <- matrix(0,dim[m],dim[n])
    for(k in 1:dim[index]){
      states <- c(chains,k)
      result <- result+compute_E_St_E_St(t,m,n,states)
    }
  }
  result
}
get_e_St_m_St_n <- function(){
  e_St_m_St_n <- list()
  for(t in 1:T_){
    e_St_m_St_n[[t]] <- list()
    for(m in 1:M){
      e_St_m_St_n[[t]][[m]] <- list()
      for(n in 1:M){
        e_St_m_St_n[[t]][[m]][[n]] <- compute_E_St_E_St(t,m,n,NULL)
      }
    }
  }
  e_St_m_St_n
}
compute_E_t_minus_t <- function(t,m,states.prev,states.curr){
  if(length(states.prev) == M){
    m1 <- states.prev[m]
    m2 <- states.curr[m]
    trans.prob <- 1.0
    sapply(1:M,function(m){
      trans.prob <<- trans.prob*P[[m]][states.prev[m],states.curr[m]]
    })
    indexes <- matrix(c(states.prev,t-1),1)
    alpha <- alpha.mat[indexes]
    indexes <- matrix(c(states.curr,t),1)
    beta <- beta.mat[indexes]
    obs.prob <- compute_obs_prob(Y[t],states.curr)
    p <- alpha*trans.prob*obs.prob*beta
    result <- p
    index[[t]][[m]][m1,m2] <<- p
  }else{
    chain_id <- length(states.prev)+1
    result <- 0
    for(k in 1:dim[chain_id]){
      curr <- c(states.curr,k)
      for(kminus in 1:dim[chain_id]){
        prev <- c(states.prev,kminus)
        result <- result+compute_E_t_minus_t(t,m,prev,curr)
      }
    }
  }
  result
}
get_e_St_minus_m_St_m <- function(){
  e_St_minus_m_St_m <- matrix(0,T_,M)
  for(t in 2:T_){
    normalizer <- 0
    index[[t]] <<- list()
    for(m in 1:M){
      index[[t]][[m]] <<- matrix(0,dim[m],dim[m])
      expect <- compute_E_t_minus_t(t,m,NULL,NULL)
      normalizer <- normalizer+expect
      e_St_minus_m_St_m[t,m] <- expect
    }
    normalizer <- normalizer[1]
    e_St_minus_m_St_m[t,] <- e_St_minus_m_St_m[t,]/normalizer
    for(m in 1:M){
      index[[t]][[m]] <<- index[[t]][[m]]/normalizer
    }
  }
  e_St_minus_m_St_m
}
############################################################################
#                             M step updates                               #
############################################################################
fill_W <- function(A,t){
  startRows <- rep(1,M)
  for(m1 in 1:M){
    startCol <- 1
    for(m2 in 1:M){
      mat <- e_St_m_St_n[[t]][[m1]][[m2]]
      endCol <- startCol+dim[m2]-1
      startRow <- startRows[m2]
      endRow <- startRow+dim[m1]-1
      A[startRow:endRow,startCol:endCol] <- mat
      startCol <- startCol+dim[m2]
      startRows[m2] <- endRow+1
    }
  }
  A
}
W_new <- function(){
  total.dim <- sum(dim)
  r1 <- matrix(0,1,total.dim)
  r2 <- matrix(0,total.dim,total.dim)
  for(t in 1:T_){
    stacked <- NULL
    for(m in 1:M){
      stacked <- c(stacked,e_St[[t]][[m]])
    }
    stacked <- t(matrix(stacked))
    r1 <- r1+Y[t]%*%stacked
    cs <- matrix(0,total.dim,total.dim)
    cs <- fill_W(cs,t)
    r2 <- r2+cs
  }
  r1%*%ginv(r2)
}
PI_new <- function(){
  result <- PI
  for(m in 1:M){
    result[[m]] <- e_St[[1]][[m]]
  }
  result
}
P_new <- function(m){
  P_m <- P[[m]]
  sapply(1:dim[m],function(i){
    sapply(1:dim[m],function(j){
      P_m[i,j] <<- P_new_helper(m,i,j)
    })
  })
  P_m
}
P_new_helper <- function(m,i,j){
  result <- marginal <- 0
  sapply(2:T_,function(t){
    sapply(1:dim[m],function(prev){
      marginal <<- marginal+index[[t]][[m]][prev,j]
    })
    result <<- result+index[[t]][[m]][i,j]
  })
  result <- ifelse(marginal > 0,result/marginal,0.0)
  result
}
C_new <- function(){
  r1 <- (Y%*%Y)
  r2 <- 0
  sapply(1:T_,function(t){
    sapply(1:M,function(m){
      start_index <- starts[m]
      end_index <- ends[m]
      W_m <- W[,start_index:end_index]
      r2 <<- r2+W_m%*%matrix(e_St[[t]][[m]])%*%t(Y[[t]])
    })
  })
  (r1-r2)/T_
}
############################################################################
#                                 runEM                                    #
############################################################################
Y <- c(1.0,1.0,1.1,1.1,1.4,1.4,1.8,1.8,2.0,2.5,2.6,2.8,3.0,3.2,3.5,3.7)
T_ <- length(Y)
dim <- c(4,4)
M <- length(dim)
D <- 1
W <- matrix(1,nrow=D,ncol=sum(dim))
C <- matrix(1,nrow=D,ncol=D)
init_transition <- function(K){
  samples <- sapply(1:M,function(m){
    samp <- sample(100,K,replace=T)
    samp <- samp/sum(samp)
    samp
  })
  matrix(c(samples),K,K)
}
P <- vector("list",length=M)
sapply(1:M,function(m){
  P[[m]] <<- init_transition(dim[m])
})
PI <- vector("list",length=M)
sapply(1:M,function(m){
  samp <- sample(100,dim[m],replace=T)
  samp <- samp/sum(samp)
  PI[[m]] <<- samp
})
start.w <- 1; end.w <- dim[1]
starts <- c(start.w)
ends <- c(end.w)
for(m in 2:M){
  end.w <- end.w+dim[m]
  ends <- c(ends,end.w)
  start.w <- ends[m-1]+1
  starts <- c(starts,start.w)
}
for(i in 1:20){
  alpha.mat <- forward()
  beta.mat <- backward()
  gamma.mat <- compute_gamma(alpha.mat,beta.mat)
  index <- list()
  e_St <- get_e_St()
  e_St_m_St_n <- get_e_St_m_St_n()
  e_St_minus_m_St_m <- get_e_St_minus_m_St_m()
  print(paste('Completed E-step iteration',i))
  W <- W_new()
  print(paste('M-step: updated W iteration',i))
  PI <- PI_new()
  print(paste('M-step: updated PI iteration',i))
  for(m in 1:M){
    P[[m]] <- P_new(m)
  }
  print(paste('M-step: updated P iteration',i))
  C <- C_new()
  print(paste('M-step: updated C iteration',i))
}
############################################################################
#                         generate predictions                             #
############################################################################
pred <- sapply(1:T_,function(t){
  mu <- 0
  for(m in 1:M){
    start_index <- starts[m]
    end_index <- ends[m]
    W_m <- W[,start_index:end_index]
    mu <- mu+W_m%*%e_St[[t]][[m]]
  }
  mu[1]
})
############################################################################
#                        visualization of results                          #
############################################################################
cs.data <- data.frame(cbind(pred,Y))
corr <- cor(cs.data$pred,cs.data$Y)
cs <- ggplot(aes(pred,Y), data=cs.data) + geom_point(color="firebrick", size=3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "blue")) +
  ylab("observed INR") +
  xlab("predicted INR") +
  ggtitle(paste("predicted INR vs. observed INR on simulated data: r =",round(corr,3))) +
  theme(axis.title.y = element_text(face="bold", size=11)) +
  theme(axis.title.x = element_text(face="bold", size=11)) +
  theme(plot.title = element_text(color="black", face="bold", size=14, hjust=0.5)) +
  scale_x_continuous(breaks = round(seq(1.0, max(max(pred),max(Y)), by = 0.5),1)) +
  scale_y_continuous(breaks = round(seq(1.0, max(max(pred),max(Y)), by = 0.5),1))
cs