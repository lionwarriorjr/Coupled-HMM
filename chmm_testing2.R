require(MASS)
require(ggplot2)
Y <- c(1.0,1.0,1.1,1.1,1.4,1.4,1.8,1.8,2.0,2.5,2.6,2.8,3.0,3.2,3.5,3.7)
T_ <- length(Y)
############################################################################
#                       Forward/Backward Recursions                        #
############################################################################
compute_obs_prob <- function(y_t,states){
  num_states <- length(states)
  mu <- 0
  for(m in 1:num_states){
    W_m <- W[,((m-1)*K+1):(m*K)]
    one_hot <- rep(0,K)
    one_hot[states[m]] <- 1
    one_hot <- matrix(one_hot)
    mu <- mu+W_m%*%one_hot
  }
  result <- (det(C)^-0.5)%*%((2*pi)^(-D/2))%*%exp(-0.5%*%t(y_t-mu)%*%solve(C)%*%(y_t-mu))
  result <- result[1]
  result
}
forward <- function(){
  alphas <- array(0,c(K,K,T_))
  for(t in 1:T_){
    for(k1 in 1:K){
      for(k2 in 1:K){
        curr <- c(k1,k2)
        obs.prob <- compute_obs_prob(Y[t],curr)
        if(t == 1){
          alphas[k1,k2,t] <- obs.prob*PI[1,k1]*PI[2,k2]
        }else{
          recur <- forward_recursion(alphas,curr,NULL,t)
          alphas[k1,k2,t] <- obs.prob*recur
        }
      }
    }
  }
  alphas
}
forward_recursion <- function(alphas,curr,prev,t){
  if(length(prev) == M){
    result <- alphas[prev[1],prev[2],t-1]
  }else{
    index <- length(prev)+1
    transmat <- P[[index]]
    result <- 0
    for(kminus in 1:K){
      next_state <- curr[index]
      trans.prob <- transmat[kminus,next_state]
      prv <- c(prev,kminus)
      result <- result+trans.prob*forward_recursion(alphas,curr,prv,t)
    }
  }
  result
}
backward <- function(){
  betas <- array(0,c(K,K,T_))
  for(t in T_:1){
    for(kminus1 in 1:K){
      for(kminus2 in 1:K){
        prev <- c(kminus1,kminus2)
        betas[kminus1,kminus2,t] <- backward_recursion(betas,NULL,prev,t)
      }
    }
  }
  betas
}
backward_recursion <- function(betas,curr,prev,t){
  if(length(curr) == M){
    result <- compute_obs_prob(Y[t],curr)
    if(t < T_){
      result <- result*betas[curr[1],curr[2],t+1]
    }
  }else{
    index <- length(curr)+1
    transmat <- P[[index]]
    result <- 0
    for(k in 1:K){
      prev_state <- prev[index]
      trans.prob <- transmat[prev_state,k]
      nxt <- c(curr,k)
      result <- result+trans.prob*backward_recursion(betas,nxt,prev,t)
    }
  }
  result
}
compute_gamma <- function(alphas,betas){
  gammas <- array(0,c(K,K,T_))
  for(t in 1:T_){
    result <- 0
    for(k1 in 1:K){
      for(k2 in 1:K){
        gammas[k1,k2,t] <- alphas[k1,k2,t]*betas[k1,k2,t]
        result <- result+alphas[k1,k2,t]*betas[k1,k2,t]
      }
    }
    gammas[,,t] <- gammas[,,t]/result
  }
  gammas
}
compute_E_St <- function(t,m){
  result <- rep(0,K)
  for(k in 1:K){
    for(n in 1:K){
      one_hot <- rep(0,K)
      one_hot[k] <- 1
      if(m == 1){
        result <- result+one_hot*gamma.mat[k,n,t]
      }else{
        result <- result+one_hot*gamma.mat[n,k,t]
      }
    }
  }
  result
}
get_e_St <- function(){
  e_St <- array(0,c(T_,M,K))
  for(t in 1:T_){
    for(m in 1:M){
      e_St[t,m,] <- compute_E_St(t,m)
    }
  }
  e_St
}
compute_E_St_E_St_diff <- function(t){
  result <- matrix(0,K,K)
  for(k1 in 1:K){
    for(k2 in 1:K){
      o1 <- rep(0,K)
      o2 <- rep(0,K)
      o1[k1] <- 1
      o2[k2] <- 1
      result <- result+matrix(o1)%*%t(matrix(o2))*gamma.mat[k1,k2,t]
    }
  }
  result
}
compute_E_St_E_St_same <- function(t,m){
  result <- matrix(0,K,K)
  for(k in 1:K){
    for(n in 1:K){
      one_hot <- rep(0,K)
      one_hot[k] <- 1
      if(m == 1){
        result <- result+matrix(one_hot)%*%t(matrix(one_hot))*gamma.mat[k,n,t]
      }else{
        result <- result+matrix(one_hot)%*%t(matrix(one_hot))*gamma.mat[n,k,t]
      }
    }
  }
  result
}
compute_E_St_E_St <- function(t,m,n){
  if(m == n){
    return (compute_E_St_E_St_same(t,m))
  }else{
    return (compute_E_St_E_St_diff(t))
  }
}
get_e_St_m_St_n <- function(){
  e_St_m_St_n <- array(rep(0,T_*M^2*K^2),c(T_,M,M,K,K))
  for(t in 1:T_){
    for(m in 1:M){
      for(n in 1:M){
        e_St_m_St_n[t,m,n,,] <- compute_E_St_E_St(t,m,n)
      }
    }
  }
  e_St_m_St_n
}
compute_E_t_minus_t <- function(t,m){
  result <- 0
  if(m == 1){
    trans.curr <- P[[1]]
    trans.other <- P[[2]]
  }else{
    trans.curr <- P[[2]]
    trans.other <- P[[1]]
  }
  for(m1 in 1:K){
    for(m2 in 1:K){
      for(n in 1:K){
        for(r in 1:K){
          trans.prob <- trans.curr[m1,m2]*trans.other[n,r]
          if(m == 1){
            alpha <- alpha.mat[m1,n,t-1]
            beta <- beta.mat[m2,r,t]
            obs.prob <- compute_obs_prob(Y[t],c(m2,r))
          }else{
            alpha <- alpha.mat[n,m1,t-1]
            beta <- beta.mat[r,m2,t]
            obs.prob <- compute_obs_prob(Y[t],c(r,m2))
          }
          p <- alpha*trans.prob*obs.prob*beta
          result <- result+p
          index[m,m1,m2,t] <<- p
        }
      }
    }
  }
  result
}
get_e_St_minus_m_St_m <- function(){
  e_St_minus_m_St_m <- matrix(0,T_,M)
  for(t in 2:T_){
    normalizer <- 0
    for(m in 1:M){
      expect <- compute_E_t_minus_t(t,m)
      normalizer <- normalizer+expect
      e_St_minus_m_St_m[t,m] <- expect
    }
    normalizer <- normalizer[1]
    e_St_minus_m_St_m[t,] <- e_St_minus_m_St_m[t,]/normalizer
    index[,,,t] <- index[,,,t]/normalizer
  }
  e_St_minus_m_St_m
}
############################################################################
#                             M step updates                               #
############################################################################
W_new <- function(){
  r1 <- matrix(0,1,M*K)
  r2 <- matrix(0,M*K,M*K)
  for(t in 1:T_){
    stacked <- t(matrix(c(e_St[t,1,],e_St[t,2,])))
    r1 <- r1+Y[t]%*%stacked
    cs <- matrix(0,M*K,M*K)
    cs[1:K,1:K] <- e_St_m_St_n[t,1,1,,]
    cs[(K+1):(M*K),1:K] <- e_St_m_St_n[t,1,2,,]
    cs[1:K,(K+1):(M*K)] <- e_St_m_St_n[t,1,2,,]
    cs[(K+1):(M*K),(K+1):(M*K)] <- e_St_m_St_n[t,2,2,,]
    r2 <- r2+cs
  }
  r1%*%ginv(r2)
}
PI_new <- function(){
  result <- matrix(0,M,K)
  for(m in 1:M){
    result[m,] <- e_St[1,m,]
  }
  result
}
P_new <- function(m){
  P_m <- P[[m]]
  sapply(1:K,function(i){
    sapply(1:K,function(j){
      P_m[i,j] <<- P_new_helper(m,i,j)
    })
  })
  P_m
}
P_new_helper <- function(m,i,j){
  result <- marginal <- 0
  sapply(2:T_,function(t){
    marginal <<- marginal+sum(index[m,,j,t])
    result <<- result+index[m,i,j,t]
  })
  result <- ifelse(marginal > 0,result/marginal,0.0)
  result
}
C_new <- function(){
  r1 <- (Y%*%Y)
  r2 <- 0
  sapply(1:T_,function(t){
    sapply(1:M,function(m){
      W_m <- W[,((m-1)*K+1):(m*K)]
      r2 <<- r2+W_m%*%matrix(e_St[t,m,])%*%t(Y[[t]])
    })
  })
  (r1-r2)/T_
}
# runEM
assign("D", 1, envir=.GlobalEnv)
assign("M", 2, envir=.GlobalEnv)
assign("K", 4, envir=.GlobalEnv)
assign("W", matrix(1,nrow=D,ncol=M*K))
assign("C", matrix(1,nrow=D,ncol=D))
init_transition <- function(){
  sample1 <- sample(100,4,replace=T)
  sample1 <- sample1/sum(sample1)
  sample2 <- sample(100,4,replace=T)
  sample2 <- sample2/sum(sample2)
  sample3 <- sample(100,4,replace=T)
  sample3 <- sample3/sum(sample3)
  sample4 <- sample(100,4,replace=T)
  sample4 <- sample4/sum(sample4)
  result <- matrix(c(sample1,sample2,sample3,sample4),nrow=K,ncol=K)
  result
}
P_1 <- init_transition()
P_2 <- init_transition()
P <- list()
P[[1]] <- P_1
P[[2]] <- P_2
assign("P", P, envir=.GlobalEnv)
sample1 <- sample(100,4,replace=T)
sample1 <- sample1/sum(sample1)
sample2 <- sample(100,4,replace=T)
sample2 <- sample2/sum(sample2)
PI <- matrix(rbind(sample1,sample2),nrow=M,ncol=K)
assign("PI", PI, envir=.GlobalEnv)
for(i in 1:200){
  alpha.mat <- forward()
  beta.mat <- backward()
  gamma.mat <- compute_gamma(alpha.mat,beta.mat)
  assign("alpha.mat", alpha.mat, envir=.GlobalEnv)
  assign("beta.mat", beta.mat, envir=.GlobalEnv)
  assign("gamma.mat", gamma.mat, envir=.GlobalEnv)
  assign("index", array(rep(0,M*K^2*T_),c(M,K,K,T_)), envir=.GlobalEnv)
  e_St <- get_e_St()
  e_St_m_St_n <- get_e_St_m_St_n()
  e_St_minus_m_St_m <- get_e_St_minus_m_St_m()
  print(paste('Completed E-step iteration',i))
  assign("e_St", e_St, envir=.GlobalEnv)
  assign("e_St_m_St_n", e_St_m_St_n, envir=.GlobalEnv)
  assign("e_St_minus_m_St_m", e_St_minus_m_St_m, envir=.GlobalEnv)
  W <- W_new()
  print(paste('M-step: updated W iteration',i))
  print(W)
  PI <- PI_new()
  print(paste('M-step: updated PI iteration',i))
  print(PI)
  P[[1]] <- P_new(1)
  P[[2]] <- P_new(2)
  print(paste('M-step: updated P iteration',i))
  print(P)
  C <- C_new()
  print(paste('M-step: updated C iteration',i))
  print(C)
}
# latent state predictions
pred <- sapply(1:T_,function(t){
  mu <- 0
  for(m in 1:M){
    W_m <- W[,((m-1)*K+1):(m*K)]
    mu <- mu+W_m%*%e_St[t,m,]
  }
  mu[1]
})
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