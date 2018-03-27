require(MASS)
require(ggplot2)

# parameters of this model are A,C,R,Q,PI_1,V_1
# example sequence 1
Y <- c(1.0,1.0,1.1,1.1,1.4,1.4,1.8,1.8,2.0,2.5,2.6,2.8,3.0,3.2,3.5,3.7)
T_ <- length(Y)
y <- lapply(Y,function(y){t(y)})

# Kalman filter forward recursions in E step
forward <- function(t){
  x <- V <- NULL
  if(t == 1){
    K <- V1%*%t(C)%*%solve(C%*%V1%*%t(C)+R)
    x <- PI+K%*%(y[[1]]-C%*%PI)
    V <- V1-K%*%C%*%V1
    Vstep <<- append(Vstep,V1)
  }else{
    prev <- forward(t-1)
    x_prev <- prev[[1]]
    V_prev <- prev[[2]]
    V_step <- A%*%V_prev%*%t(A)+Q
    K <- V_step%*%t(C)%*%solve(C%*%V_step%*%t(C)+R)
    x_step <- A%*%x_prev
    x <- x_step+K%*%(y[[t]]-C%*%x_step)
    V <- V_step-K%*%C%*%V_step
    Vstep <<- append(Vstep,V_step)
  }
  xinf <<- append(xinf,x)
  Vinf <<- append(Vinf,V)
  return (list(x,V))
}

# Kalman smoothing backward recursions in E step
backward1 <- function(t){
  if(t == T_){
    xsmooth <<- append(xsmooth,xinf[[T_]])
    Vsmooth <<- append(Vsmooth,Vinf[[T_]])
    return (list(xinf[[T_]],Vinf[[T_]]))
  }else{
    J <- Vinf[[t]]%*%t(A)%*%solve(Vstep[[t+1]])
    Jsmooth <<- append(Jsmooth,J)
    nxt <- backward1(t+1)
    x_next <- nxt[[1]]
    V_next <- nxt[[2]]
    x <- xinf[[t]]+J%*%(x_next-A%*%xinf[[t]])
    V <- Vinf[[t]]+J%*%(V_next-Vstep[[t+1]])%*%t(J)
    xsmooth <<- append(xsmooth,x)
    Vsmooth <<- append(Vsmooth,V)
    return (list(x,V))
  }
}

# Kalman smoothing backward recursions in E step
backward2 <- function(t){
  V <- NULL
  if(t+1 == T_){
    K <- Vstep[[T_]]%*%t(C)%*%solve(C%*%Vstep[[T_]]%*%t(C)+R)
    KC <- K%*%C
    V <- (diag(nrow(KC))-KC)%*%A%*%Vinf[[T_-1]]
  }else{
    V_next <- backward2(t+1)
    V <- Vinf[[t+1]]%*%t(Jsmooth[[t]])+Jsmooth[[t+1]]%*%(V_next-A%*%Vinf[[t+1]])%*%t(Jsmooth[[t]])
  }
  Vdouble <<- append(Vdouble,V)
  return (V)
}

C_new <- function(){
  dot1 <- 0
  sapply(1:T_,function(t){
    dot1 <<- dot1 + y[[t]]%*%t(xsmooth[[t]])
  })
  dot2 <- solve(Reduce('+',P))
  dot1%*%dot2
}

R_new <- function(){
  result <- 0
  sapply(1:T_,function(t){
    result <<- result + y[[t]]%*%t(y[[t]])-C%*%xsmooth[[t]]%*%t(y[[t]])
  })
  result/T_
}

A_new <- function(){
  p1 <- Reduce('+',Pdouble)
  p2 <- Reduce('+',P)-P[[length(P)]]
  p1%*%solve(p2)
}

Q_new <- function(){
  p1 <- Reduce('+',P)-P[[1]]
  p2 <- Reduce('+',Pdouble)
  result <- p1-A%*%p2
  result/(T_-1)
}

assign("A", matrix(1.0), envir = .GlobalEnv)
assign("C", matrix(1.0), envir = .GlobalEnv)
assign("Q", matrix(0.1), envir = .GlobalEnv)
assign("R", matrix(0.1), envir = .GlobalEnv)
assign("PI", matrix(1.0), envir = .GlobalEnv)
assign("V1", matrix(0.1), envir = .GlobalEnv)
# runEM
for(i in 1:1000){
  print(i)
  assign("xinf", list(), envir = .GlobalEnv)
  assign("Vinf", list(), envir = .GlobalEnv)
  assign("Vstep", list(), envir = .GlobalEnv)
  assign("xsmooth", list(), envir = .GlobalEnv)
  assign("Vsmooth", list(), envir = .GlobalEnv)
  assign("Jsmooth", list(), envir = .GlobalEnv)
  assign("Vdouble", list(), envir = .GlobalEnv)
  assign("P", list(), envir = .GlobalEnv)
  assign("Pdouble", list(), envir = .GlobalEnv)
  forward(T_)
  backward1(1)
  backward2(1)
  for(t in 1:T_){
    P <- append(P,Vsmooth[[t]]+xsmooth[[t]]%*%t(xsmooth[[t]]))
  }
  for(t in 1:(T_-1)){
    p <- Vdouble[[t]]+xsmooth[[t+1]]%*%t(xsmooth[[t]])
    Pdouble <- append(Pdouble,p)
  }
  C <- C_new()
  R <- R_new()
  A <- A_new()
  Q <- Q_new()
  PI <- xsmooth[[1]]
  V1 <- P[[1]]-xsmooth[[1]]%*%t(xsmooth[[1]])
}
pred <- NULL
# latent state predictions
sapply(1:T_,function(t){
  pred <<- c(pred,C%*%(A%*%xsmooth[[t]]+rnorm(1,0,Q))+rnorm(1,0,R))
})
# evaluation of latent states on selected process with scalar output
sapply(1:T_,function(t){
  MASS::ginv(C)%*%y[[t]]
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