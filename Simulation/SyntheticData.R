library(MASS)

C<-2; R<-2; n<-500; m<-20
lbd<-0; ubd<-2;
n.grids<-100
grids<-seq(lbd,ubd,length.out=n.grids+1)
grids.mid<-seq(lbd+(ubd-lbd)/n.grids/2,ubd-(ubd-lbd)/n.grids/2,by=(ubd-lbd)/n.grids)


Z <- list()
for(c in 1:C){
  Z[[c]] <- array(runif(51*2*R,-1,1), dim=c(51,2,R))
}
mu.x <- function(t, c, r){
  output <- 1
  for(k in 0:50){
    output <- output + (-1)^{k+1}*(k+1)^{-2}*(Z[[c]][k+1,1,r]*cos(k*pi*t)+Z[[c]][k+1,2,r]*sin(k*pi*t))
  }
  return(output)
}

p.x <- 6; p.y <- 3; p.z <- 4

Eigen.x <- function(t, k){
  t*0+sqrt(2)/2*(k==0)+sin(k*pi*t)*(k>0)
}
eta.x <- array(runif(p.x*C*R)*0.3,dim=c(C,R,p.x))

sigma.x <- array(0, dim=c(R*p.x,R*p.x,C))
for(c in 1:C){
  for(r1 in 1:(R-1)){
    for(r2 in (r1+1):R){
      for(k1 in 1:p.x){
        for(k2 in 1:p.x){
          e <- runif(1,-1,1)*sqrt(eta.x[c,r1,k1]*eta.x[c,r2,k2])
          sigma.x[(r1-1)*p.x+k1,(r2-1)*p.x+k2,c] <- e
          sigma.x[(r2-1)*p.x+k2,(r1-1)*p.x+k1,c] <- e
        }
      }
    }
  }
}
for(c in 1:C){
  for(r in 1:R){
    sigma.x[((r-1)*p.x+1):(r*p.x),((r-1)*p.x+1):(r*p.x),c] <- diag(eta.x[c,r,])
  }
}

Eigen.y <- list()
funy.1 <- function(t){
  sqrt(2)/2+t*0
}
funy.2 <- function(t){
  sqrt(2)*sin(pi*t)*(t<1)
}
funy.3 <- function(t){
  sqrt(2)*sin(pi*t)*(t>=1)
}

Eigen.y[[1]] <- funy.1
Eigen.y[[2]] <- funy.2
Eigen.y[[3]] <- funy.3

Eigen.z <- list()
funz.1 <- function(t){
  2*sin(4*pi*t)*(t<1/2)
}
funz.2 <- function(t){
  2*sin(4*pi*t)*(t>=1/2)*(t<1)
}
funz.3 <- function(t){
  2*sin(4*pi*t)*(t>=1)*(t<3/2)
}
funz.4 <- function(t){
  2*sin(4*pi*t)*(t>=3/2)
}

Eigen.z[[1]] <- funz.1
Eigen.z[[2]] <- funz.2
Eigen.z[[3]] <- funz.3
Eigen.z[[4]] <- funz.4


eta.y <- matrix(runif(R*p.y)*0.3,R,p.y)
eta.z <- matrix(runif(R*p.z)*0.3,R,p.z)*0.1

########################################
Sim_Poisson <- function(lbd,ubd,intensity,intensity.max=NULL){
  ###Get the rhomax if not specified
  if(is.null(intensity.max)){
    grids <- seq(lbd,ubd,l=1000)  
    intensity.max=max(intensity(grids))  
  }
  intensity.max <- max(intensity.max,15)
  ###Simulate total number of points
  count <- rpois(1,intensity.max*(ubd-lbd)) 
  if(count==0) process <- numeric(0) else{
    process_nothining <- sort(runif(count,lbd,ubd))  
    prob <- intensity(process_nothining)/intensity.max
    process <- process_nothining[prob>runif(count)]
  }
  
  process  
}

rho <- 0.8
xi.y <- array(rnorm(m*p.y*R,0,1),dim=c(m,p.y,R))
for(r in 1:R){
  for(j in 1:(m-1)){
    xi.y[j+1,1,r] <- xi.y[j,1,r]*rho+rnorm(1)*sqrt(1-rho^2)
  }
  xi.y[,,r] <- xi.y[,,r]%*%diag(sqrt(eta.y[r,]))
}
id <- timestamp <- days <- types <- labels <- c()
for(c in 1:C){
  xi.x <- array(mvrnorm(n,rep(0,R*p.x),sigma.x[,,c]), dim=c(n,p.x,R))
  
  for(i in 1:n){
    print(i)
    for(j in 1:m){
      for(r in 1:R){
        xi.z <- rnorm(p.z,0,1)*sqrt(eta.z[r,])
        intensity <- function(t){
          tmp.x <-tmp.y<-tmp.z<- 0
          for(k in 1:p.x)
            tmp.x <- tmp.x+Eigen.x(t,k-1)*xi.x[i,k,r]
          for(k in 1:p.y)
            tmp.y <- tmp.y+Eigen.y[[k]](t)*xi.y[j,k,r]
          for(k in 1:p.z)
            tmp.z <- tmp.z+Eigen.z[[k]](t)*xi.z[k]
          exp(mu.x(t,c,r)+tmp.x+tmp.y+tmp.z)
        }
        sample <- Sim_Poisson(lbd,ubd,intensity,intensity.max=NULL)
        if(length(sample)!=0){
          ct <- length(sample)
          id <- c(id, rep(i,ct))
          timestamp <- c(timestamp, sample)
          days <- c(days, rep(j,ct))
          types <- c(types, rep(r,ct))
        }
      }
    }
  } 
  labels <- c(labels, rep(c,n))
}

dataset = data.frame(accountNumber=id, timeStamps=timestamp, days=days, category=types)
write.csv(dataset, file = "dataset.csv", row.names = FALSE)

labelset <- data.frame(accountNumber=1:(C*n), label=labels)
write.csv(labelset, file = "labels.csv", row.names = FALSE)
