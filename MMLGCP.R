library(R.matlab)
library(spatstat)
library(compiler)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(fossil)

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

Sim_onelevel <- function(n,C,lbd,ubd,mu.x,Eigen.x,eta.x,prob=rep(1/C,C),n.grids=100){
  if(length(n)==1){
    id <- sample(1:C, size=n, prob=prob, replace=TRUE)
    N <- n
  }else{
    N <- length(n)
    id <- sample(n, size=N, replace=FALSE)
  }
  logIntensity <- Process <- list()
  l <- (ubd-lbd)/n.grids
  grids <- seq(lbd+l/2,ubd-l/2,l=n.grids)
  for(c in 1:C){
    idc <- which(id==c)
    nc <- length(idc)
    p.x <- length(eta.x[[c]])
    xi.x <- matrix(rnorm(p.x*nc,0,1),p.x,nc)
    xi.x <- (xi.x-rowMeans(xi.x))/sqrt(diag(var(t(xi.x))))
    xi.x <- t(xi.x)%*%diag(sqrt(eta.x[[c]]))
    
    for(i in 1:nc){
      intensity <- function(t){
        tmp.x <- mu.x[[c]](t)
        for(k in 1:p.x)
          tmp.x <- tmp.x+Eigen.x[[c]][[k]](t)*xi.x[i,k]
        exp(tmp.x)
      }
      logIntensity[[idc[i]]] <- list(log(intensity(grids)))
      s <- Sim_Poisson(lbd,ubd,intensity,intensity.max=NULL)
      Process[[idc[i]]] <- s
    }
  }
 
  list(Process=Process,id=id,xi.x=xi.x,logIntensity=logIntensity)
}


MSMPP <- function(Prcess, C, lbd, ubd, bwd, kern="epanechnikov",
                  n.grids = 100, n.sampling = 2000, max.iter = 50, reltol = 1e-4,
                  id.true = NULL){
  N <- length(Process)
  n <- sapply(Process,length)
  l <- (ubd-lbd)/n.grids
  grids <- seq(lbd+l/2,ubd-l/2,by=l)
  breaks <- seq(lbd,ubd,by=l)
  win <- as.owin(c(lbd,ubd,0,1))
  Kh <- function(t){
    dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
  }
  edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
  edge_sq <- outer(edge,edge,FUN="*")
  
  Items <- list(a=list(), b=matrix(0,N,n.grids), ids=list(), v=list(), y=list())
  for(i in 1:N){
    process <- unlist(Process[[i]])
    cp <- crosspairs(ppp(x=process,y=process*0,window=win),ppp(x=grids,y=grids*0,window=win),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), n.grids))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    A <- tmp2-Matrix::t(tmp)%*%tmp
    Items$a[[i]] <- A/edge_sq/N
    Items$b[i,] <- tmp1/edge/N
    s <- c(process, grids)
    Items$ids[[i]] <- cut(s, breaks = breaks, labels = FALSE)
    Items$v[[i]] <- hist(process, breaks = breaks, plot = FALSE)$counts + 1
    Items$v[[i]] <- l/Items$v[[i]][Items$ids[[i]]]
    Items$y[[i]] <- 1/Items$v[[i]] * c(rep(1,n[i]),rep(0,n.grids))
  }
  
  Pars <- list(pi = rep(1/C,C), mu.x = rep(list(0),C), A = rep(list(0),C),
               Eigen.x = rep(list(0),C), eta.x = rep(list(0),C), p.x = rep(0,C))
  #omega <- kmeans(Items$b, C)$cluster
  omega <- sample(rep(1:C,N), size = N, replace=FALSE)
  for(i in 1:N){
    Pars$mu.x[[omega[i]]] = Pars$mu.x[[omega[i]]] + Items$b[i,]
    Pars$A[[omega[i]]] = Pars$A[[omega[i]]] + Items$a[[i]]
  }
  for(c in 1:C){
    Pars$A[[c]] <- as.matrix(log(Pars$A[[c]]*Pars$pi[c]/(Pars$mu.x[[c]]%*%t(Pars$mu.x[[c]]))))
    Pars$mu.x[[c]] <- log(Pars$mu.x[[c]]/Pars$pi[c]) - diag(Pars$A[[c]])/2
    Eigen <- eigen(Pars$A[[c]])
    Pars$eta.x[[c]] <- Eigen$values[which(Eigen$values>(0.01*Eigen$values[1]))]*l
    Pars$p.x[c] <- length(Pars$eta.x[[c]])
    Pars$Eigen.x[[c]] <- Eigen$vectors[,1:Pars$p.x[c]]/sqrt(l)
  }
  
  loglikelihood <- c()
  randindex <- c()
  id.x <- rep(0,N)
  delta2 <- 0
  for(ii in 1:max.iter){
    print(ii)
    llh <- 0
    New <- list(pi = rep(1/C,C), mu.x = rep(list(0),C), A = rep(list(0),C),
                Eigen.x = rep(list(0),C), eta.x = rep(list(0),C), p.x = rep(0,C))
    for(i in 1:N){
      f.c = matrix(0, C, n.sampling)
      for(c in 1:C){
        xi.x <- matrix(rnorm(n.sampling*Pars$p.x[c],0,1),Pars$p.x[c],n.sampling)
        xi.x <- (xi.x-rowMeans(xi.x))/sqrt(diag(var(t(xi.x))))
        xi.x <- t(xi.x)%*%diag(sqrt(Pars$eta.x[[c]]))
        phi.x <- Pars$Eigen.x[[c]][Items$ids[[i]],]
        X <- Pars$mu.x[[c]][Items$ids[[i]]] + phi.x%*%t(xi.x)
        f.c[c,] = colSums(Items$v[[i]]*(Items$y[[i]]*X-exp(X)))
      }
      omega <- rowMeans(exp(f.c - max(f.c)))
      omega <- Pars$pi*omega
      llh <- llh + log(sum(omega)) + max(f.c)
      omega <- omega/sum(omega)
      id.x[i] <- which(omega==max(omega)) 
      New$pi <- New$pi + omega
      for(c in 1:C){
        New$mu.x[[c]] <- New$mu.x[[c]] + omega[c]*Items$b[i,]
        New$A[[c]] <- New$A[[c]] + omega[c]*Items$a[[i]]
      }
    }
    New$pi = New$pi/N
    delta1 <- 0
    for(c in 1:C){
      New$A[[c]] <- as.matrix(log(New$A[[c]]*New$pi[c]/(New$mu.x[[c]]%*%t(New$mu.x[[c]]))))
      New$mu.x[[c]] <- log(New$mu.x[[c]]/New$pi[c]) - diag(New$A[[c]])/2
      Eigen <- eigen(New$A[[c]])
      New$eta.x[[c]] <- Eigen$values[which(Eigen$values>(0.01*Eigen$values[1]))]*l
      New$p.x[c] <- length(New$eta.x[[c]])
      New$Eigen.x[[c]] <- Eigen$vectors[,1:New$p.x[c]]/sqrt(l)
      delta1 <- max(delta1, abs(Pars$mu.x[[c]]-New$mu.x[[c]]))
    }
    
    loglikelihood <- c(loglikelihood, llh)
    if(!is.null(id.true)){
      randindex <- c(randindex, rand.index(id.true, id.x))
    }
    print(loglikelihood)
    print(delta1)
    print(Pars$pi)
    print(randindex)
    
    r1 <- max(c(mu.x[[1]](grids),mu.x[[2]](grids)))
    r2 <- min(c(mu.x[[1]](grids),mu.x[[2]](grids)))
    plot(x=1:100,y=mu.x[[1]](grids),col='red',type='l',ylim=c(r2,r1))
    lines(x=1:100,y=mu.x[[2]](grids),col='red')
    #lines(x=1:100,y=mu.x[[3]](grids),col='red')
    lines(x=1:100,y=Pars$mu.x[[1]])
    lines(x=1:100,y=Pars$mu.x[[2]])
    #lines(x=1:100,y=Pars$mu.x[[3]])
    
    if(ii > 2){
      delta2 <- loglikelihood[ii] - loglikelihood[ii-1]
    }
    Pars <- New
  }
}


Plot_Covariance_Pattern <- function(R,lbd,ubd,n.dummy=100,limits=c(-0.2,0.2),n_color=2){
  l <- (ubd-lbd)/2/n.dummy
  coords <- expand.grid(seq(lbd+l,ubd-l,2*l),seq(lbd+l,ubd-l,2*l))
  
  output=data.frame('Covariance'=as.vector(R),'lon'=coords[,1],'lat'=coords[,2])
  ggplot(output, aes(lon, lat)) +
    geom_raster(aes(fill = Covariance)) +
    scale_fill_gradientn(colours = rainbow(n_color),limits=limits)
}


Lam<-function(t,v,Theta,xi){
  output=rep(0,length(t))
  for(i in 1:length(t)){
    output[i]=t(B_Spline(t[i]))%*%(v+Theta%*%xi)
  }
  return(exp(output))
}

f<-function(data,v,Theta,xi,T,nd=20){
  t_dummy=seq(T/(2*nd),T-1/(2*nd), length.out=nd)
  mu=sum(Lam(t_dummy,v,Theta,xi))*T/nd
  if(setequal(data,integer(0))){
    return(exp(-mu))
  }
  output=exp(-mu)*prod(Lam(data,v,Theta,xi))
  return(output)
}


Init_W<-function(n,c,B){
  output=array(0,c(n,c,B))
  for(i in 1:n){
    output[i,sample(1:c,1),1]=1
  }
  return(output)
}


Sample_W<-function(f_old,pi_hat){
  C=length(f_old)
  p=f_old*pi_hat
  p=p/sum(p)
  ic=sample(1:C,size=1,prob=p)
  output=rep(0,C)
  output[ic]=1
  return(output)
}

Sample_Xi<-function(f_old,xi,v,Theta,Sigma,T,data){
  pX=length(Sigma);xi_new=rep(0,pX)
  for(i in 1:pX){
    xi_new[i]=rnorm(1,0,sqrt(Sigma[i]))
  }
  alpha=min(1,f(data,v,Theta,xi_new,T)/f_old)
  if(runif(1)<alpha){
    return(xi_new)
  }else{
    return(xi)
  }
}

Update_Pi_Sigma<-function(w,xi){
  C=dim(w)[2];n=dim(xi)[1];pX=dim(xi)[2];B=dim(xi)[3]
  pi=rep(0,C);Sigma=array(0,dim=c(pX,pX,C));
  
  for(i in 1:n){
    for(b in 1:B){
      c=which(w[i,,b]==0)
      Sigma[,,c]=Sigma[,,c]+xi[i,,b]%*%t(xi[i,,b])
    }
  }
  
  for(c in 1:C){
    pi[c]=sum(as.vector(w[,c,]))
    Sigma[,,c]=Sigma[,,c]/pi[c]
  }
  
  output=list()
  output$pi=pi/(n*B)
  output$Sigma=Sigma
  return(output)
}




Update_Coefficients<-function(data,Theta,v,xi,w,T,
                              lr=0.01,n_it=50,nd=20,epsilon=0.01){
  shape=dim(w);n=shape[1];C=shape[2];B=shape[3];pX=dim(xi)[2]
  
  para=matrix(0,pX*(pX+1),C)
  for(c in 1:C){
    para[,c]=c(v[,c],as.vector(Theta[,,c]))
  }
  
  t_dummy=seq(T/(2*nd),T-1/(2*nd), length.out=nd)
  kappa=T/nd
  
  S_dummy=matrix(0,pX,0)
  for(t in t_dummy){
    S_dummy=cbind(S_dummy,B_Spline(t))
  }
  
  for(it in 1:n_it){
    para_old=para
    U=matrix(0,pX*(pX+1),C)
    D=array(0,c(pX*(pX+1),pX*(pX+1),C))
    
    for(i in 1:n){
      SumS=rep(0,pX)
      if(!setequal(data[[i]]$N,integer(0))){
        for(t in data[[i]]$N){
          SumS=SumS+B_Spline(t)
        }
      }
      
      for(b in 1:B){
        c=which(w[i,,b]==1)
        SumX=SumS;X_dummy=S_dummy
        for(l in 1:pX){
          SumX=c(SumX,xi[i,l,b]*SumS)
          X_dummy=rbind(X_dummy,xi[i,l,b]*S_dummy)
        }
        U[,c]=U[,c]-SumX
        intsty_dummy=kappa*exp(colSums(X_dummy*para[,c]))
        for(k in 1:nd){
          u=X_dummy[,k]*intsty_dummy[k]
          U[,c]=U[,c]+u
          D[,,c]=D[,,c]+u%*%t(X_dummy[,k])
        }
      }
    }
    for(c in 1:C){
      para[,c]=para[,c]-lr*solve(D[,,c],U[,c])
    }
    
    v_hat=para[1:pX,]
    Theta_hat=array(para[(pX+1):nrow(para),],dim=c(pX,pX,C))
    
    if(max(abs(para-para_old))<epsilon){
      break
    }
  }
  output=list(v_hat=v_hat,Theta_hat=Theta_hat,
              F=Cal_F(data,Theta_hat,v_hat,xi,w,T))
  return(output)
}


Cal_F<-function(data,Theta,v,xi,w,T,nd=20){
  shape=dim(w);n=shape[1];C=shape[2];B=shape[3];pX=dim(xi)[2]
  
  t_dummy=seq(T/(2*nd),T-1/(2*nd), length.out=nd)
  kappa=T/nd
  
  S_dummy=matrix(0,pX,0)
  for(t in t_dummy){
    S_dummy=cbind(S_dummy,B_Spline(t))
  }
  
  F=0
  for(i in 1:n){
    SumS=rep(0,pX)
    if(!setequal(data[[i]]$N,integer(0))){
      for(t in data[[i]]$N){
        SumS=SumS+B_Spline(t)
      }
    }
    for(b in 1:B){
      c=which(w[i,,b]==1)
      F=F+t(SumS)%*%(v[,c]+Theta[,,c]%*%xi[i,,b])
      for(k in 1:nd){
        F=F-kappa*exp(t(S_dummy[,k])%*%(v[,c]+Theta[,,c]%*%xi[i,,b]))
      }
    }
  }
  return(F/(n*B))
}

Cal_P<-function(w,p){
  shape=dim(w);n=shape[1];C=shape[2];B=shape[3]
  P=0
  for(i in 1:n){
    for(b in 1:B){
      c=which(w[i,,b]==1)
      P=P+log(p[c])
    }
  }
  return(P/(n*B))
}

Cal_G<-function(w,xi,Sigma){
  shape=dim(w);n=shape[1];C=shape[2];B=shape[3];pX=dim(xi)[2]
  G=0
  for(i in 1:n){
    for(b in 1:B){
      c=which(w[i,,b]==1)
      G=G-1/2*sum(log(Sigma[,c]))-1/2*sum(xi[i,,b]^2/Sigma[,c])
    }
  }
  return(G/(n*B))
}

Cal_L<-function(data,pi,Theta,v,Sigma,T,nd=20){
  n=length(data);pX=nrow(v);C=ncol(v)
  
  B=1000
  xi=array(0,dim=c(B,pX,C))
  for(c in 1:C){
    for(j in 1:pX){
      xi[,j,c]=rnorm(B,0,sqrt(Sigma[j,c]))
    }
  }
  
  t_dummy=seq(T/(2*nd),T-1/(2*nd), length.out=nd)
  kappa=T/nd
  
  S_dummy=matrix(0,pX,0)
  for(t in t_dummy){
    S_dummy=cbind(S_dummy,B_Spline(t))
  }
  
  l=0
  for(i in 1:n){
    L=0
    SumS=rep(0,pX)
    if(!setequal(data[[i]]$N,integer(0))){
      for(t in data[[i]]$N){
        SumS=SumS+B_Spline(t)
      }
    }
    for(b in 1:B){
      temp=rep(0,C)
      for(c in 1:C){
        temp[c]=temp[c]+t(SumS)%*%(v[,c]+Theta[,,c]%*%xi[b,,c])
        for(k in 1:nd){
          temp[c]=temp[c]-kappa*exp(t(S_dummy[,k])%*%(v[,c]+Theta[,,c]%*%xi[b,,c]))
        }
      }
      L=L+sum(exp(temp)*pi/B)
    }
    l=l+log(L)
  }
  return(l)
}

C<-2
R<-1
n.grids<-100
grids<-seq(lbd,ubd,length.out=n.grids+1)
grids.mid<-seq(lbd+(ubd-lbd)/n.grids/2,ubd-(ubd-lbd)/n.grids/2,by=(ubd-lbd)/n.grids)

mu<-list(); R.x<-list(); Q<-list()
alpha<-0
phi_k<-function(s,k){
  f=rep(0,n.grids)+(((k-1)< s)&(k>=s))*sqrt(2)*sin(pi*s)
  return(f)
}

eta.x <- Sigma.x <- list()
for(c in 1:C){
  mu[[c]] <- R.x[[c]] <- eta.x[[c]] <- Sigma.x[[c]] <- list()
  for(r in 1:R){
    mu[[c]][[r]] <- rep(1+runif(1,-0.5,0.5),n.grids)
    R.x[[c]][[r]] <- matrix(0,n.grids,n.grids)
    for(l in 1:50){
      mu[[c]][[r]] <- mu[[c]][[r]]+runif(1,-0.5,0.5)*(-1)^(l+1)*l^(-2)*cos(l*pi*grids.mid)+
        runif(1,-0.5,0.5)*(-1)^(l+1)*l^(-2)*sin(l*pi*grids.mid)
    }
    eta.x[[c]][[r]] <- runif(10,0,0.3)
    for(l in 1:10){
      R.x[[c]][[r]] <- R.x[[c]][[r]]+eta.x[[c]][[r]][l]*phi_k(grids.mid,l)%*%t(phi_k(grids.mid,l))
    }
  }
  if(R>1){
    Q[[c]] <- matrix(rep(list(),R*R),R,R)
    Sigma.x[[c]] <- matrix(rep(list(),R*R),R,R)
    for(r1 in 1:R){
      for(r2 in r1:R){
        Q[[c]][r1,r2][[1]] <- matrix(0,n.grids,n.grids)
        Sigma.x[[c]][r1,r2][[1]] <- matrix(runif(100,-1,1),10,10)
        for(k1 in 1:10){
          for(k2 in 1:10){
            Q[[c]][r1,r2][[1]] <- Q[[c]][r1,r2][[1]]+
              Sigma.x[[c]][r1,r2][[1]][k1,k2]*sqrt(eta.x[[c]][[r1]][[k1]]*eta.x[[c]][[r2]][[k2]])*phi_k(grids.mid,k1)%*%t(phi_k(grids.mid,k2))
          }
        }
      }
    }
  }
}


