
# Helper function for VY
add.na.matrix <- function(A,B){
  dim1 <- dim(B)[1]
  dim2 <- dim(B)[2]

  for(i0 in 1:dim1){
    for(i1 in 1:dim2){
      if(is.na(B[i0,i1])==F){
        A[i0,i1] <- B[i0,i1]
      }

    }
  }

  A
}

#' makevy function for jags markov-model, currently breaks with bs items
#' @param ydat y data in item order, anx1, anx2, ... sed10
#' @param tord the tord of the ymat
#' @param p the number of items per factor
#' @param q the number of factors
#' @return the N by k vy matrix
#' @export
# Makes Concurrent variance matrix

makevy <- function(ydat, tord, p=10, q=10){ # p = # of items per factor 9 = # factor.
  vy <- as.matrix(ydat)
  N  <- dim(ydat)[1]

  ydat0 <- matrix(NA,p,q)

  for(i in 1:N){#i <- 1
    ydati <- ydat[i,]
    ordi  <- tord[i,]

    y02 <- matrix(NA,p,q)
    for(k0 in 1:(q*p)){#k0<-1
      wo <- which(ordi==k0)
      y01 <- rep(NA,q*p)
      y01[wo] <- unlist(ydati[wo])

      y02 <- add.na.matrix(y02,matrix(y01,p,q,byrow=F))

      vy[i,ordi==k0] <- mean(apply(y02,2,var,na.rm=T),na.rm=T)
    }


  }

  vy[is.nan(vy)] <- 0
  vy[is.na(vy)]  <- 0

  colnames(vy) <- colnames(ydat)

return(data.frame(vy))

}


#' Person fit index for jags markov model, currently breaks with bs items
#' @param ydat y data in item order, anx1, anx2, ... sed10
#' @param tord the tord of the ymat
#' @param p the number of items per factor
#' @param q the number of factors
#' @return the N by k pf matrix
#' @export

cdfun <- function(ydat,ord,p=10,q=10){

  modx <- '
  anx =~ anx1+anx2+anx3+anx4+anx5+anx6+anx7+anx8+anx9+anx10
  ang =~ ang1+ang2+ang3+ang4+ang5+ang6+ang7+ang8+ang9+ang10
  fre =~ fre1+fre2+fre3+fre4+fre5+fre6+fre7+fre8+fre9+fre10
  gre =~ greg1+greg2+greg3+greg4+greg5+greg6+greg7+greg8+greg9+greg10
  img =~ img1+img2+img3+img4+img5+img6+img7+img8+img9+img10
  art =~ art1+art2+art3+art4+art5+art6+art7+art8+art9+art10
  trs =~ trst1+trst2+trst3+trst4+trst5+trst6+trst7+trst8+trst9+trst10
  alt =~ alt1+alt2+alt3+alt4+alt5+alt6+alt7+alt8+alt9+alt10
  ord =~ ord1+ord2+ord3+ord4+ord5+ord6+ord7+ord8+ord9+ord10
  sed =~ sed1+sed2+sed3+sed4+sed5+sed6+sed7+sed8+sed9+sed10
  '
  upsilon <- ydat
  N  <- dim(ydat)[1]

  sem0 <- lavaan::sem(modx,ydat,meanstructure =T)
  sigma0 <- lavaan::fitted(sem0)

  for(i in 1:N){#i <- 1
    ydati <- ydat[i,]
    ordi  <- ord[i,]
    nom.ydati<- qp <- c()
    for(k0 in 1:(q*p)){#k0<-1
      nom.ydatik <- names(ydati)[ordi==k0]
      nom.ydati <- c(nom.ydati,nom.ydatik)
      # m-distances
      if(length(nom.ydati)==1){
        mu0 <- mean(ydat[,nom.ydati])
        sig0 <- sd(ydat[,nom.ydati])^2
      }else{
        mu0 <- apply(ydat[,nom.ydati],2,mean)
        sig0 <- cov(ydat[,nom.ydati])
      }
      md1 <- mahalanobis(ydati[,nom.ydati],center=mu0,cov=sig0)
      md2 <- mahalanobis(ydati[,nom.ydati],center=sigma0$mean[nom.ydati],cov=sigma0$cov[nom.ydati,nom.ydati])
      if(length(nom.ydati)==1){
        cd1 <- -0.5*(k0*log(2*pi)+log(sig0)+md1)
        cd2 <- -0.5*(k0*log(2*pi)+log(sigma0$cov[nom.ydati,nom.ydati])+md2)
      }else{
        cd1 <- -0.5*(k0*log(2*pi)+log(det(sig0))+md1)
        cd2 <- -0.5*(k0*log(2*pi)+log(det(sigma0$cov[nom.ydati,nom.ydati]))+md2)
      }
      upsilon[i,ordi==k0] <- -2*(cd2-cd1)
    }
  }
  return(upsilon)
}


p = 11

#xx = cdfun(ydat = ydat,ord = tord)


library(ggplot2)

d <- t(xx[p,iord[p,]])

ggplot() +
  geom_point(aes(x = 1:100,y = d)) +
  geom_vline(aes(xintercept = cnrmturk::treatment.timing[p]))

#' Person fit index for jags markov model, currently breaks with bs items
#' @param ydat y data in item order, anx1, anx2, ... sed10
#' @param tord the tord of the ymat
#' @param p the number of items per factor
#' @param q the number of factors
#' @return Returns the string function which provides the number of concurrent items which are in the middle category of 3 or 4
#' @export

stfun <- function(ydat,tord,p=10,q=10){
  N  <- dim(ydat)[1]
  st <- matrix(0,N,p*q)
  for(i in 1:N){#i <- 1
    ydati <- ydat[i,]
    ordi  <- tord[i,]
    y01 <- c()
    for(k0 in 1:(q*p)){#k0<-1
      wo <- which(ordi==k0)
      y01 <- c(y01,unlist(ydati[wo]))
      if(k0>5){
        #if(all(abs(y01[1:6+k0-6])==0.5)){
        #  vy[i,ordi==k0] <- 1
        #}
        st[i,ordi==k0] <- sum(abs(y01[1:6+k0-6])%in%c(3,4))
      }
    }
  }
  st[is.nan(st)] <- NA
  return(st)
}



