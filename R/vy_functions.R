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

#' Make moving variance function correct
#' @param ymat y data in item order, anx1, anx2, ... sed10
#' @param iord the iord of the ymat
#' @param tord the tord of the ymat
#' @param startLag how many items to skip before starting moving variance calculation, defaults to 4, should probably be higher.
#' @return Returns the moving variance matrix in item presentation order.
#' @export
makevy2 <- function(ymat,iord,tord,startLag = 4 ){
  skip <- startLag + 1
  holder <- list()

  for(j in 1:dim(ymat)[1]){
    iy <-  ymat[j,iord[j,]]
    vyi <- rep(NA,length(iy))
    for(i in skip:length(iy)){
      if(i == skip){
        vyi[skip]<- var(as.numeric(iy[1:i]))
      }else{
        vyi[i]<- var(as.numeric(iy[1:i]))
      }
    }
    holder[[j]]<- vyi[tord[j,]]
  }
  # This part scales the entire matrix while removing the starlag and adding it back as 0 after
  vyx <- do.call("rbind",holder)
  VYZ <- scale(c(vyx))
  vy <- matrix(VYZ, nrow = dim(vyx)[1])
  vy[is.na(vy)] <-  0
  return(vy)
}


#' makevy function for jags markov-model, currently breaks with bs items
#' @param ydat y data in item order, anx1, anx2, ... sed10
#' @param tord the tord of the ymat
#' @param p the number of items per factor
#' @param q the number of factors
#' @return the N by k vy matrix
#' @export
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
cdfun <- function(ydat,tord,p=10,q=10,static = FALSE){

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
    ordi  <- tord[i,]
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

  if(static == FALSE){
  # This part scales the entire matrix
  cd <- (upsilon - mean(as.matrix(upsilon)))/sd(as.matrix(upsilon))
  return(cd)
  }else if (static == TRUE){
  return(cd)
  }
}


#' Repeated string values for jags markov model, currently breaks with bs items
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


#' Extracts computes and formats Jags model to plot format data
#' @param RDSjags JAGS output model as list with inputs
#' @param cd Defaults to TRUE, was the CD function included in the markov model
#' @return Returns the plotting DF of a post-proceed JAGS model
#' @export
getCprobs <- function(RDSjags,cd = TRUE){
out <- RDSjags
SumDat <- out$fit$BUGSoutput$summary
cond <- data.frame("cond"= out$Cond,"ID2" = names(out$Cond))
Treatr <- data.frame("Treat"= out$Treat)
Treat <- Treatr
for(i in 1:nrow(Treatr)){
  if(is.na(Treatr[i,])){
    Treat[i,] <- 100
  }else{
    Treat[i,] <- Treatr[i,]
  }
}

Treat <- data.frame(cond,Treat)


N <- out$input$N
Y <- out$input$y
Treat$ID <- 1:N
Treat$ID2 <-  row.names(Treat)
K <- 100 # K item
index <- expand.grid(1:N,1:K)

C <- data.frame("Est_class" = SumDat[paste0("C[",index[,1],",",index[,2],"]"),"mean"],
                "ID" = index[,1],
                "MO" = index[,2])

matC <- matrix(NA,N,K)
for(i in 1:N){
  for(j in 1:K){
    matC[i,j] <- C[C$ID == i & C$MO == j,"Est_class"]
  }
}

iord <- out$input$ord
holder <- list()
for(i in 1:N){
  holder[[i]] <- matC[i,iord[i,]]
}

C <- do.call(rbind,holder)
holder <- list()
for(i in 1:N){
  holder[[i]] <- data.frame("class" = C[i,],
                            "t" = 1:K,
                            "ID" = rep(i,K),
                            "ID2" = rep(Treat[i,"ID2"],K))
}

CC <- do.call("rbind",holder)

PC1 <- data.frame("PC_prob" = SumDat[paste0("PC1[",index[,1],",",index[,2],",1]"),"mean"],
                  "ID" = index[,1],
                  "MO" = index[,2])

matP <- matrix(NA,N,K)
for(i in 1:N){
  for(j in 1:K){
    matP[i,j] <- PC1[PC1$ID == i & PC1$MO == j,"PC_prob"]
  }
}

iord <- out$input$ord
holder <- list()
for(i in 1:N){
  holder[[i]] <- matP[i,iord[i,]]
}

PC1 <- do.call(rbind,holder)
holder <- list()
for(i in 1:N){
  holder[[i]] <- data.frame("classP" = PC1[i,],
                            "t" = 1:K,
                            "ID" = rep(i,K),
                            "ID2" = rep(Treat[i,"ID2"],K))
}

PCC1 <- do.call("rbind",holder)

# Now formatting the variance and cd function

matVY <- out$input$vy
holder <- list()
for(i in 1:N){
  holder[[i]] <- matVY[i,iord[i,]]
}

VYY1 <- do.call("rbind",holder)

holder <- list()
for(i in 1:N){
  holder[[i]] <- data.frame("VY1" = VYY1[i,],
                            "t" = 1:K,
                            "ID" = rep(i,K),
                            "ID2" = rep(Treat[i,"ID2"],K))
}

VYX <- do.call("rbind",holder)

if(cd == TRUE){
  matCD <- as.matrix(out$input$cd)
  holder <- list()
  for(i in 1:N){
    holder[[i]] <- matCD[i,iord[i,]]
  }

  CDD1 <- as.matrix(do.call("rbind",holder))

  holder <- list()
  for(i in 1:N){
    holder[[i]] <- data.frame("CD1" = as.vector(CDD1[i,]),
                              "t" = 1:K,
                              "ID" = rep(i,K),
                              "ID2" = rep(Treat[i,"ID2"],K))
  }

  CDX <- do.call("rbind",holder)

}
# Observed data

holder <- list()
for(i in 1:N){
  holder[[i]] <- Y[i,iord[i,]]
}

Y1 <- as.matrix(do.call("rbind",holder))

holder <- list()
for(i in 1:N){
  holder[[i]] <- data.frame("Yobs" = Y1[i,],
                            "t" = 1:K,
                            "ID" = rep(i,K),
                            "ID2" = rep(Treat[i,"ID2"],K))
}

YY1 <- do.call("rbind",holder)

# Putting it together in a plotting data frame
pdat <- merge(x = merge(x = CC,
                        y = PCC1,
                        by = c("ID2","t")),
              y = Treat,
              by = "ID2")

pdat2 <- merge(x = pdat,y = YY1, by = c("ID2","t"))

pdat3 <- merge(x = pdat2,y = VYX, by = c("ID2","t"))[,c("ID2",
                                                        "t",
                                                        "VY1",
                                                        "Yobs",
                                                        "class",
                                                        "classP",
                                                        "Treat",
                                                        "cond")]

if(cd == TRUE){
  pdat3 <- merge(x = pdat3,y = CDX, by = c("ID2","t"))[,c("ID2",
                                                          "t",
                                                          "VY1",
                                                          "CD1",
                                                          "Yobs",
                                                          "class",
                                                          "classP",
                                                          "Treat",
                                                          "cond")]
}



pdat3$N <- rep(N,nrow(pdat3))

# before or after
index <- unique(pdat3$ID2)
holder <- list()
for(i in 1:length(index)){
  casei <- pdat3[pdat3$ID2 == index[i],]
  TT <- casei$Treat[1]
  casei$BIN <- ifelse(casei$t < TT,
                      "before",
                      "after")
  holder[[i]] <- casei
}
pdat4 <- do.call("rbind",holder)

return(pdat4)
}

# getCprobs <- function(RDSjags, cd = TRUE){
#   out <- RDSjags
#   SumDat <- out$fit$BUGSoutput$summary
#   Treat <- data.frame("Treat"= out$Treat)
#   N <- out$input$N
#   Y <- out$input$y
#   Treat$ID <- 1:N
#   Treat$ID2 <-  row.names(Treat)
#   K <- 100 # K item
#   index <- expand.grid(1:N,1:K)
#
#   C <- data.frame("Est_class" = SumDat[paste0("C[",index[,1],",",index[,2],"]"),"mean"],
#                   "ID" = index[,1],
#                   "MO" = index[,2])
#
#   matC <- matrix(NA,N,K)
#   for(i in 1:N){
#     for(j in 1:K){
#       matC[i,j] <- C[C$ID == i & C$MO == j,"Est_class"]
#     }
#   }
#
#   iord <- out$input$ord
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- matC[i,iord[i,]]
#   }
#
#   C <- do.call(rbind,holder)
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- data.frame("class" = C[i,],
#                               "t" = 1:K,
#                               "ID" = rep(i,K),
#                               "ID2" = rep(Treat[i,"ID2"],K))
#   }
#
#   CC <- do.call("rbind",holder)
#
#   PC1 <- data.frame("PC_prob" = SumDat[paste0("PC1[",index[,1],",",index[,2],",1]"),"mean"],
#                     "ID" = index[,1],
#                     "MO" = index[,2])
#
#   matP <- matrix(NA,N,K)
#   for(i in 1:N){
#     for(j in 1:K){
#       matP[i,j] <- PC1[PC1$ID == i & PC1$MO == j,"PC_prob"]
#     }
#   }
#
#   iord <- out$input$ord
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- matP[i,iord[i,]]
#   }
#
#   PC1 <- do.call(rbind,holder)
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- data.frame("classP" = PC1[i,],
#                               "t" = 1:K,
#                               "ID" = rep(i,K),
#                               "ID2" = rep(Treat[i,"ID2"],K))
#   }
#
#   PCC1 <- do.call("rbind",holder)
#
#   # Now formatting the variance and cd function
#
#   matVY <- out$input$vy
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- matVY[i,iord[i,]]
#   }
#
#   VYY1 <- do.call("rbind",holder)
#
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- data.frame("VY1" = VYY1[i,],
#                               "t" = 1:K,
#                               "ID" = rep(i,K),
#                               "ID2" = rep(Treat[i,"ID2"],K))
#   }
#
#   VYX <- do.call("rbind",holder)
#
#   if(cd == TRUE){
#     matCD <- as.matrix(out$input$cd)
#     holder <- list()
#     for(i in 1:N){
#       holder[[i]] <- matCD[i,iord[i,]]
#     }
#
#     CDD1 <- as.matrix(do.call("rbind",holder))
#
#     holder <- list()
#     for(i in 1:N){
#       holder[[i]] <- data.frame("CD1" = as.vector(CDD1[i,]),
#                                 "t" = 1:K,
#                                 "ID" = rep(i,K),
#                                 "ID2" = rep(Treat[i,"ID2"],K))
#     }
#
#     CDX <- do.call("rbind",holder)
#
#   }
#   # Observed data
#
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- Y[i,iord[i,]]
#   }
#
#   Y1 <- as.matrix(do.call("rbind",holder))
#
#   holder <- list()
#   for(i in 1:N){
#     holder[[i]] <- data.frame("Yobs" = Y1[i,],
#                               "t" = 1:K,
#                               "ID" = rep(i,K),
#                               "ID2" = rep(Treat[i,"ID2"],K))
#   }
#
#   YY1 <- do.call("rbind",holder)
#
#   # Putting it together in a plotting data frame
#   pdat <- merge(x = merge(x = CC,
#                           y = PCC1,
#                           by = c("ID2","t")),
#                 y = Treat,
#                 by = "ID2")
#
#   pdat2 <- merge(x = pdat,y = YY1, by = c("ID2","t"))
#
#
#   pdat3 <- merge(x = pdat2,y = VYX, by = c("ID2","t"))[,c("ID2",
#                                                           "t",
#                                                           "VY1",
#                                                           "Yobs",
#                                                           "class",
#                                                           "classP",
#                                                           "Treat")]
#
#   if(cd == TRUE){
#     pdat3 <- merge(x = pdat3,y = CDX, by = c("ID2","t"))[,c("ID2",
#                                                             "t",
#                                                             "VY1",
#                                                             "CD1",
#                                                             "Yobs",
#                                                             "class",
#                                                             "classP",
#                                                             "Treat")]
#   }
#
#
#
#   pdat3$N <- rep(N,nrow(pdat3))
#
#   # before or after
#   index <- unique(pdat3$ID2)
#   holder <- list()
#   for(i in 1:length(index)){
#     casei <- pdat3[pdat3$ID2 == index[i],]
#     TT <- casei$Treat[1]
#     casei$BIN <- ifelse(casei$t < TT,
#                         "before",
#                         "after")
#     holder[[i]] <- casei
#   }
#   pdat4 <- do.call("rbind",holder)
#
#   return(pdat4)
# }
