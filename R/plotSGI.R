
MAPK.screen.as.array <- function(data, Anno) {
  concentration = c(0,10,20,40,80,100,120,140)
  D = array(0.0,dim=c(ncol(data),6,6,8,8))
  ds = 0
  for (dataSlot in colnames(data)) {
    ds = ds + 1
    d <- data[[dataSlot]]

    Genes = c("Ptp69D","CG42327","drk","Ras85D","Gap1","CG13197")
    I = list()
    I[[1]] = 2:9
    I[[2]] = 10:17
    I[[3]] = 18:25
    I[[4]] = 26:33
    I[[5]] = 34:41
    I[[6]] = 42:49

    for (i in 1:6) {
      for (j in 1:6) {
        for (a in 1:8) {
          for (b in 1:8) {
            A = I[[i]][a]
            B = I[[j]][b]
            J = which((Anno$targetID1 == A) & (Anno$targetID2 == B))[1]
            D[ds,i,j,a,b] = d[J]
          }
        }
      }
    }
  }
  dimnames(D) <- list(channel=colnames(data), gene1=Genes, gene2=Genes, concentration1=concentration, concentration2=concentration)
  return(D)
}

MAPK.cv.TPS <- function(A, range.df = 6:56, channel=1) {
  concentration = c(0,10,20,40,80,100,120,140)
  Genes = dimnames(A)[[2]]
  res = list()
  d = c(dim(A)[1:3],length(concentration),length(concentration)) 
  res[["data"]] = A
  res[["smoothed"]] = array(0,dim=d)
  res[["smoothedfactor"]] = array(0,dim=d)
  res[["expected"]] = array(0,dim=d)
  res[["diff"]] = array(0,dim=d)
  res[["diffzscore"]] = array(0,dim=d)

  n.in = 8

  DF <- matrix(0, nrow=6, ncol=6)
  CVerror = list()
  CVerrorSD = list()
  for (i in 1:6) {
    CVerror[[i]] = list()
    CVerrorSD[[i]] = list()
    for (j in 1:6) {
      d.in = concentration
      z = as.vector(A[channel,i,j,,])
      name = sprintf("T = %s  Q = %s",Genes[i],Genes[j])
      x = rep(d.in,each=n.in)
      y = rep(d.in,times=n.in)
      data = data.frame(x=x,y=y)
      F = c(1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 1, 8, 1, 2, 3, 
4, 5, 6, 7, 3, 4, 5, 6, 7, 8, 1, 2, 7, 8, 1, 2, 3, 4, 5, 6, 4, 
5, 6, 7, 8, 1, 2, 3, 6, 7, 8, 1, 2, 3, 4, 5, 5, 6, 7, 8, 1, 2, 
3, 4)
      res = matrix(0.0,nrow=64,ncol=56)
      for (f in 1:8) {
        Train = which(F != f)
        Test = which(F == f)
        for (df in range.df) {
          cat("f=",f," df=",df,"\r")
          tps = Tps(data[Train,],z[Train],df=df)
          znew = predict(tps,data[Test,])
          res[Test,df] = (znew - z[Test])^2
        }
      }
      res2 = matrix(0.0,nrow=8,ncol=56)
      for (df in range.df) {
        res2[,df] = tapply(res[,df],F,mean)
      }
      CVerror[[i]][[j]] = apply(res2,2,mean)[range.df]
      CVerrorSD[[i]][[j]] = apply(res2,2,sd)[range.df]
      ## DF[i,j] = (range.df)[which.min(CVerror[[i]][[j]])]
      DF[i,j] = (range.df)[which(CVerror[[i]][[j]] <= (CVerror[[i]][[j]]+CVerrorSD[[i]][[j]])[which.min(CVerror[[i]][[j]])])[1]]
      cat("\n",name, ": DF = ", DF[i,j], "\n")
    }
  }
  dimnames(DF) = list(template = Genes, query = Genes)
  res = list(DF = DF, CVerror = CVerror, CVerrorSD = CVerrorSD)
  return(res)
}

MAPK.estimate.TPS <- function(A, DF, n.out=8, channel=1) {
  concentration = c(0,10,20,40,80,100,120,140)
  Genes = dimnames(A)[[2]]
  res = list()
  res[["data"]] = A
  d = c(dim(A)[1:3],n.out,n.out)
  res[["smoothed"]] = array(0,dim=d)
  res[["smoothedfactor"]] = array(0,dim=d)
  res[["expected"]] = array(0,dim=d)
  res[["diff"]] = array(0,dim=d)
  res[["diffzscore"]] = array(0,dim=d)
  res[["Genes"]] = Genes

  n.in = 8

  for (i in 1:6) {
    for (j in 1:6) {
      ## cat(channel," ",i," ",j,"\n")
      d.in = concentration
      d.out = seq(0,140, length.out=n.out)
      z = as.vector(A[channel,i,j,,])
      name = sprintf("T = %s  Q = %s",Genes[i],Genes[j])
      x = rep(d.in,each=n.in)
      y = rep(d.in,times=n.in)
      
      tps = Tps(data.frame(x=x,y=y),z,df=DF[i,j])
      xnew = rep(d.out,each=n.out)
      ynew = rep(d.out,times=n.out)
      znew = predict(tps,data.frame(x=xnew,y=ynew))
      zpred = predict(tps,data.frame(x=x,y=y))
      
      Z = matrix(z,nrow=n.in,ncol=n.in)
      Znew = matrix(znew,nrow=n.out,ncol=n.out)
      res[["smoothed"]][channel,i,j,,] = znew
      
      mainmain = Znew[1,1]
      res[["smoothedfactor"]][channel,i,j,,] = res[["smoothed"]][channel,i,j,,] - mainmain
      main1 = Znew[1,] - mainmain
      main2 = Znew[,1] - mainmain
      SD = median(abs(matrix(zpred,nrow=n.in,ncol=n.in) - Z))
      res[["expected"]][channel,i,j,,] = mainmain + rep(main1,each=n.out) + rep(main2,times=n.out)
      res[["diff"]][channel,i,j,,] = Znew - res[["expected"]][channel,i,j,,]
      res[["diffzscore"]][channel,i,j,,] = res[["diff"]][channel,i,j,,] / SD
    }
  }
  return(res)
}

MAPK.plot.TPS.all <- function(TPSmodel, range=c(-6,6), fill=c("cornflowerblue","cornflowerblue","black","#777700","yellow"), channel=1) {
  n.in = 8
  n.out = dim(TPSmodel$expected)[4]
  Genes = TPSmodel$Genes

  concentration = c("0","10","20","50","80","100","140","160")

  zcutoff=0
  
  Z = TPSmodel[["diffzscore"]][channel,,,,]
  for (i in 1:6) {
    for (j in 1:6) {
      Z[i,j,,] = t(Z[i,j,,])
    }
  }
  Z = Z[,,,rev(1:n.out)]
  for (i in 1:6) {
    Z[i,i,,] = 0.0
  }
  nt = n.out
  ct = 140 / (nt-1)
  D = data.frame(ZInteraction = as.vector(Z),
    Gene1 = factor(rep(Genes,times=6*n.out*n.out),levels=Genes),
    Gene2 = factor(rep(rep(Genes,each=6),times=n.out*n.out),levels=Genes),
    c1 = rep(rep(1:n.out,each=6*6),times=n.out),
    c2 = rep(1:n.out,each=6*6*n.out))
  D$ZInteraction[D$ZInteraction > range[2]] = range[2]
  D$ZInteraction[D$ZInteraction < range[1]] = range[1]
  lp = levelplot(ZInteraction ~ c1*c2 | Gene2+Gene1,
    data = D,
    at=seq(range[1],range[2],length.out=257),
    col.regions=colorRampPalette(fill)(256),
    col = "white",
    range=range,
    xlab="dsRNA concentration",
    ylab="dsRNA concentration",
    main="Genetic Interaction Surfaces",
    scales=list(x=list(at=c(0,50,100)/ct+1,labels=c("0","50","100"),cex=1,alternating=2),y=list(at=nt-(c(0,50,100)/ct),labels=c("0","50","100"),cex=1)),
    panel = function(at, contour, range, col.regions, labels, ...) {
      panel.levelplot(at = seq(range[1],range[2],length.out=257), contour=FALSE,col.regions=col.regions,...)
      panel.contourplot(at = seq(range[1],range[2],length.out=9), labels=list(labels=seq(range[1],range[2],length.out=9), cex=0.7, col="white"), contour=TRUE,col.regions=NA,...)
    }
    )
  lp
}

MAPK.plot.TPS.single <- function(gene1, gene2, TPSmodel, range=c(-6,6), fill=c("cornflowerblue","cornflowerblue","black","#777700","yellow"), channel=1) {
  n.in = 8
  n.out = dim(TPSmodel$expected)[4]

  Genes = TPSmodel$Genes
  i = match(gene1, Genes)
  j = match(gene2, Genes)

  zcutoff=0
  
  frange = range(TPSmodel[["smoothedfactor"]][channel,,,,],na.rm=TRUE)
  name = sprintf("T = %s  Q = %s",Genes[i],Genes[j])
  x = rep(seq(1,n.in,length.out=n.out),each=n.out)
  y = rep(seq(1,n.in,length.out=n.out),times=n.out)

  ZInteraction = TPSmodel[["diffzscore"]][channel,i,j,,]
  ZInteraction[ZInteraction > range[2]] = range[2]
  ZInteraction[ZInteraction < range[1]] = range[1]
  col = colorRampPalette(fill)(513)
  Col = matrix(col[257],nrow=n.out,ncol=n.out)
  Col[ZInteraction > zcutoff] = col[floor(257 + 256 * (ZInteraction[ZInteraction > zcutoff] - zcutoff) / (range[2] - zcutoff))]
  Col[ZInteraction < -zcutoff] = col[floor(257 + 256 * (ZInteraction[ZInteraction < -zcutoff] + zcutoff) / (range[2] - zcutoff))]

  par(mar = c(1, 1, 0, 0) + 0.1)
  cp=contourLines(seq(0,140, length.out=n.out),seq(0,140, length.out=n.out),t(ZInteraction)[,rev(1:ncol(ZInteraction))],levels=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))
  nt = dim(ZInteraction)[1]
  ct = 140 / (nt-1)
  lp = levelplot(t(ZInteraction)[,rev(1:ncol(ZInteraction))],
    at = seq(range[1],range[2],length.out=257),
    col.regions=colorRampPalette(fill)(256),
    colorkey=NULL,
    col = "white",
    range=range,
    scales=list(x=list(at=c(0,50,100)/ct+1,labels=c("0","50","100"),cex=2,alternating=2),y=list(at=nt-(c(0,50,100)/ct),labels=c("0","50","100"),cex=2)),
    xlab=list(label=sprintf("%s [ng]",Genes[j]),cex=2),
    ylab=list(label=sprintf("%s [ng]",Genes[i]),cex=2),
    panel = function(at, contour, range, col.regions, labels, ...) {
      z = ZInteraction[,]
      id = sapply(z, function(x) { sum(x > at) } )
      id[id < 1] = 1
      id[id > length(col.regions)] = length(col.regions)
      C = col2rgb(col.regions)/256
      Z = array(0.0,dim=c(dim(z)[1],dim(z)[2],3))
      Z[,,1] = C[1,id]
      Z[,,2] = C[2,id]
      Z[,,3] = C[3,id]
      grid.raster(Z, width = unit(1, "npc"), height = unit(1, "npc"))
      panel.contourplot(at = seq(range[1],range[2],length.out=9), labels=list(labels=seq(range[1],range[2],length.out=9), col="white"), contour=TRUE,col.regions=NA,...)
    }
    )
  lp
}

