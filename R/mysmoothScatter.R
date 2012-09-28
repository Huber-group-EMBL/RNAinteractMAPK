
MAPK.smooth.scatter <- function(x,y,n=75, nrpoints = 100, col="blue", pch=20, size=unit(0.3,"char"), cex=1.2, colramp=colorRampPalette(c("white","blue","green","yellow","red"))(256),
                            xlab="", ylab="", respect=FALSE) {
  if (respect) {
    rx = ry = range(c(x,y))
  } else {
    rx = range(x)
    ry = range(y)
  }
  D = kde2d(x,y,n=n, lims = c(rx, ry))

  xid = sapply(x, function(e) { sum(e > D$x) } )
  yid = sapply(y, function(e) { sum(e > D$y) } )
  xid[xid < 1] = 1
  yid[yid < 1] = 1
  xid[xid > n] = n
  yid[yid > n] = n
  d = D$z[xid+(yid-1)*n]
  Index = order(d)[1:nrpoints]
  
  C = col2rgb(colramp)
  nc = ncol(C)
  at = seq(0,1.0001*max(D$z),length.out=nc+1)

  colid = sapply(D$z, function(e) { sum(e > at) } )
  colid[colid < 1] = 1
  colid[colid > nc] = nc
  colid = matrix(colid, nrow=nrow(D$z),ncol=ncol(D$z))
  colid = t(colid)
  colid = colid[nrow(colid):1,]
  Img = array(0, dim=c(nrow(D$z),ncol(D$z),3))
  Img[,,1] = C[1,colid]/256
  Img[,,2] = C[2,colid]/256
  Img[,,3] = C[3,colid]/256

  vp=viewport(name="smoothscatter", layout = grid.layout(3,3,
                                      widths = unit(c(6,93,1),c("lines","null","lines")),
                                      heights = unit(c(1,93,6),c("lines","null","lines")),
                                      respect=TRUE)
                    )
  pushViewport(vp)
  pushViewport(viewport(name="plotsmoothscatter", layout.pos.row=2, layout.pos.col=2,clip=FALSE, xscale = rx, yscale = ry))
  grid.raster(Img,interpolate = TRUE)
  grid.xaxis(gp=gpar(cex=cex))
  grid.yaxis(gp=gpar(cex=cex))
  grid.rect()
  grid.lines(unit(rx,"native"),unit(c(0,0),"native"),gp=gpar(col="gray"))
  grid.lines(unit(c(0,0),"native"),unit(ry,"native"),gp=gpar(col="gray"))
  grid.points(unit(x[Index],"native"),unit(y[Index],"native"),pch=pch,size=size,gp=gpar(col=col))
  for (i in 1:length(xlab)) {
    grid.text(xlab[i],unit(0.5,"npc"),unit(-i-2,"lines"),gp=gpar(cex=cex))
  }
  for (i in 1:length(ylab)) {
    grid.text(ylab[i],unit(-3-length(ylab)+i,"lines"),unit(0.5,"npc"),rot=90,gp=gpar(cex=cex))
  }
  popViewport()
  popViewport()
}

