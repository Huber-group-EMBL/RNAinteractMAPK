
MAPK.plot.heatmap.raster <- function(X, subset=NULL, hc.row = NULL, hc.col = NULL, pi.max=NULL) {
  dd.row = as.dendrogram(hc.row)
  dd.col = as.dendrogram(hc.col)
  ord.row = order.dendrogram(dd.row)
  ord.col = order.dendrogram(dd.col)

  rbcol <- colorRampPalette(c("cornflowerblue","black","yellow"))(513)
  C = t(col2rgb(rbcol))/255

  pi.min = mad(X,center=0.0)
  if (is.null(pi.max)) {
    pi.max = 3*pi.min
  }
  pi = abs(upperTriangle(X))
  pi = pi[(pi >= pi.min)]
  quant = quantile(pi, probs = seq(0,1,length.out=257))
  at = c(-quant,quant)

  X = t(X[ord.row, ord.col])
  diag(X) = 0.0
  if (!is.null(subset)) {
    X = X[subset,]
  }
  colid = sapply(X, function(e) { sum(e > at) } )
  colid[colid < 1] = 1
  colid[colid > length(rbcol)] = length(rbcol)
  colid = matrix(colid, nrow=nrow(X),ncol=ncol(X))
  colid = colid[nrow(colid):1,]
  Img = array(0, dim=c(nrow(X),ncol(X),3))
  Img[,,1] = C[colid,1]
  Img[,,2] = C[colid,2]
  Img[,,3] = C[colid,3]

  N = nrow(X)
  M = ncol(X)
  K = 3
  Img2 = array(0, dim=c(K*N,K*M,3))
  for (i in (1:K-1)) {
    for (j in (1:K-1)) {
      for (k in 1:3) {
        Img2[1:N*K-i,1:M*K-j,k] = Img[1:N,1:M,k]
      }
    }
  }


  grid.raster(Img2,interpolate=FALSE)
  grid.rect()
}

