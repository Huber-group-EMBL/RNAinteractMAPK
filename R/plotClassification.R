
# The code for the ternary plot is adapted from the function ternaryplot in the package vcd
# Author of the original code: David Meyer David.Meyer@R-project.org
# References: M. Friendly (2000), Visualizing Categorical Data. SAS Institute, Cary, NC.
# This code is specialized for the publication "Mapping Signalling Networks by RNAi ..." in Nat. Methods
# It is highly recommended to use the original code by David Meyer.
MAPK.ternary.plot <-function (x, scale = 1,
                            dimnames = NULL, dimnames_position = c("corner", "edge", "none"),
                            dimnames_color = "black", id = NULL, id_color = "black", 
                            coordinates = FALSE, grid = TRUE, grid_color = "darkgray",
                            labels = c("inside", "outside", "none"),
                            labels_color = "darkgray", border = "black", 
                            bg = "white", pch = 19, cex = 1, prop_size = FALSE, col = "red", 
                            main = "ternary plot", newpage = TRUE, pop = TRUE, ...)
{
    labels <- match.arg(labels)
    if (grid == TRUE) 
        grid <- "dotted"
    if (coordinates)
        id <- paste("(", round(x[, 1] * scale, 1), ",", round(x[, 
            2] * scale, 1), ",", round(x[, 3] * scale, 1), ")", 
            sep = "")
    dimnames_position <- match.arg(dimnames_position)
    if (is.null(dimnames) && dimnames_position != "none") 
        dimnames <- colnames(x)
    if (is.logical(prop_size) && prop_size) 
        prop_size <- 3
    if (ncol(x) != 3) 
        stop("Need a matrix with 3 columns")
    if (any(x < 0)) 
        stop("X must be non-negative")
    s <- rowSums(x)
    if (any(s <= 0)) 
        stop("each row of X must have a positive sum")
    x <- x/s
    top <- sqrt(3)/2
    if (newpage) 
        grid.newpage()
    xlim <- c(-0.03, 1.03)
    ylim <- c(-1, top)
    pushViewport(viewport(width = unit(1, "snpc"),name="ternary"))
    if (!is.null(main)) 
        grid.text(main, y = 0.93, gp = gpar(cex=cex, fontstyle = 1))
    pushViewport(viewport(width = 0.8, height = 0.8, xscale = xlim, 
        yscale = ylim, name = "plot"))
    eps <- 0.01
    grid.polygon(c(0, 0.5, 1), c(0, top, 0), gp = gpar(fill = bg, 
        col = border), ...)
    if (dimnames_position == "corner") {
        grid.text(x = c(0, 1, 0.5), y = c(-0.06, -0.06, top + 0.06),
                  label = dimnames,
                  gp = gpar(cex=cex,col=dimnames_color))
    }
    if (dimnames_position == "edge") {
        shift <- eps * if (labels == "outside") 
            8
        else 0
        grid.text(x = 0.25 - 2 * eps - shift, y = 0.5 * top + 
            shift, label = dimnames[2], rot = 60, gp = gpar(cex=cex,col = dimnames_color))
        grid.text(x = 0.75 + 3 * eps + shift, y = 0.5 * top + 
            shift, label = dimnames[1], rot = -60, gp = gpar(cex=cex,col = dimnames_color))
        grid.text(x = 0.5, y = -0.02 - shift, label = dimnames[3], 
            gp = gpar(cex=cex,col = dimnames_color))
    }
    if (is.character(grid)) 
        for (i in 1:4 * 0.2) {
            grid.lines(c(1 - i, (1 - i)/2), c(0, 1 - i) * top, 
                gp = gpar(lty = grid, col = grid_color))
            grid.lines(c(1 - i, 1 - i + i/2), c(0, i) * top, 
                gp = gpar(lty = grid, col = grid_color))
            grid.lines(c(i/2, 1 - i + i/2), c(i, i) * top, gp = gpar(lty = grid, 
                col = grid_color))
            if (labels == "inside") {
                grid.text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 * 
                  top, label = i * scale, gp = gpar(col = labels_color), 
                  rot = 120)
                grid.text(x = 1 - i + i/4 + eps, y = i/2 * top - 
                  eps, label = (1 - i) * scale, gp = gpar(col = labels_color), 
                  rot = -120)
                grid.text(x = 0.5, y = i * top + eps, label = i * 
                  scale, gp = gpar(col = labels_color))
            }
            if (labels == "outside") {
                grid.text(x = (1 - i)/2 - 6 * eps, y = (1 - i) * 
                  top, label = (1 - i) * scale, gp = gpar(cex=cex,col = labels_color))
                grid.text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 - 
                  i) * top + 5 * eps, label = i * scale, rot = -120, 
                  gp = gpar(cex=cex,col = labels_color))
                grid.text(x = i + eps, y = -0.05, label = (1 - 
                  i) * scale, vjust = 1, rot = 120, gp = gpar(cex=cex,col = labels_color))
            }
        }
    xp <- x[, 2] + x[, 3]/2
    yp <- x[, 3] * top
    size = unit(if (prop_size) 
        prop_size * (s/max(s))
    else cex, "lines")
    res = data.frame(xp=xp,yp=yp,col=col,size=convertWidth(size,"lines",valueOnly=TRUE),s=s/max(s),a=x[,1],b=x[,2],c=x[,3],label=row.names(x),stringsAsFactors=FALSE)
    if (!is.null(id)) 
        grid.text(x = xp, y = unit(yp - 0.015, "snpc") - 0.5 * 
            size, label = as.character(id), gp = gpar(col = id_color, 
            cex = cex))
    if (pop) 
        popViewport(2)
    else upViewport(2)
    return(res)
}

MAPK.getXY <- function(a, b, c, d = 0.03, shift = 0.2, toleft=TRUE) {
  if (toleft) {
    s=a+c
    a=a/s
    c=c/s
    a=a*(1-b)+b/2
    c=c*(1-b)+b/2
  }
  l=0.05-d
  r=1+d
  I=1:length(a)
  for (i in 1:length(a)) {
    M = matrix(0,nr=length(I),nc=2)
    M[,1] = a[I] - l
    M[,2] = r - a[I]
    
    m = rowMin(M)
    J = I[which.min(m)]
    if (a[J] - l < r - a[J]) {
      if (a[J] - l < d) {
        a[J] = l+d
      }
      l = a[J]
    } else {
      if (r - a[J] < d) {
        a[J] = r-d
      }
      r = a[J]
    }
    I = I[I != J]
  }
  c=1-a

  b = rep(-shift,length(a))
  a=a*0.9+0.05
  c=c*0.9+0.05
  s = a+c
  a = a * (1+shift) / s
  c = c * (1+shift) / s
  if (toleft) {
    x = b+c/2
    y = c*sqrt(3)/2
  } else {
    x = a+c/2
    y = c*sqrt(3)/2
  }
  return(list(x=x,y=y))
}

MAPK.plot.classification <- function(posterior,
                               classes = NULL, classnames = NULL,
                               col = "darkgray", y=NULL, classcol=NULL,
                               main="predicted classification probabilities",
                               pop=TRUE,
                               threshText = 0.3,
                               textToLeft = NULL,
                               textToRight = NULL) {
  n <- nrow(posterior) 
  if (is.null(classes)) {
    classes <- colnames(posterior)[1:3]
  }
  if (is.null(classnames)) {
    classnames <- classes
  }
  if (!is.null(y)) {
    if (is.null(col)) {
      col = rep("lightgray",n)
    } else {
      col = rep(col[1],n)
    }
    if (is.null(classcol)) { classcol <- rainbow(3) }
    if (is.null(names(classcol))) { names(classcol) <- classes }
    col[y==classes[1]] = classcol[classes[1]]
    col[y==classes[2]] = classcol[classes[2]]
    col[y==classes[3]] = classcol[classes[3]]
  }

  res = MAPK.ternary.plot(posterior[,classes],pch=20,
    labels="outside",
    col=col,cex=1.8,prop_size=TRUE,
    main=main,pop=FALSE,
    dimnames = classnames,
    dimnames_color = classcol)

  I = rev(order(res$size))
  res = res[I,]
  grid.points(res$xp, res$yp, pch = 21, gp = gpar(fill = res$col),
              default.units = "snpc", size = unit(0.7*res$size,"lines"),
              vp=vpPath("ternary","plot"))

    A = which(((res$b <= res$a) & (res$s > threshText) & !(res$label %in% textToRight)) | (res$label %in% textToLeft))
    B = which(((res$b > res$a) & (res$s > threshText) & !(res$label %in% textToLeft)) | (res$label %in% textToRight))
    if (length(A) > 0) {
      XY=MAPK.getXY(res$a[A],res$b[A],res$c[A],d=0.06,shift=0.1)
      grid.text(res$label[A],XY$x, XY$y,
                gp = gpar(cex=1.8,col = res$col[A]),
                default.units = "snpc",
                just=c("right","bottom"),
                vp=vpPath("ternary","plot"))
      grid.polyline(c(res$xp[A],XY$x),c(res$yp[A],XY$y),id=rep(1:length(A),2),
                    vp=vpPath("ternary","plot"))
    }
    if (length(B) > 0) { 
      XY=MAPK.getXY(res$b[B],res$a[B],res$c[B],d=0.06,shift=0.1,toleft=FALSE)
      grid.text(res$label[B],XY$x, XY$y,
                gp = gpar(cex=1.8,col = res$col[B]),
                default.units = "snpc",
                just=c("left","bottom"),
                vp=vpPath("ternary","plot"))
      grid.polyline(c(res$xp[B],XY$x),c(res$yp[B],XY$y),id=rep(1:length(B),2),
                    vp=vpPath("ternary","plot"))
    }
}

