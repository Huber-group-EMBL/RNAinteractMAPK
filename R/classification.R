
MAPK.getCV <- function(X, y) {

  if (any(!is.finite(X))) {
    if (any(is.finite(X))) {
      M = mean(X[is.finite(X)])
    } else {
      M = 0
    }
    X[!is.finite(X)] = M
  }

  if (!is.factor(y)) {
    y = factor(y)
  }
  y = factor(y)

  ypred = y
  ypred[] = NA
  posterior = matrix(NA,nr=length(y),nc=nlevels(y))
  colnames(posterior) = levels(y)

  index = rep(0,ncol(X))
  weight = rep(0,ncol(X))
  row.names(posterior) = row.names(X)
  for (i in 1:length(y)) {
    model <-sda(X[-i,,drop=FALSE],y[-i],stop=-2)
    pp = predict(model, X[i,,drop=FALSE])
    ypred[i] <- pp$class
    posterior[i,] <- pp$posterior
    print(model$varIndex)
    index[model$varIndex] = index[model$varIndex] + 1
  }
  res <- list()
  res$posterior = posterior
  res$class = ypred
  res$index = index
  return(res)
}

MAPK.cv.classifier <- function(sgi, traingroups) {
  stopifnot( is( sgi, "RNAinteract" ) )

  format="targetMatrix"
  mixTemplateQuery = TRUE
  withoutgroups=c("neg", "pos")
  
  D <- getData(sgi, type="pi", format=format, mixTemplateQuery = mixTemplateQuery, withoutgroups=withoutgroups,drop=FALSE)
  y <- rep(NA, dim(D)[1])
  names(y) <- dimnames(D)[[1]]
  for (i in 1:length(traingroups)) {
    y[dimnames(D)[[1]] %in% traingroups[[i]]] <- names(traingroups)[i]
  }
  Dtrain = D[!is.na(y),,,,drop=FALSE]
  y = y[!is.na(y)]
  levels = unique(y)
  y = factor(y,levels = levels)

  res <- list()
  res$singleResults <- list()
  res$y = y

  for (s in getScreenNames(sgi)) {
    res$singleResults[[s]] <- list()
    for (c in getChannelNames(sgi)) {
      resCV <- MAPK.getCV(X=Dtrain[,,s,c], y=y)
      res$singleResults[[s]][[c]] <- list()
      res$singleResults[[s]][[c]]$CVposterior <- resCV$posterior
      res$singleResults[[s]][[c]]$CVclass <- resCV$class
      res$singleResults[[s]][[c]]$CVindex <- resCV$index
    }
  }

  CVposterior = res$singleResults[[1]][[1]]$CVposterior
  CVposterior[] = 0
  for (i in 1:length(res$singleResults)) {
    for (j in 1:length(res$singleResults[[i]])) {
      CVposterior[] = CVposterior[] + res$singleResults[[i]][[j]]$CVposterior
    }
  }
  for (i in 1:nrow(CVposterior)) {
    s=sum(CVposterior[i,])
    CVposterior[i,] = CVposterior[i,] / s
  }

  res$CVposterior = CVposterior

  invisible(res)
}

MAPK.predict.classification <- function(sgi, traingroups) {
  stopifnot( is( sgi, "RNAinteract" ) )

  format="targetMatrix"
  mixTemplateQuery = TRUE
  withoutgroups=c("neg", "pos")
  
  D <- getData(sgi, type="pi", format=format, mixTemplateQuery = mixTemplateQuery, withoutgroups=withoutgroups,drop=FALSE)
  y <- rep(NA, dim(D)[1])
  names(y) <- dimnames(D)[[1]]
  for (i in 1:length(traingroups)) {
    y[dimnames(D)[[1]] %in% traingroups[[i]]] <- names(traingroups)[i]
  }
  Dnew = D[is.na(y),,,,drop=FALSE]
  ynew = y[is.na(y)]
  Dtrain = D[!is.na(y),,,,drop=FALSE]
  y = y[!is.na(y)]
  levels = unique(y)
  y = factor(y,levels = levels)
  ynew = factor(ynew, levels = levels)

  res <- list()
  res$singleResults <- list()

  for (s in getScreenNames(sgi)) {
    res$singleResults[[s]] <- list()
    for (c in getChannelNames(sgi)) {
      model <-sda(Dtrain[,,s,c],y,stop=-2)
      posteriornew <- predict(model, Dnew[,,s,c])$posterior
      res$singleResults[[s]][[c]]$posteriornew <- posteriornew
    }
  }

  posteriornew = res$singleResults[[1]][[1]]$posteriornew
  posteriornew[] = 0
  for (i in 1:length(res$singleResults)) {
    for (j in 1:length(res$singleResults[[1]])) {
      posteriornew[] = posteriornew[] + res$singleResults[[i]][[j]]$posteriornew
    }
  }
  for (i in 1:nrow(posteriornew)) {
    s=sum(posteriornew[i,])
    posteriornew[i,] = posteriornew[i,] / s
  }

  res$posteriornew = posteriornew

  invisible(res)
}

