
MAPK.report.gene.lists.paper <- function(sgi, sgilimma, sgi3T2, screen = "mean") {
  stopifnot( is( sgi, "RNAinteract" ) )

  controls = sgi@targets$Symbol[sgi@targets$group %in% c("pos","neg")]

  zz <- 0
  ZZ <- sgi@S * sgi@C
  ## for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {

      q.value <- getData(sgi, screen = screen, channel = c, type="q.value", format="targetMatrix")
      p.value <- getData(sgi, screen = screen, channel = c, type="p.value", format="targetMatrix")
      q.value.limma <- getData(sgilimma, screen = screen, channel = c, type="q.value", format="targetMatrix")
      p.value.limma <- getData(sgilimma, screen = screen, channel = c, type="p.value", format="targetMatrix")
      gene1 <- matrix(row.names(p.value),nr=dim(p.value)[1], nc=dim(p.value)[2])
      gene2 <- matrix(colnames(p.value),nr=dim(p.value)[2], nc=dim(p.value)[1])
      gene2 <- t(gene2)
      maint <- getMain(sgi, design="template", screen=screen, channel=c, summary="target")
      mainq <- getMain(sgi, design="query", screen=screen, channel=c, summary="target")
      main = RNAinteract:::invtransform(sgi, 0.5 * (maint + mainq), channel=c)
      main1 <- matrix(main,nr=dim(p.value)[1], nc=dim(p.value)[2])
      main2 <- matrix(main,nr=dim(p.value)[2], nc=dim(p.value)[1])
      main2 <- t(main2)
      NI <- getData(sgi, screen = screen, channel = c, type="ni.model", format="targetMatrix", do.inv.trafo=TRUE)
      D <- getData(sgi, screen = screen, channel = c, type="data", format="targetMatrix", do.trafo=FALSE)
      PI <- getData(sgi, screen = screen, channel = c, type="pi", format="targetMatrix", do.inv.trafo=TRUE)

      UT <- upper.tri(p.value)
      m1 <- match(gene1[UT], sgi@targets$Symbol)
      m2 <- match(gene2[UT], sgi@targets$Symbol)
      GL <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value.ttest = q.value[UT], p.value.ttest = p.value[UT],
                       q.value.limma = q.value.limma[UT], p.value.limma = p.value.limma[UT],
                       main1 = main1[UT], main2 = main2[UT],
                       neg = getMainNeg(sgi, screen=screen, channel=c, do.inv.trafo=TRUE),
                       NI = NI[UT], Measured = D[UT], pi = PI[UT],
                       FBgn1 = sgi@targets$GID[m1], FBgn2 = sgi@targets$GID[m2],
                       CG1 = sgi@targets$AnnotationSymbol[m1], CG2 = sgi@targets$AnnotationSymbol[m2],
                       Name1 = sgi@targets$Name[m1], Name2 = sgi@targets$Name[m2])

      ## main.log = 0.5 * (maint + mainq)
      ## main1.log <- matrix(main.log,nr=dim(p.value)[1], nc=dim(p.value)[2])
      ## main2.log <- matrix(main.log,nr=dim(p.value)[2], nc=dim(p.value)[1])
      ## main2.log <- t(main2.log)
      ## NI.log <- getData(sgi, screen = s, channel = c, type="ni.model", format="targetMatrix", do.inv.trafo=FALSE)
      ## D.log <- getData(sgi, screen = s, channel = c, type="data", format="targetMatrix", do.trafo=TRUE)
      ## PI.log <- getData(sgi, screen = s, channel = c, type="pi", format="targetMatrix", do.inv.trafo=FALSE)

      ## GL.log <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value = q.value[UT], p.value = p.value[UT],
      ##                  main1 = main1.log[UT], main2 = main2.log[UT],
      ##                  neg = getMainNeg(sgi, screen=s, channel=c, do.inv.trafo=FALSE),
      ##                  NI = NI.log[UT], Measured = D.log[UT], pi = PI.log[UT])

      isctrl1 <- matrix(names(maint) %in% controls,nr=dim(p.value)[1], nc=dim(p.value)[2])
      isctrl2 <- matrix(names(maint) %in% controls,nr=dim(p.value)[2], nc=dim(p.value)[1])
      isctrl2 <- t(isctrl2)
      isctrl <- isctrl1[UT] | isctrl2[UT]

      pv <- GL$p.value.ttest
      pv[isctrl] <- pv[isctrl] + 2
      I <- order(pv)
      GL <- GL[I,]
      ## GL.log <- GL.log[I,]
      GL.text <- GL
      for (i in 3:10) {
        GL.text[[i]] <- sprintf("%0.3f",GL.text[[i]])
      }

      con <- file(sprintf("Tab3_%s.txt", c), "wb")
      t <- switch(c,
                  nrCells = "# Supplementary Table 3: Genetic interaction data based on number of cells",
                  area = "# Supplementary Table 3: Genetic interaction data based on area",
                  intensity = "# Supplementary Table 3: Genetic interaction data based on intensity") 
      writeLines(t,con)
      writeLines("# ",con)																
      writeLines("# gene1:\tGene Symbol of first gene",con)
      writeLines("# gene2:\tGene Symbol of second gene",con)
      writeLines("# q.value:\tq-value computed by the method of Storey and Tibshirani (2003) using p-values from column p.value",con)
      writeLines("# p.value:\tp-value computed by a t-statistic with n-1 degrees of freedom (t-test). Variance is computed as the maximum of local and pooled variance. These p-values are more conservative than p.value.limma. These p-values are used in Horn et al. 2011. ",con)
      writeLines("# q.value.limma:\tq-value computed by the method of Storey and Tibshirani (2003) using p-values from column p.value.limma  ",con)
      writeLines("# p.value.limma:\tp-value computed by a moderated t-test with the bioconductor-package limma.",con)
      writeLines("# main1:\tThe main effect (single knock-down effect) of the first gene as a multiplicative factor.",con)
      writeLines("# main2:\tThe main effect (single knock-down effect) of the second gene as a multiplicative factor.",con)
      writeLines("# neg:\tThe phenotypic measurement of the negative control after correction of systematic measurement errors (e.g. time, plate, edge effects). This is used as the baseline for 'no perturbation'.",con)
      writeLines("# NI:\tThe predicted value for the pairwise dsRNA treatment for non-interacting (NI) gene assuming a multiplicative model.",con)
      writeLines("# Measured:\tThe phenotypic measurement of the pairwise dsRNA after correction of systematic measurement errors (e.g. time, plate, edge effects).",con)
      writeLines("# pi:\tThe pairwise interaction score: ratio of Measured over NI.",con)
      writeLines("# FBgn1:\tFlybase identifier of the first gene",con)
      writeLines("# FBgn2:\tFlybase identifier of the second gene",con)
      writeLines("# CG1:\tCG-identifier of the first gene",con)
      writeLines("# CG2:\tCG-identifier of the second gene",con)																
      writeLines("# Name1:\tName of the first gene",con)														
      writeLines("# Name2:\tName of the second gene",con)																
      writeLines("# ",con)																
      close(con)

      write.table(GL, file = sprintf("Tab3_%s.txt",c), append = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

    }

  c = "nrCells"
  q.value <- getData(sgi3T2, screen = screen, channel = c, type="q.value", format="targetMatrix")
  p.value <- getData(sgi3T2, screen = screen, channel = c, type="p.value", format="targetMatrix")
  gene1 <- matrix(row.names(p.value),nr=dim(p.value)[1], nc=dim(p.value)[2])
  gene2 <- matrix(colnames(p.value),nr=dim(p.value)[2], nc=dim(p.value)[1])
  gene2 <- t(gene2)

  UT <- upper.tri(p.value)
  m1 <- match(gene1[UT], sgi@targets$Symbol)
  m2 <- match(gene2[UT], sgi@targets$Symbol)
  GL <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value.T2 = q.value[UT], p.value.T2 = p.value[UT],
                   FBgn1 = sgi@targets$GID[m1], FBgn2 = sgi@targets$GID[m2],
                   CG1 = sgi@targets$AnnotationSymbol[m1], CG2 = sgi@targets$AnnotationSymbol[m2],
                   Name1 = sgi@targets$Name[m1], Name2 = sgi@targets$Name[m2])
  
      ## main.log = 0.5 * (maint + mainq)
      ## main1.log <- matrix(main.log,nr=dim(p.value)[1], nc=dim(p.value)[2])
      ## main2.log <- matrix(main.log,nr=dim(p.value)[2], nc=dim(p.value)[1])
      ## main2.log <- t(main2.log)
      ## NI.log <- getData(sgi, screen = s, channel = c, type="ni.model", format="targetMatrix", do.inv.trafo=FALSE)
      ## D.log <- getData(sgi, screen = s, channel = c, type="data", format="targetMatrix", do.trafo=TRUE)
      ## PI.log <- getData(sgi, screen = s, channel = c, type="pi", format="targetMatrix", do.inv.trafo=FALSE)

      ## GL.log <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value = q.value[UT], p.value = p.value[UT],
      ##                  main1 = main1.log[UT], main2 = main2.log[UT],
      ##                  neg = getMainNeg(sgi, screen=s, channel=c, do.inv.trafo=FALSE),
      ##                  NI = NI.log[UT], Measured = D.log[UT], pi = PI.log[UT])

  isctrl1 <- matrix(names(maint) %in% controls,nr=dim(p.value)[1], nc=dim(p.value)[2])
  isctrl2 <- matrix(names(maint) %in% controls,nr=dim(p.value)[2], nc=dim(p.value)[1])
  isctrl2 <- t(isctrl2)
  isctrl <- isctrl1[UT] | isctrl2[UT]

  pv <- GL$p.value.T2
  pv[isctrl] <- pv[isctrl] + 2
  I <- order(pv)
  GL <- GL[I,]
  ## GL.log <- GL.log[I,]
  GL.text <- GL
  for (i in 3:10) {
    GL.text[[i]] <- sprintf("%0.3f",GL.text[[i]])
  }

  con <- file("Tab3_HotellingT2.txt", "wb")
  t <- switch(c,
              nrCells = "Supplementary Table 3: Genetic interaction data based on number of cells",
              area = "Supplementary Table 3: Genetic interaction data based on area",
              intensity = "Supplementary Table 3: Genetic interaction data based on intensity") 
  writeLines(t,con)
  writeLines("",con)																
  writeLines("gene1:\tGene Symbol of first gene",con)
  writeLines("gene2:\tGene Symbol of second gene",con)
  writeLines("q.value:\tq-value computed by the method of Storey and Tibshirani (2003) using p-values from column p.value",con)
  writeLines("p.value:\tp-value computed by a Hotelling T2 test. ",con)
  writeLines("FBgn1:\tFlybase identifier of the first gene",con)
  writeLines("FBgn2:\tFlybase identifier of the second gene",con)
  writeLines("CG1:\tCG-identifier of the first gene",con)
  writeLines("CG2:\tCG-identifier of the second gene",con)																
  writeLines("Name1:\tName of the first gene",con)														
  writeLines("Name2:\tName of the second gene",con)																
  writeLines("",con)																
  close(con)

  write.table(GL, file = "Tab3_HotellingT2.txt", append = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
  
  invisible(NULL)
}



