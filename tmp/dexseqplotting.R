plotDEXSeq2 <- function (object, geneID, FDR = 0.1, fitExpToVar = "condition", 
          norCounts = FALSE, expression = TRUE, splicing = FALSE, displayTranscripts = FALSE, 
          names = FALSE, legend = FALSE, color = NULL, color.samples = NULL, 
          transcriptDb = NULL, additionalAnnotation = NULL, maxRowsMF = 2400, 
          ...) 
{
  stopifnot(is(object, "DEXSeqResults") | is(object, "DEXSeqDataSet"))
  if (!fitExpToVar %in% colnames(object@modelFrameBM)) {
    stop(sprintf("The value of the parameter fitExpToVar,'%s', is not a column name of the 'colData' DataFrame from the DEXSeqDataSet object.", 
                 fitExpToVar))
  }
  op <- sum(c(expression, splicing, norCounts))
  if (op == 0) {
    stop("Please indicate what would you like to plot\n")
  }
  if (!is.null(transcriptDb)) {
    stopifnot(is(transcriptDb, "TxDb"))
  }
  if (!is.null(additionalAnnotation)) {
    stopifnot(is(additionalAnnotation, "GRangesList"))
  }
  if (is(object, "DEXSeqResults")) {
    sampleData <- object@sampleData
    genomicData <- object$genomicData
    rt <- which(object$groupID == geneID)
    count <- t(t(object$countData[rt, ])/sampleData$sizeFactor)
    transcripts <- object$transcripts[rt]
    each <- object$padj[rt]
  }else {
    sampleData <- sampleAnnotation(object)
    genomicData <- rowRanges(object)
    rt <- which(mcols(object)$groupID == geneID)
    transcripts <- genomicData$transcripts[rt]
    mcols(genomicData) <- NULL
    count <- featureCounts(object, normalized = TRUE)[rt, 
    ]
    each <- rep(1, length.out = length(rt))
  }
  if (sum(count) == 0) {
    warning("No read counts falling in this gene, there is nothing to plot.")
    return()
  }
  if (FDR > 1 | FDR < 0) {
    stop("FDR has to be a numeric value between 0 - 1")
  }
  rango <- seq(along = rt)
  intervals <- (0:nrow(count))/nrow(count)
  numcond <- length(unique(sampleData[[fitExpToVar]]))
  numexons <- nrow(count)
  exoncol <- ifelse(each <= FDR, "#F219ED", "#CCCCCC")
  exoncol[is.na(exoncol)] <- "white"
  colorlines <- ifelse(each <= FDR, "#F219ED60", "#B3B3B360")
  colorlines[is.na(colorlines)] <- "#B3B3B360"
  colorlinesB <- ifelse(each <= FDR, "#9E109B", "#666666")
  colorlinesB[is.na(colorlinesB)] <- "#666666"
  if (length(unlist(start(genomicData))) > 0) {
    sub <- data.frame(start = start(genomicData[rt]), end = end(genomicData[rt]), 
                      chr = as.character(seqnames(genomicData[rt])), strand = as.character(strand(genomicData[rt])))
    rownames(sub) <- rownames(object)[rt]
    if (!is.null(additionalAnnotation)) {
      additionalHits <- findOverlaps(additionalAnnotation, 
                                     range(genomicData[rt]))
      additionalAnnotation <- additionalAnnotation[queryHits(additionalHits)]
      if (length(additionalAnnotation) == 0) {
        additionalAnnotation <- NULL
      }
    }
    rel <- (data.frame(sub$start, sub$end)) - min(sub$start)
    rel <- rel/max(rel[, 2])
    trans <- unique(unlist(transcripts))
    trans <- trans[!is.na(trans)]
    numberOfTrans <- length(trans) + length(additionalAnnotation)
    if ((displayTranscripts & !is.null(unlist(transcripts))) | 
        !is.null(additionalAnnotation)) {
      if (numberOfTrans > 40) {
        warning("This gene contains more than 40 transcripts annotated, only the first 40 will be plotted\n")
      }
      if (!displayTranscripts) {
        numberOfTrans <- numberOfTrans - length(trans)
      }
      mat <- seq_len(3 + min(numberOfTrans, 40))
      hei <- c(8, 1, 1.5, rep(1.5, min(numberOfTrans, 40)))
    }
    else {
      mat <- 1:3
      hei <- c(5, 1, 1.5)
    }
    if (op > 1) {
      hei <- c(rep(hei[1], op - 1), hei)
      mat <- c(mat, length(mat) + seq(along = op))
    }
    hei <- c(hei, 0.2)
    mat <- c(mat, length(mat) + 1)
    layout(matrix(mat), heights = hei)
    par(mar = c(2, 4, 4, 2))
  }else if (op > 1) {
    par(mfrow = c(op, 1))
  }
  if (is.null(color)) {
    if (numcond < 10) {
      color <- suppressWarnings(brewer.pal(numcond, "Set1")[seq_len(numcond)])
    } else {
      color <- rgb(colorRamp(brewer.pal(5, "Set1"))(seq(0, 
                                                        1, length.out = numcond)), maxColorValue = 255, 
                   alpha = 175)
    }
  }
  names(color) <- sort(unique(as.character((sampleData[[fitExpToVar]]))))
  if (expression | splicing) {
    if (!is(object, "DEXSeqResults")) {
      stop("To visualize beta estimates, please provide a DEXSeqResults object")
    }
    effects <- DEXSeq:::getEffectsForGene(geneID, object, maxRowsMF, 
                                 fitExpToVar)
    if (is.null(effects[["splicing"]])) {
      return()
    }
    if (!all(rownames(sub) %in% rownames(effects[["splicing"]]))) {
      return()
    }
  }
  if (expression) {
    coeff <- effects[["expression"]]
    # coeff <- coeff[rownames(sub), ]
    coeff <- exp(coeff)
    ylimn <- c(0, max(coeff, na.rm = TRUE))
    coeff <- DEXSeq:::vst(coeff, object)
    DEXSeq:::drawPlot(matr = coeff, ylimn, object, intervals, rango, 
             textAxis = "Expression", rt = rt, color = rep(color[colnames(coeff)], 
                                                           each = numexons), colorlines = colorlines, ...)
  }
  if (splicing) {
    coeff <- effects[["splicing"]]
    # coeff <- coeff[rownames(sub), ]
    coeff <- exp(coeff)
    ylimn <- c(0, max(coeff, na.rm = TRUE))
    coeff <- DEXSeq:::vst(coeff, object)
    DEXSeq:::drawPlot(matr = coeff, ylimn, object, intervals, rango, 
             textAxis = "Exon usage", rt = rt, color = rep(color[colnames(coeff)], 
                                                           each = numexons), colorlines = colorlines)
  }
  if (norCounts) {
    ylimn <- c(0, max(count, na.rm = TRUE))
    count <- DEXSeq:::vst(count, object)
    if (is.null(color.samples)) {
      colorcounts <- rep(color[as.character(sampleData[[fitExpToVar]])], 
                         each = numexons)
    }
    else {
      colorcounts <- rep(color.samples, each = numexons)
    }
    DEXSeq:::drawPlot(matr = count, ylimn, object, intervals, rango, 
             textAxis = "Normalized counts", rt = rt, color = colorcounts, 
             colorlines = colorlines, ...)
  }
  if (length(unlist(start(genomicData))) > 0) {
    par(mar = c(0, 4, 0, 2))
    plot.new()
    segments(apply((rbind(rel[rango, 2], rel[rango, 1])), 
                   2, median), 0, apply(rbind(intervals[rango], intervals[rango + 
                                                                            1] - ((intervals[rango + 1] - intervals[rango]) * 
                                                                                    0.2)), 2, median), 1, col = colorlinesB)
    par(mar = c(1.5, 4, 0, 2))
    drawGene(min(sub$start), max(sub$end), tr = sub, exoncol = exoncol, 
             names, trName = "Gene model", cex = 0.8)
    if (length(unlist(transcripts)) > 0) {
      i <- 1
      if (displayTranscripts) {
        for (i in seq_len(min(length(trans), 40))) {
          logicexons <- sapply(transcripts, function(x) {
            length(which(x == trans[i]))
          })
          tr <- reduce(IRanges(sub$start[logicexons == 
                                           1], sub$end[logicexons == 1]))
          if (is.null(transcriptDb)) {
            tr <- as.data.frame(tr)[, c("start", "end")]
            drawGene(min(sub$start), max(sub$end), tr = tr, 
                     exoncol = "black", names, trName = trans[i], 
                     cex = 0.8)
          }
          else {
            codingRanges <- select(transcriptDb, keys = trans[i], 
                                   columns = c("CDSSTART", "CDSEND"), keytype = "TXNAME")
            if (is.na(any(codingRanges$CDSSTART))) {
              tr <- as.data.frame(tr)[, c("start", "end")]
              drawGene(min(sub$start), max(sub$end), 
                       tr = tr, exoncol = NULL, names, trName = trans[i], 
                       cex = 0.8, miny = 0.25, maxy = 0.75)
            }
            else {
              codingRanges <- IRanges(codingRanges$CDSSTART, 
                                      codingRanges$CDSEND)
              utrRanges <- setdiff(tr, codingRanges)
              drawGene(min(sub$start), max(sub$end), 
                       tr = as.data.frame(codingRanges)[, c("start", 
                                                            "end")], exoncol = "black", names, 
                       trName = trans[i], cex = 0.8, drawNames = FALSE, 
                       drawIntronLines = FALSE)
              if (length(utrRanges) > 0) {
                drawGene(min(sub$start), max(sub$end), 
                         tr = as.data.frame(utrRanges)[, c("start", 
                                                           "end")], exoncol = NULL, names, trName = trans[i], 
                         cex = 0.8, drawNames = FALSE, drawIntronLines = FALSE, 
                         newPanel = FALSE, miny = 0.25, maxy = 0.75)
              }
              drawGene(min(sub$start), max(sub$end), 
                       tr = as.data.frame(tr)[, c("start", "end")], 
                       exoncol = "black", names, trName = trans[i], 
                       cex = 0.8, newPanel = FALSE, drawExons = FALSE)
            }
          }
        }
      }
      if (!is.null(additionalAnnotation)) {
        for (j in seq_along(additionalAnnotation)) {
          tr <- as.data.frame(additionalAnnotation[[j]])[, 
                                                         c("start", "end")]
          drawGene(min(sub$start), max(sub$end), tr = tr, 
                   exoncol = "darkred", names, trName = names(additionalAnnotation)[j], 
                   cex = 0.8, introncol = "darkred")
          i <- i + 1
          if (i > 40) 
            break
        }
      }
    }
    axis(1, at = round(seq(min(sub$start), max(sub$end), 
                           length.out = 10)), labels = round(seq(min(sub$start), 
                                                                 max(sub$end), length.out = 10)), pos = 0, lwd.ticks = 0.2, 
         padj = -0.7, ...)
  }
  if (legend) {
    mtext(paste(geneID), side = 3, adj = 0.25, 
          padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
    posforlegend <- seq(0.7, 0.9, length.out = numcond)
    for (i in seq(along = color)) {
      mtext(names(color[i]), side = 3, adj = posforlegend[i], 
            padj = 1.5, line = 0, outer = TRUE, col = color[i], 
            ...)
    }
  } else {
    mtext(paste(geneID), side = 3, adj = 0.5, 
          padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
  }
}
