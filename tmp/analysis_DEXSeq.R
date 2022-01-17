# DEXSeq analysis ----

library(DEXSeq)



# Creating sample table
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3", 
                 "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",  
                "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end", 
               "single-end", "single-end", "paired-end", "paired-end" ) )

runDEXseq <- function(dds){
  feat_group <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  svs <- as.character(design(dds))[2] %>% 
    stringr::str_split("\\+ ", n = 2, simplify = TRUE) %>% 
    .[2]
  des <- as.formula(
    paste0("~ sample + exon +", groupCol, ":exon + ", svs)
  )
  sampleData <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% #row.names = .$id) %>%
    dplyr::select(all_of(groupCol), starts_with("sv"))
  #   rename(sample = id)
  # sampleData <- sampleData%>% 
  #   bind_rows(sampleData) %>% 
  #   mutate(sample = factor(sample),
  #          exon = factor(c(rep("this", nrow(sampleData)),
  #                          rep("others", nrow(sampleData))))) %>% 
  #   as.data.frame(row.names = as.character(.$sample))
  
    
  
  # Convert to DEXSeq object
  dxd <- DEXSeqDataSet(
    countData = assay(dds),
    sampleData = sampleData,
    design = des,
    featureID = feat_group[, 1],
    groupID = feat_group[, 2]
    )
  
}

res_dexseq <- DEXSeq(dxd)

# Normalization
dxd <- DEXSeq::estimateSizeFactors( dxd )

# Estimate dispersion
dxd <- DEXSeq::estimateDispersions( dxd )

res_dexseq <- DEXSeq::DEXSeqResults( dxd )

# Shrinkage diagostic
plotDispEsts( dxd )




DEXSeqDataSet <- function( countData, sampleData, design= ~ sample + exon + condition:exon , featureID, groupID, featureRanges=NULL, transcripts=NULL, alternativeCountData=NULL)
{
  ### Checking inputs ###
  if( !(is( countData, "matrix" ) | is( countData, "data.frame" )) )
    stop( "Unexpected input: the parameter 'countData' must be either a matrix or a data.frame", 
          call.=FALSE )
  countData <- as.matrix( countData )
  if( !( is( featureID, "character" ) | is( featureID, "factor" ) ) )
    stop( "Unexpected input: the parameter 'featureID' must be either a character or a factor", 
          call.=FALSE )
  if( !( is( groupID, "character" ) | is( groupID, "factor" ) ) )
    stop( "Unexpected input: the parameter 'groupID' must be either a character or a factor", 
          call.=FALSE )
  if( !is( sampleData, "data.frame" ) )
    stop( "Unexpected input: the parameter 'sampleData' must be a data.frame", 
          call.=FALSE )
  rowNumbers <- nrow( countData )
  if( length(groupID) != rowNumbers )
    stop( "Unexpected length of 'groupID' parameter, it must be the same as the number of rows of countData", 
          call.=FALSE )
  if( length(featureID) != rowNumbers )
    stop( "Unexpected length of 'featureID' parameter, it must be the same as the number of rows of countData", 
          call.=FALSE )
  if( nrow(sampleData) != ncol(countData) )
    stop( "Unexpected number of rows of the 'sampleData' parameter, it must be the same as the number of columns of countData", 
          call.=FALSE )
  
  modelFrame <- cbind(
    sample = as.factor(rownames(sampleData)), sampleData )
  modelFrame <- rbind( cbind(modelFrame, exon = "this"),
                       cbind(modelFrame, exon = "others"))
  modelFrame$exon <- as.factor(modelFrame$exon)
  rownames(modelFrame) <- NULL
  colData <- DataFrame( modelFrame )
  
  if( !"exon" %in% all.vars( design ) )
    stop("The design formula does not specify an interaction contrast with the variable 'exon'", 
         call.=FALSE )
  
  allVars <- all.vars(design)
  if( any(!allVars %in% colnames( colData )) ){
    notPresent <- allVars[!allVars %in% colnames( colData ) ]
    notPresent <- paste(notPresent, collapse=",")
    stop(sprintf("The variables '%s', present in the design formula must be columns of 'sampleData'", 
                 notPresent ),
         call.=FALSE )
  }
  
  if( any( grepl(" |:", groupID ) | grepl(" |:", featureID) ) ) {
    warning("empty spaces or ':' characters were found either in your groupIDs or in your featureIDs, these will be removed from the identifiers")
    groupID <- gsub(" |:", "", groupID)
    featureID <- gsub(" |:", "", featureID)
  }
  
  rownames( countData ) <- paste( groupID, featureID, sep=":" )
  forCycle <- split( seq_len(nrow( countData )), as.character( groupID ) )
  
  if( is.null(alternativeCountData) ){
    others <- lapply( forCycle, function(i){
      sct <- countData[i, , drop = FALSE]
      rs <- t( vapply( seq_len(nrow(sct)), function(r) colSums(sct[-r, , drop = FALSE]), numeric(ncol(countData) ) ) )
      rownames(rs) <- rownames(sct)
      rs })
    others <- do.call(rbind, others)
  }else{
    stopifnot( identical(dim(countData), dim(alternativeCountData)) )
    stopifnot( identical( colnames(countData), colnames(alternativeCountData)))
    others <- alternativeCountData
    rownames( others ) <- paste( groupID, featureID, sep=":" )
    
  }
  
  stopifnot( all( rownames(countData) %in% rownames(others) ) )
  others <- others[rownames(countData),]
  nCountData <- cbind( countData, others )
  colnames(nCountData) <- NULL
  
  if( !is.null(featureRanges) ){
    stopifnot(is(featureRanges, "GRanges") ||
                is(featureRanges, "GRangesList"))
    se <- SummarizedExperiment( nCountData, colData=colData, rowRanges=featureRanges )
  }else{
    se <- SummarizedExperiment( nCountData, colData=colData )
  }
  
  names(assays(se))[1] = "counts"
  mcols( se )$featureID <- featureID
  mcols( se )$groupID <- groupID
  mcols( se )$exonBaseMean <- rowMeans( countData )
  mcols( se )$exonBaseVar <- matrixStats::rowVars( countData )
  
  if( !is.null(transcripts) ){
    mcols(se)$transcripts <- transcripts
  }
  
  rownames(se) <- paste( groupID, featureID, sep=":")
  rse <- as( se, "RangedSummarizedExperiment" )
  mcols(rse) <- mcols(se)
  
  dds <- DESeqDataSet( rse, design, ignoreRank=TRUE )
  
  modelFrame <- makeBigModelFrame(dds)
  
  dxd <- new( "DEXSeqDataSet", dds, modelFrameBM=modelFrame )
  return(dxd)
}

makeBigModelFrame <- function(object){
  groupID <- mcols(object)$groupID
  featureID <- mcols(object)$featureID
  sampleData <- as.data.frame(colData(object)[colData(object)$exon == "this",])
  numExonsPerGene <- table(groupID)
  maxGene <- names(which.max(numExonsPerGene))
  rows <- mcols(object)$groupID %in%  maxGene
  numExons <- sum( rows )
  exonCol <-
    rep(factor(featureID[rows]), nrow(sampleData))
  modelFrame <- data.frame(
    sample=rep( sampleData$sample, each=numExons),
    exon = exonCol )
  varNames <- colnames(sampleData)[!colnames(sampleData) %in% c("sample", "exon")]
  for( i in varNames ){
    modelFrame[[i]] <- rep( sampleData[[i]], each=numExons )
  }
  modelFrame$dispersion <- NA
  if( is.null(object$sizeFactor )){
    modelFrame$sizeFactor <- NA
  }
  modelFrame$count <- NA
  modelFrame
}

DEXSeq <- function( object,
                    fullModel=design(object),
                    reducedModel = ~ sample + exon,
                    BPPARAM=BiocParallel::MulticoreParam(workers=4),
                    fitExpToVar="group_nr", quiet=TRUE ){
  stopifnot(is( object, "DEXSeqDataSet") )
  # Temporary hack for backward compatibility with "old" DEXSeqDataSet
  # objects. Remove once all serialized DEXSeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  message("Estimating size factors")
  object <- estimateSizeFactors( object )
  message("Estimating dispersions")
  object <- estimateDispersions( object, formula=fullModel, BPPARAM=BPPARAM, quiet=FALSE)
  message("Testing for DEU")
  object <- testForDEU( object, reducedModel=reducedModel, fullModel=fullModel, BPPARAM=BPPARAM )
  message("Estimating Fold Changes")
  object <- estimateExonFoldChanges( object, fitExpToVar=fitExpToVar )
  message("Extracting results")
  res <- DEXSeqResults( object )
  return(res)
}



estimateDispersions.DEXSeqDataSet <-
  function( object, fitType=c("parametric","local","mean"),
            maxit=100, niter=10, quiet=FALSE, formula=design(object), BPPARAM=SerialParam())
  {
    # Temporary hack for backward compatibility with "old" DEXSeqDataSet
    # objects. Remove once all serialized DEXSeqDataSet objects around have
    # been updated.
    if (!.hasSlot(object, "rowRanges"))
      object <- updateObject(object)
    if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
      stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
    }
    if (!is.null(dispersions(object))) {
      if (!quiet) message("you had estimated dispersions, replacing these")
      mcols(object) <- mcols(object)[,!(mcols(mcols(object))$type %in% c("intermediate","results"))]
    }
    stopifnot(length(maxit)==1)
    fitType <- match.arg(fitType, choices=c("parametric","local","mean"))
    
    allVars <- all.vars(formula)
    if( any(!allVars %in% colnames( colData(object) )) ){
      notPresent <- allVars[!allVars %in% colnames( colData(object) )]
      notPresent <- paste(notPresent, collapse=",")
      stop(sprintf("the variables '%s' of the parameter 'formula' are not specified in the columns of the colData", notPresent ) )
    }
    
    if( is( BPPARAM, "SerialParam" ) ){
      numParts <- 1L
    }else{
      numParts <- BPPARAM$workers
    }
    
    splitParts <- sort( rep( seq_len( numParts ), length.out=nrow(object) ) )
    splitObject <- split( object, splitParts )
    
    modelMatrix <- rmDepCols(
      model.matrix(formula, as.data.frame(colData(object))))
    
    splitObject <- bplapply( splitObject,
                             function(x, ... ){
                               #          library(DEXSeq)
                               estimateDispersionsGeneEst(x,
                                                          maxit=maxit, quiet=quiet,
                                                          modelMatrix = modelMatrix,
                                                          niter = niter)},
                             maxit=maxit, quiet=quiet,
                             modelMatrix=modelMatrix,
                             niter=niter,
                             BPPARAM=BPPARAM )
    
    mergeObject <- do.call( rbind, splitObject )
    matchedNames <- match( rownames(object), rownames(mergeObject))
    mcols(object) <- mcols( mergeObject )[matchedNames,]
    assays(object) <- assays(mergeObject[matchedNames,])
    
    mcols(object)$baseMean <- rowMeans( featureCounts(object, normalized=TRUE) )
    mcols(object)$baseVar <- mcols(object)$exonBaseVar
    mcols(object)$allZero <- unname( rowSums( featureCounts(object)) == 0 |
                                       rowSums(counts(object, normalized = TRUE)[, colData(object)$exon == "others"]) ==0 )
    
    object <- estimateDispersionsFit(object, fitType=fitType, quiet=quiet)
    
    dispPriorVar <- estimateDispersionsPriorVar(object, modelMatrix=modelMatrix)
    
    splitObject <- split( object, splitParts )
    
    splitObject <- bplapply( splitObject,
                             function(x, ... ){
                               #          library(DEXSeq)
                               estimateDispersionsMAP(x,
                                                      maxit=maxit,
                                                      quiet=quiet,
                                                      modelMatrix=modelMatrix,
                                                      dispPriorVar=dispPriorVar)
                             },
                             maxit=maxit, quiet=quiet,
                             modelMatrix=modelMatrix,
                             dispPriorVar=dispPriorVar,
                             BPPARAM=BPPARAM )
    
    mergeObject <- do.call( rbind, splitObject )
    matchedNames <- match( rownames(object), rownames(mergeObject) )
    mcols(object) <- mcols( mergeObject )[matchedNames,]
    mcols(object)$baseMean <- unname( rowMeans( counts(object, normalized=TRUE) ) )
    mcols(object)$baseVar <- unname( rowVars( counts(object, normalized=TRUE) ) )
    mcols(object)$dispersion <- pmin( mcols(object)$dispersion, ncol(object) )
    object
  }