###########################################
## MULTI-seq sample classification suite ##
###########################################

###################
## classifyCells ##
###################
classifyCells <- function(barTable, q) {
  require(KernSmooth)

  ## Normalize Data: Log2 Transform, mean-center
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }

  ## Pre-allocate memory to save processing time/memory
  n_BC <- ncol(barTable.n) # Number of barcodes
  n_cells <- nrow(barTable.n) # Numer of cells
  bc_calls <- rep("Negative",n_cells) # Barcode classification for each cell
  names(bc_calls) <- rownames(barTable.n)

  ## Generate thresholds for each barcode:
  ## 1. Gaussian KDE with bad barcode detection, outlier trimming
  ## 2. Define local maxima for GKDE
  ## 3. Split maxima into low and high subsets, adjust low if necessary
  ## 4. Threshold and classify cells according to user-specified inter-maxima quantile

  for (i in 1:n_BC) {
    ## Step 1: GKDE
    model <- tryCatch( { approxfun(bkde(barTable.n[,i], kernel="normal")) },
                       error=function(e) { print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) } )
    if (class(model) == "character") { next }
    x <-  seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999), length.out=100)

    ## Step 2: Local maxima definition
    extrema <- localMaxima(model(x))

    if (length(extrema) <= 1) {
      print(paste0("No threshold found for ", colnames(barTable.n)[i],"..."))
      next
    }

    ## Step 3: Select maxima
    ## Assumes negative cells are largest mode
    ## Assumes positive cells are highest extreme -- favors negatives over doublets
    low.extreme <- extrema[which.max(model(x)[extrema])]
    high.extreme <- max(extrema)
    if (low.extreme == high.extreme) {print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) }

    ## Step 4: Threshold and classify cells
    thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
    cell_i <- which(barTable.n[,i] >= thresh)
    n <- length(cell_i)
    if (n == 0) { next } # Skip to next barcode if no cells classified
    bc_calls[cell_i] <- sapply(bc_calls[cell_i],
                               FUN = function(x) {
                                 if (x == "Negative") {
                                   return(colnames(barTable.n)[i])
                                 } else {
                                   return("Doublet")
                                 } })
  }

  return(bc_calls)

}

################
## findThresh ##
################
findThresh <- function(call.list) {
  require(reshape2)

  res <- as.data.frame(matrix(0L, nrow=length(call.list), ncol=4))
  colnames(res) <- c("q","pDoublet","pNegative","pSinglet")
  q.range <- unlist(strsplit(names(call.list), split="q="))
  res$q <- as.numeric(q.range[grep("0", q.range)])
  nCell <- length(call.list[[1]])

  for (i in 1:nrow(res)) {
    temp <- table(call.list[[i]])
    if ( "Doublet" %in% names(temp) == TRUE ) { res$pDoublet[i] <- temp[which(names(temp) == "Doublet")] }
    if ( "Negative" %in% names(temp) == TRUE ) { res$pNegative[i] <- temp[which(names(temp) == "Negative")] }
    res$pSinglet[i] <- sum(temp[which(names(temp) != "Doublet" & names(temp) != "Negative")])
  }

  res <- melt(res, id.vars="q")
  res[,4] <- res$value/nCell
  colnames(res)[2:4] <- c("Subset","nCells","Proportion")

  extrema <- res$q[localMaxima(res$Proportion[which(res$Subset == "pSinglet")])]
  return(list("extrema" = extrema, "res" = res)
}

#################
## localMaxima ##
#################
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

######################
## findReclassCells ##
######################
findReclassCells <- function(barTable, neg.cells) {
  require(KernSmooth)

  ## Normalize Data: Log2 Transform, mean-center
  print("Normalizing barode data...")
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }
  barTable.n_neg.cells <- barTable.n[neg.cells, ]

  ## Pre-allocate memory to save processing time/memory
  print("Pre-allocating data structures...")
  n_BC <- ncol(barTable.n)
  bcs <- colnames(barTable.n)
  n_neg <- length(neg.cells)
  bc.mats <- list()
  q.range <- seq(0.01, 0.99, by=0.02)
  ref.mat <- as.data.frame(matrix(0L, nrow=n_neg, ncol=length(q.range)))
  rownames(ref.mat) <- neg.cells
  colnames(ref.mat) <- paste("q", q.range, sep="_")
  for (i in 1:n_BC) {
    bc.mats[[i]] <- ref.mat
  }
  names(bc.mats) <- bcs

  ## Iterate through BCs, modeling distribution using GKDE
  ## Define low and high maxima before testing all negative cells across all q
  print("Determining classifications for negative cells across all BCs, all q...")
  for (i in 1:n_BC) {
    ## Step 1: GKDE
    model <- approxfun(bkde(barTable.n[,i], kernel="normal"))
    x <- seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999), length.out=100)

    ## Step 2: Local maxima definition
    extrema <- localMaxima(model(x))

    ## Step 3: Select maxima, avoiding noisy extremes (if possible)
    ## Assumes negative cells are largest mode
    ## Assumes positive cells are highest extreme -- favors negatives over doublets
    low.extreme <- extrema[which.max(model(x)[extrema])]
    high.extreme <- max(extrema)

    ## Iterate through q, tracking classiication status for negative cells
    col.ind <- 0
    for (q in q.range) {
      col.ind <- col.ind + 1
      thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
      cell_i <- which(barTable.n_neg.cells[,i] >= thresh)
      bc.mats[[i]][cell_i, col.ind] <- 1
    }
  }

  ## Determine classification status across q.range for negative cells across all barcodes
  print("Computing classification stability...")
  bc.mat.final <- Reduce('+',bc.mats)

  ## Compute classification stability ~ i.e., number of q for which cell has 1 classification
  res <- as.data.frame(matrix(0L, nrow=n_neg, ncol=2))
  colnames(res) <- c("ClassStability","Reclassification")
  rownames(res) <- neg.cells
  res$Reclassification <- "Negative"
  for (cell in neg.cells) {
    res[cell,"ClassStability"] <- length(which(bc.mat.final[cell, ] == 1))
  }

  ## Extract classifications for cells above cluster.stability reclassification threshold
  # Define reclassifiable cells
  print(paste("Extracting rescued classifications...",sep=""))
  cells2reclass <- neg.cells[which(res$ClassStability >= 1)]

  # Find q at which reclassifiable cells have single classification, extract corresponding bc
  for (cell in cells2reclass) {
    q.ind <- max(which(bc.mat.final[cell, ] == 1))

    for (bar in 1:n_BC) {
      temp <- bc.mats[[bar]][cell, q.ind]
      if (temp == 1) { res[cell, "Reclassification"] <- bcs[bar] }
    }
  }

  # Return only cells with an available classification
  res <- res[-which(res$Reclassification == "Negative"), ]
  return(res)

}

#################
## rescueCells ##
#################
rescueCells <- function(barTable, classifications, reclassifications) {

  ## Extract subset of data with equal numbers of high-confidence classifications
  smallest.bc <- min(table(classifications)) # Find least-frequent BC classification
  bcs <- colnames(barTable)
  conf.cells <- NULL
  for (bc in bcs) { # Get n cellIDs for each bc, where n is number of cells in least-frequent BC classification
    temp <- names(classifications)[which(classifications == bc)]
    temp <- sample(temp, smallest.bc, replace=FALSE)
    conf.cells <- c(conf.cells, temp)
  }
  conf.dat.classifications <- classifications[conf.cells]

  ## Normalize data including high-confidence and negative cells prior to k-means
  kmeans.dat <- barTable[c(conf.cells, rownames(reclassifications)), ]
  kmeans.dat <- as.data.frame(log2(kmeans.dat))
  for (i in 1:ncol(kmeans.dat)) {
    ind <- which(is.finite(kmeans.dat[,i]) == FALSE)
    kmeans.dat[ind,i] <- 0
    kmeans.dat[,i] <- kmeans.dat[,i]-mean(kmeans.dat[,i])
  }

  ## Initialize k-means results output
  conf.res.temp <- as.data.frame(matrix(0L, nrow=length(conf.cells), ncol=3))
  colnames(conf.res.temp) <- c("GT","GT.adj","KMEANS")
  rownames(conf.res.temp) <- conf.cells
  conf.res.temp$GT <- factor(conf.dat.classifications, levels=colnames(barTable))
  conf.res.temp$GT.adj <- as.numeric(conf.res.temp$GT)

  neg.res.list <- list()
  for (iter in 1:smallest.bc) {
    neg.res.temp <- reclassifications
    neg.res.temp$Reclassification <- factor(neg.res.temp$Reclassification, levels=colnames(barTable))
    neg.res.temp[,"Reclass.adj"] <- as.numeric(neg.res.temp$Reclassification)
    neg.res.temp[, "KMEANS"] <- rep(0L)
    neg.res.list[[iter]] <- neg.res.temp
  }

  ## Initial match rate results output
  conf.match.rate <- rep(0L, smallest.bc)
  neg.match.rate <- as.data.frame(matrix(0L, nrow=max(neg.res.temp$ClassStability), ncol=3))
  colnames(neg.match.rate) <- c("ClassStability", "MatchRate_mean","MatchRate_sd")
  neg.match.rate$ClassStability <- 1:nrow(neg.match.rate)

  ## Perform k-means iterations and match rate computations
  for (iter in 1:smallest.bc) {
    ## Set centers
    print(iter)
    n_BC <- ncol(barTable)
    centers.ind <- seq(1, smallest.bc*n_BC, by=smallest.bc) + (iter-1)
    centers <- kmeans.dat[centers.ind, ]

    ## Run k-means, store results for confident and negative cells
    kmeans_temp <- kmeans(kmeans.dat, centers = centers)
    conf.res.temp$KMEANS <- kmeans_temp$cluster[1:length(conf.cells)]
    neg.res.list[[iter]]$KMEANS <- kmeans_temp$cluster[(length(conf.cells)+1):length(kmeans_temp$cluster)]

    ## Compute match rate for confident cells
    conf.match.rate[iter] <- 1-(length(which(conf.res.temp$GT.adj != conf.res.temp$KMEANS))/nrow(conf.res.temp))
  }

  ## Compute match rate statistics across classification stabilities
  for (cs in 1:nrow(neg.match.rate)) {
    ind <- which(neg.res.list[[1]]$ClassStability == cs)
    match.rate.means <- NULL

    for (iter in 1:smallest.bc) {
      neg.res.temp <- neg.res.list[[iter]][ind, ]
      x <- 1-(length(which(neg.res.temp$Reclass.adj != neg.res.temp$KMEANS))/nrow(neg.res.temp))
      match.rate.means <- c(match.rate.means,x)
    }

    neg.match.rate$MatchRate_mean[cs] <- mean(match.rate.means)
    neg.match.rate$MatchRate_sd[cs] <- sd(match.rate.means)

  }

  res <- rbind(c(100,mean(conf.match.rate),sd(conf.match.rate)),neg.match.rate)
  return(res)

}

#############
## barTSNE ##
#############
barTSNE <- function(barTable) {
  require(Rtsne)

  ## Normalize barcode count matrix
  n_BC <- ncol(barTable)
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:n_BC) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }

  ## Run tSNE
  tsne.res <- Rtsne(barTable.n, dims=2, initial_dims=n_BC, verbose=TRUE, check_duplicates=FALSE, max_iter=2500)
  tsne.embedding <- as.data.frame(tsne.res$Y)
  colnames(tsne.embedding) <- c("TSNE1","TSNE2")
  tsne.embedding[,3:(n_BC+2)] <- barTable.n
  rownames(tsne.embedding) <- rownames(barTable)

  return(tsne.embedding)
}

###########
## findQ ##
###########
findQ <- function(res, extrema) {
  temp <- res[which(res$q %in% extrema & res$Subset == "pSinglet"), ]
  x <- as.numeric(temp$q[which.max(temp$Proportion)])
  return(x)
}
