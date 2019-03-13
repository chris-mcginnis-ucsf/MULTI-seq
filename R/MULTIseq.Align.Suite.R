
###############################################################################################################################
## MULTIseq.align aligns barcode reads to a reference in 10000-cell buckets, returning a cellIDs x ref barcode count matrix  ##
## 'read.table' is a cellID x 3 dataframe with Cell, UMI, and Sample columns                                                 ##
## 'read.table' must be pre-subsetted to include only cell barcode reads present in cellIDs                                  ##
## 'cellIDs' is a vector of sequenced cellIDs desired for alignment; 'ref' is a vector of reference sample barcode sequences ##                                                    ##
###############################################################################################################################
MULTIseq.align <- function(readTable, cellIDs, ref) {
  require(stringdist)

  ## Bucket cell IDs
  print("Bucketing cell IDs...")
  cellIDs_bucketed <- bucket_cellIDs(cellIDs)

  ## Bucket read.table
  print("Bucketing read tables...")
  readTable_bucketed <- bucket_readTable(readTable, cellIDs_bucketed)

  ## Initialize output
  bar.table <- as.data.frame(matrix(0L, ncol=length(ref)+2, nrow=length(cellIDs)))
  colnames(bar.table) <- c(paste("Bar", 1:length(ref),sep=""), "nUMI", "nUMI_total")
  rownames(bar.table) <- cellIDs

  ## Iterate through cell ID and readTable buckets
  for (bucket in 1:length(cellIDs_bucketed)) {
    print(paste("Aligning bucket #",bucket,"...",sep=""))
    t0 <- Sys.time()
    print(Sys.time())

    cellIDs.temp <- cellIDs_bucketed[[bucket]]
    readTable.temp <- readTable_bucketed[[bucket]]

    for (cell in cellIDs.temp) {
      ## Extract umis and tags from reads with specific cellID
      r1.ind <- which(readTable.temp$Cell == cell)
      umis <- readTable.temp$UMI[r1.ind]
      tags <- readTable.temp$Sample[r1.ind]

      ## Record total number of reads
      bar.table[cell, "nUMI_total"] <- length(umis)

      ## Find tags with HD <= 1 from any reference tag
      tag.dists <- stringdistmatrix(a=tags, b=ref)

      tag.hds <- apply(tag.dists, 1, min)
      tag.ind <- which(tag.hds <= 1)
      if (length(tag.ind) == 0) { next } ## Skip cellIDs with 0 aligned read (bug fix)
      tag.dists <- tag.dists[tag.ind, ]

      ## Remove duplicated UMIs for reads with reference-matching tags
      umis <- umis[tag.ind]
      umi.ind <- which(duplicated(umis) == FALSE)

      ## Subset by unique UMIs, record total number of unique UMIs, fill count table with HD <= 1
      if (length(umi.ind) == 1) {
        if (is.null(ncol(tag.dists)) == TRUE) { tag.dists <- tag.dists[umi.ind] }
        if (is.null(ncol(tag.dists)) == FALSE) { tag.dists <- tag.dists[umi.ind, ] }
        bar.table[cell,"nUMI"] <- length(umi.ind)

        for (tag in 1:length(tag.dists)) {
          bar.table[cell,tag] <- length(which(tag.dists[tag] <= 1))
        }
      }

      if (length(umi.ind) > 1) {
        tag.dists <- tag.dists[umi.ind, ]
        bar.table[cell,"nUMI"] <- length(umi.ind)

        for (tag in 1:ncol(tag.dists)) {
          bar.table[cell,tag] <- length(which(tag.dists[,tag] <= 1))
        }
      }
    }

    print(paste("Done aligning bucket #",bucket,"...",sep=""))
    print(Sys.time()-t0)

  }

  return(bar.table)

}

##########################################################################################
## bucket_cellIDs is an internal function, sets up 10000-member cellID buckets          ##
## 'cellIDs' is a character vector of cellIDs desired in the final barcode count matrix ##
##########################################################################################
bucket_cellIDs <- function(cellIDs) {

  bucket.ind <- seq(1, length(cellIDs), by=10000)
  cellID.list <- list()

  for (bucket in 1:length(bucket.ind)) {
    if (bucket != length(bucket.ind)) {
      cellIDs.temp <- cellIDs[bucket.ind[bucket]:(bucket.ind[bucket+1]-1)]
      cellID.list[[bucket]] <- cellIDs.temp
    }

    if (bucket == length(bucket.ind)) {
      cellIDs.temp <- cellIDs[bucket.ind[bucket]:length(cellIDs)]
      cellID.list[[bucket]] <- cellIDs.temp
    }
  }

  return(cellID.list)

}

#############################################################################################
## bucket_readTable is an internal function, splits read data according to cellID buckets  ##
## 'readTable' is a parsed read table, as provided by 'MULTIseq.preProcess'                ##
## 'cellIDs_bucketed' is a bucketed list of cellIDs, as provided by 'bucket_cellIDs'       ##
#############################################################################################
bucket_readTable <- function(readTable, cellIDs_bucketed) {

  readTable.list <- list()

  for (bucket in 1:length(cellIDs_bucketed)) {
    cellIDs.temp <- cellIDs_bucketed[[bucket]]
    ind <- which(readTable$Cell %in% cellIDs.temp)
    readTable.list[[bucket]] <- readTable[ind, ]
  }

  return(readTable.list)

}

#########################################################################################################
## MULTIseq.preProcess reads in raw barcode FASTQs and allocates reads into cell barcode,              ##
## UMI, and sample barcode subsets. Function parses FASTQs to only include desired cellIDs             ##
## 'R1' is an R1 FASTQ file; 'R2' is an R2 FASTQ file; 'cellIDs' is a character vecotro of cellIDs     ##
#########################################################################################################
MULTIseq.preProcess <- function(R1, R2, cellIDs, chemistry = "V3") {
  require(ShortRead)

  print("Reading in R1...")
  r1 <- readFastq(R1)
  gc()

  print("Reading in R2...")
  r2 <- readFastq(R2)
  gc()

  print("Assembling read table...")
  if (chemistry == "V2") {
    readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
                       as.data.frame(subseq(sread(r1),17,26)),
                       as.data.frame(subseq(sread(r2),1,8)))
    colnames(readTable) <- c("Cell","UMI","Sample")
    gc()
    ind <- which(readTable$Cell %in% cellIDs)
    readTable <- readTable[ind, ]
    return(readTable)
  }

  if (chemistry == "V3") {
    readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
                       as.data.frame(subseq(sread(r1),17,28)),
                       as.data.frame(subseq(sread(r2),1,8)))
    colnames(readTable) <- c("Cell","UMI","Sample")
    gc()
    ind <- which(readTable$Cell %in% cellIDs)
    readTable <- readTable[ind, ]
    return(readTable)
  }

}

#########################################################################################################
## MULTIseq.preProcess_BD reads in raw barcode FASTQs and allocates reads into cell barcode,           ##
## UMI, and sample barcode subsets. Function parses FASTQs to only include desired cellIDs             ##
## 'R1' is an R1 FASTQ file; 'R2' is an R2 FASTQ file; 'cellIDs' is a character vecotro of cellIDs     ##
#########################################################################################################

# MULTIseq.preProcess_BD <- function(R1, R2, cellIDs) {
#   require(ShortRead)
#
#   print("Reading in R1...")
#   r1 <- readFastq(R1)
#   gc()
#   print("Reading in R2...")
#   r2 <- readFastq(R2)
#   gc()
#   print("Assembling read table...")
#   readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
#                      as.data.frame(subseq(sread(r1),17,26)),
#                      as.data.frame(subseq(sread(r2),1,75)))
#   colnames(readTable) <- c("Cell","UMI","Sample")
#   gc()
#   ind <- which(readTable$Cell %in% cellIDs)
#   readTable <- readTable[ind, ]
#
#   return(readTable)
# }


#########################################################################################################
## MULTIseq.preProcess_allCells reads in raw barcode FASTQs and allocates reads into cell barcode,     ##
## UMI, and sample barcode subsets Function buckets data for pre-processing all present barcodes       ##
## 'R1' is an R1 FASTQ file; 'R2' is an R2 FASTQ file; 'cellIDs' is the 10X cell barcode whitelist     ##
#########################################################################################################
MULTIseq.preProcess_allCells <- function(R1, R2, whitelist, chemistry = "V3") {
  require(ShortRead)

  print("Reading in R1...")
  r1 <- readFastq(R1)
  gc()

  print("Reading in R2...")
  r2 <- readFastq(R2)
  gc()

  print("Assembling read table...")
  if (chemistry == "V2") {
    readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
                       as.data.frame(subseq(sread(r1),17,26)),
                       as.data.frame(subseq(sread(r2),1,8)))
    colnames(readTable) <- c("Cell","UMI","Sample")
    gc()
  }

  if (chemistry == "V3") {
    readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
                       as.data.frame(subseq(sread(r1),17,28)),
                       as.data.frame(subseq(sread(r2),1,8)))
    colnames(readTable) <- c("Cell","UMI","Sample")
    gc()
  }

  print("Bucketing cell IDs...")
  cellIDs_bucketed <- bucket_cellIDs(whitelist)

  ## Initialize logicals for storing present Cell IDs, cell-associated reads
  present.cells <- logical(length(whitelist))
  cell.linked.reads <- logical(nrow(readTable))

  ## Iterate through cell ID buckets, finding reads associated with each cell
  x <- length(cellIDs_bucketed)
  for (bucket in 1:x) {
    print(paste("Starting bucket ",bucket,"/",x,sep=""))
    t0 <- Sys.time()
    print(Sys.time())

    cellID.counter <- 10000*(bucket-1)
    cellIDs.temp <- cellIDs_bucketed[[bucket]]

    read.ind <- which(readTable$Cell %in% whitelist)
    cell.linked.reads[read.ind] <- TRUE

    cell.ind <- which(whitelist %in% readTable$Cell)
    present.cells[(cell.ind+cellID.counter)] <- TRUE
  }

  print("Parsing final cell IDs...")
  cellIDs.parsed <- whitelist[present.cells]
  assign(x="cellIDs.parsed", value=cellIDs.parsed, envir = .GlobalEnv)

  print("Parsing final read table...")
  readTable.parsed <- readTable[cell.linked.reads, ]
  return(readTable.parsed)
}

#########################################################################################################
## MULTIseq.preProcess_allCells_BD reads in raw barcode FASTQs and allocates reads into cell barcode,  ##
## UMI, and sample barcode subsets Function buckets data for pre-processing all present barcodes       ##
## 'R1' is an R1 FASTQ file; 'R2' is an R2 FASTQ file; 'cellIDs' is the 10X cell barcode whitelist     ##
#########################################################################################################

# MULTIseq.preProcess_allCells_BD <- function(R1, R2, whitelist) {
#   require(ShortRead)
#
#   print("Reading in R1...")
#   r1 <- readFastq(R1)
#   print("Reading in R2...")
#   r2 <- readFastq(R2)
#   print("Assembling read table...")
#   readTable <- cbind(as.data.frame(subseq(sread(r1),1,16)),
#                      as.data.frame(subseq(sread(r1),17,26)),
#                      as.data.frame(subseq(sread(r2),1,96)))
#   colnames(readTable) <- c("Cell","UMI","Sample")
#
#   print("Bucketing cell IDs...")
#   cellIDs_bucketed <- bucket_cellIDs(whitelist)
#
#   ## Initialize logicals for storing present Cell IDs, cell-associated reads
#   present.cells <- logical(length(whitelist))
#   cell.linked.reads <- logical(nrow(readTable))
#
#   ## Iterate through cell ID buckets, finding reads associated with each cell
#   x <- length(cellIDs_bucketed)
#   for (bucket in 1:x) {
#     print(paste("Starting bucket ",bucket,"/",x,sep=""))
#     t0 <- Sys.time()
#     print(Sys.time())
#
#     cellID.counter <- 10000*(bucket-1)
#     cellIDs.temp <- cellIDs_bucketed[[bucket]]
#
#     read.ind <- which(readTable$Cell %in% whitelist)
#     cell.linked.reads[read.ind] <- TRUE
#
#     cell.ind <- which(whitelist %in% readTable$Cell)
#     present.cells[(cell.ind+cellID.counter)] <- TRUE
#   }
#
#   print("Parsing final cell IDs...")
#   cellIDs.parsed <- whitelist[present.cells]
#   assign(x="cellIDs.parsed", value=cellIDs.parsed, envir = .GlobalEnv)
#
#   print("Parsing final read table...")
#   readTable.parsed <- readTable[cell.linked.reads, ]
#   return(readTable.parsed)
# }



##############################################################################
## 'alignRate' computes the proportion of reads that align to the reference ##
##############################################################################
alignRate <- function(readTable, cellIDs, ref) {
  require(stringdist)

  print("Subsetting readTable...")
  ind <- which(readTable$Cell %in% cellIDs)
  readTable <- readTable[ind, ]

  ## Compute align rate
  print("Computing alignment rate...")
  align.rate <- rep(0L, length(cellIDs))
  counter <- 0
  for (cell in cellIDs) {
    counter <- counter + 1
    print(counter)
    ## Extract umis and tags from reads with specific cellID
    r1.ind <- which(readTable$Cell == cell)
    umis <- readTable$UMI[r1.ind]
    tags <- readTable$Sample[r1.ind]

    ## Compute HD for all tags relative to reference
    tag.dists <- stringdistmatrix(a=tags, b=ref)

    ## Find rate at which reads align to reference
    tag.hds <- apply(tag.dists, 1, min)
    x <- (length(which(tag.hds <= 1))/length(tag.hds))*100
    align.rate[counter] <- x
  }

  return(align.rate)

}

