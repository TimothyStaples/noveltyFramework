novelty <- function(x, refSet, method, nSize = 1, # basic arguments
                    largeMatSwitch = TRUE, largeMatThreshold = 5000, largeMatCores=1, largeMatnSize = 10, # large data options
                    ...){

  require(vegan)

  if(class(refSet)[1] != "matrix"){

    error.trap <- tryCatch(as.matrix(refSet), error=function(e) e)
    if(inherits(error.trap, "error")){
      stop("Reference set must be a matrix or matrix-like object")
    }

    refSet <- as.matrix(refSet)}

  if(nrow(refSet)==1){
    warning("Reference set is of size 1: novelty will be nonsensical")
    }

  if(class(x)[1] != "matrix"){
    x = matrix(x, ncol=length(x))
  }

  if(ncol(x) != ncol(refSet)){
    stop("Target and reference set must have the same number of columns")
  }

  if(!mode(x) %in% c("numeric", "double", "logical") |
     !mode(refSet) %in% c("numeric", "double", "logical")){
  stop("Dissimilarity can only be calculated using binary or numeric values")
  }

  if(nSize > nrow(refSet)){
    warning(paste0("nSize has been corrected to ", nrow(refSet), " (size of reference set)"))
    nSize <- nrow(refSet)
  }

  # set up noveldist function

  if(largeMatSwitch & (nrow(refSet) + ifelse(is.null(dim(x)), 1, nrow(x)) > largeMatThreshold)){
    message("Implementing truncated calculation for large data")
    require(parallel)

    xMat <- do.call("rbind", mclapply(1:nrow(x), function(n){
              tar <- x[n,]
              ranks <- order(rowSums(t(abs(t(refSet)-tar))))

              suppressMessages(
              sort(as.matrix(noveldist(rbind(tar, refSet[ranks,][1:largeMatnSize,]),
                                       method=method, ...))[,1])[2:(nSize+1)]
              )

    }, mc.cores=largeMatCores))

    if(ncol(xMat)==1){return(as.vector(xMat))}
    return(xMat)
    }

  xMat <- suppressMessages(
              as.matrix(noveldist(rbind(x, refSet), method=method, ...))[-(1:nrow(x)), 1:nrow(x)]
  )

  if(class(xMat)[1]=="matrix"){

    return(t(apply(xMat, 2, sort)[1:nSize, 1:nrow(x)]))
  } else {
    return(sort(xMat)[1:nSize])

  }

  return(t(xMat))

  }
