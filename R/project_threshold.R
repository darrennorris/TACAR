
#' Project density dependent population dynamics
#'
#' @description
#' Project density dependent dynamics of a population matrix
#' projection model.
#'
#' @param A A matrix, or list of matrices.
#' @param vector A numeric vector describiing the age/stage distribution.
#' @param time The number of projection intervals.
#' @param return.vec If TRUE, returns the time series of demographic (st)age
#' vectors.
#' @param threshold The population threshold.
#' @param A_DD A single high density population matrix.
#'
#' @details
#' A version of popbio::project that will project forward
#' until reaching a threshold population size, at which point
#' it will switch to a specified high-density population matrix.
#'
#' @author Chrissy Hernandez <chrissy.hernandez@biology.ox.ac.uk>
#' @returns A list with projection results.
#'
#' @export project_threshold
#'
#' @examples
#' \dontrun{
#' # Try out the DD projection function when we use a single matrix:
#' nohunt_threshold<- project_threshold(pop_mat_100_none, vector=pop_n,
#'                                     time=100, threshold=1000, A_DD=nohunt_highD)
#'}
#'
project_threshold<- function(A, vector=NULL, time=100, return.vec=TRUE, threshold=NULL, A_DD=NULL){
  # This function will project forward according to the sequence of matrices in
  # A until reaching a threshold population size, at which point it will switch
  # to a specified high-density population matrix

  if (is.list(A) & length(A) == 1) {
    A <- A[[1]] }
  if (is.matrix(A)) {
    M1 <- A
    dim(A) <- c(dim(A), 1)
    dimnames(A)[[1]] <- dimnames(M1)[[1]]
    dimnames(A)[[2]] <- dimnames(M1)[[2]]
    dimnames(A)[[3]] <- NULL
  }
  if (is.list(A) & length(A) > 1) {
    numA <- length(A)
    alldim <- sapply(A, dim)
    if (!diff(range(alldim)) == 0) {
      stop("all matrices in A must be square and have the same dimension as each other")
    }
    dimA <- mean(alldim)
    L <- A
    A <- numeric(dimA * dimA * numA)
    dim(A) <- c(dimA, dimA, numA)
    for (i in 1:numA) {
      A[, , i] <- L[[i]]
    }
    dimnames(A)[[1]] <- dimnames(L[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(L[[1]])[[2]]
    dimnames(A)[[3]] <- names(L)
  }

  order <- dim(A)[1]
  stagenames <- dimnames(A)[[2]]
  if (is.null(stagenames)) {
    stagenames <- paste("S", as.character(1:order), sep = "") }
  nmat <- dim(A)[3]

  # If we're using the DD part of this function:
  if (!is.null(threshold)){
    if (is.null(A_DD)){
      stop("If you supply a population threshold, then you must supply a high-density projection matrix.") }
    # For now, assume that we'll supply a single DD matrix - can adjust to work
    # with a sequence of matrices to enable stochastic simulations
    if (is.list(A_DD) & length(A_DD) == 1) {
      A_DD <- A_DD[[1]] }
    if (is.matrix(A_DD)) {
      M1 <- A_DD
      dim(A_DD) <- c(dim(A_DD), 1)
      dimnames(A_DD)[[1]] <- dimnames(M1)[[1]]
      dimnames(A_DD)[[2]] <- dimnames(M1)[[2]]
      dimnames(A_DD)[[3]] <- NULL
    }
    if (!is.list(A_DD)){
      if (dim(A_DD)[1]!=dim(A_DD)[2]){
        stop("The high-density projection matrix must be square.")
      }
      order_DD<- dim(A_DD)[1]
      if (order_DD!=order){
        stop("The high-density projection matrix must be the same order as the matrices in A.")
      }
    } else {
      numA <- length(A_DD)
      alldim <- sapply(A_DD, dim)
      if (!diff(range(alldim)) == 0) {
        stop("all matrices in A_DD must be square and have the same dimension as each other")
      }
      dimA <- mean(alldim)
      L <- A_DD
      A_DD <- numeric(dimA * dimA * numA)
      dim(A_DD) <- c(dimA, dimA, numA)
      for (i in 1:numA) {
        A_DD[, , i] <- L[[i]]
      }
      dimnames(A_DD)[[1]] <- dimnames(L[[1]])[[1]]
      dimnames(A_DD)[[2]] <- dimnames(L[[1]])[[2]]
      dimnames(A_DD)[[3]] <- names(L)
    }


  }

  # Initialize an output object:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())

  # Initial vector:
  if (is.null(vector)){
    n0<- rep(1/order, order)
  } else{
    n0 <- vector
  }

  # Set the order that the matrices will be applied:
  if (dim(A)[3]==1){
    MC<- rep(1, time)
  } else if (dim(A)[3]>=time){
    MC<- 1:time
  } else { # if A has fewer matrices than the requested number of timesteps,
    # then loop through them until you reach the number of requested timesteps
    MC<- rep(1:nmat, ceiling(time/nmat))[1:time]
  }

  # Set the order that the DD matrices will be applied:
  if (dim(A_DD)[3]==1){
    MC_DD<- rep(1, time)
  } else if (dim(A_DD)[3]>=time){
    MC_DD<- 1:time
  } else { # if A has fewer matrices than the requested number of timesteps,
    # then loop through them until you reach the number of requested timesteps
    nmat_DD<- dim(A_DD)[3]
    MC_DD<- rep(1:nmat_DD, ceiling(time/nmat_DD))[1:time]
  }

  # Initialize the output for population vector:
  Vec <- matrix(0, ncol = order, nrow = time + 1)
  Pop <- numeric(time + 1)
  dimnames(Vec)[[2]] <- stagenames
  Vec[1, ] <- vector
  Pop[1] <- sum(vector)
  for (i in 1:time) {
    if (Pop[i]>threshold & dim(A_DD)[3]==1){
      Vec[(i + 1), ] <- A_DD[,,1] %*% Vec[i, ]
      Pop[i + 1] <- sum(Vec[(i + 1), ])
    } else if (Pop[i]>threshold & dim(A_DD)[3]>1){
      Vec[(i+1), ]<- A_DD[,,MC_DD[i]] %*% Vec[i, ]
      Pop[i + 1] <- sum(Vec[(i + 1), ])
    }  else {
      Vec[(i + 1), ] <- A[, , MC[i]] %*% Vec[i, ]
      Pop[i + 1] <- sum(Vec[(i + 1), ])
    }
  }
  if (is.null(dim(Pop)))
    dim(Pop) <- time + 1
  out$pop <- Pop
  if (return.vec)
    out$vec <- Vec
  out$mat <- A

  return(out)

}
