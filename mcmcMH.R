mcmcMH<- function (target, init.theta, proposal.sd = NULL, n.iterations, 
                   covmat = NULL, limits = list(lower = NULL, upper = NULL), 
                   adapt.size.start = NULL, adapt.size.cooling = 0.99, adapt.shape.start = NULL, 
                   adapt.shape.stop = NULL, print.info.every = n.iterations/1000, 
                   verbose = FALSE, max.scaling.sd = 50) 
{
  theta.current <- init.theta
  theta.propose <- init.theta
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  theta.names <- names(init.theta)
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  if (is.null(covmat.proposal)) {
    if (is.null(proposal.sd)) {
      proposal.sd <- init.theta/10
    }
    covmat.proposal <- matrix(diag(proposal.sd[theta.names]^2, 
                                   nrow = length(theta.names)), nrow = length(theta.names), 
                              dimnames = list(theta.names, theta.names))
  }
  else {
    covmat.proposal <- covmat.proposal[theta.names, theta.names]
  }
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  }
  else {
    lower.proposal <- lower.proposal[theta.names]
  }
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  }
  else {
    upper.proposal <- upper.proposal[theta.names]
  }
  covmat.proposal.init <- covmat.proposal
  adapting.size <- FALSE
  adapting.shape <- 0
  theta.estimated.names <- names(which(diag(covmat.proposal) > 
                                         0))
  target.theta.current <- target(theta.current)
  if (class(target.theta.current) == "numeric") {
    target.theta.current <- list(log.density = target.theta.current, 
                                 trace = theta.current)
  }
  if (!is.null(print.info.every)) {
    message(Sys.time(), ", Init: ", printNamedVector(theta.current[theta.estimated.names]), 
            ", target: ", target.theta.current[["log.density"]])
  }
  trace <- matrix(ncol = length(target.theta.current[["trace"]]) + 
                    2, nrow = n.iterations, 0)
  colnames(trace) <- c(theta.estimated.names, "log.density", "Accepted")
  acceptance.rate <- 0
  scaling.sd <- 1
  scaling.multiplier <- 1
  covmat.empirical <- covmat.proposal
  covmat.empirical[, ] <- 0
  theta.mean <- theta.current
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }
  start_iteration_time <- Sys.time()
  for (i.iteration in seq_len(n.iterations)) {
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start && 
        (is.null(adapt.shape.start) || acceptance.rate * 
         i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration - 
                                                      adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd, max.scaling.sd))
      covmat.proposal.new <- scaling.sd^2 * covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] < 
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
    }
    else if (!is.null(adapt.shape.start) && acceptance.rate * 
             i.iteration >= adapt.shape.start && (adapting.shape == 
                                                  0 || is.null(adapt.shape.stop) || i.iteration < adapting.shape + 
                                                  adapt.shape.stop)) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        adapting.shape <- i.iteration
      }
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    }
    else if (adapting.shape > 0) {
      message("\n---> Stop adapting shape of covariance matrix")
      adapting.shape <- -1
    }
    if (i.iteration%%ceiling(print.info.every) == 0) {
      state.mcmc <- target.theta.current$trace
      message(Sys.time(), ", Iteration: ", i.iteration, 
              "/", n.iterations, ", acceptance rate: ",
              sprintf("%.3f", acceptance.rate), appendLF = FALSE)
      seq_save <- seq(from=0,to=n.iterations, length.out = 5)

      
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd), 
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier), 
                appendLF = FALSE)
      }
      message(", state: ", printNamedVector(state.mcmc))
      message(", logdensity: ", target.theta.current$log.density)
      
    }
    if (any(diag(covmat.proposal)[theta.estimated.names] < 
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names, theta.estimated.names])
      stop("non-positive definite covmat", call. = FALSE)
    }
    if (length(theta.estimated.names) > 0) {

      theta.propose[theta.estimated.names] <- as.vector(rtmvnorm(1, 
                                                                 mean = theta.current[theta.estimated.names], 
                                                                 sigma = covmat.proposal,
                                                                 lower = lower.proposal[theta.estimated.names], 
                                                                 upper = upper.proposal[theta.estimated.names]))
    }

    target.theta.propose <- target(theta.propose)
    if (class(target.theta.propose) == "numeric") {
      target.theta.propose <- list(log.density = target.theta.propose, 
                                   trace = theta.propose)
    }
    if (!is.finite(target.theta.propose$log.density)) {
      log.acceptance <- -Inf
    }
    else {
      log.acceptance <- target.theta.propose$log.density - 
        target.theta.current$log.density
      log.acceptance <- log.acceptance + dtmvnorm(x = theta.current[theta.estimated.names],
                                                  mean = theta.propose[theta.estimated.names],
                                                  sigma = covmat.proposal
                                                  , lower = lower.proposal[theta.estimated.names],
                                                  upper = upper.proposal[theta.estimated.names],
                                                  log = TRUE)
      log.acceptance <- log.acceptance - dtmvnorm(x = theta.propose[theta.estimated.names],
                                                  mean = theta.current[theta.estimated.names],
                                                  sigma = covmat.proposal, lower = lower.proposal[theta.estimated.names],
                                                  upper = upper.proposal[theta.estimated.names],
                                                  log = TRUE)
     }
      
      
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names], 
              ", target: ", target.theta.propose[["log.density"]], 
              ", acc prob: ", exp(log.acceptance), ", ", appendLF = FALSE)
    }

    if (is.accepted <- (log(runif(1)) < log.acceptance)) {
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    }
    else if (verbose) {
      message("rejected")
    }
    trace[i.iteration, ] <- c(target.theta.current[["trace"]], 
                              target.theta.current[["log.density"]],
                              is.accepted)
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    }
    
    else {
      acceptance.rate <- acceptance.rate + (is.accepted - 
                                              acceptance.rate)/i.iteration
    }
    
    # else if (i.iteration < 1000) {
    #   acceptance.rate <- sum(trace[,"Accepted"])/i.iteration
    #     #acceptance.rate + (is.accepted - acceptance.rate)/i.iteration
    # }
    # else if (i.iteration > 1000){
    #   acceptance.rate <- sum(trace[(i.iteration-1000):i.iteration,"Accepted"])/1000
    # }
    if (adapting.shape >= 0) {
      tmp <- updateCovmat(covmat.empirical, theta.mean, 
                          theta.current, i.iteration)
      covmat.empirical <- tmp$covmat
      theta.mean <- tmp$theta.mean
    }
  }
  return(list(trace = trace, acceptance.rate = acceptance.rate, 
              covmat.empirical = covmat.empirical))
}

updateCovmat<- function (covmat, theta.mean, theta, i) 
{
  if (is.null(names(theta))) {
    stop("Argument ", sQuote("theta"), " must be named.", 
         .call = FALSE)
  }
  if (is.null(names(theta.mean))) {
    stop("Argument ", sQuote("theta.mean"), " must be named.", 
         .call = FALSE)
  }
  if (is.null(rownames(covmat))) {
    stop("Argument ", sQuote("covmat"), " must have named rows.", 
         .call = FALSE)
  }
  if (is.null(colnames(covmat))) {
    stop("Argument ", sQuote("covmat"), " must have named columns.", 
         .call = FALSE)
  }
  covmat <- covmat[names(theta), names(theta)]
  theta.mean <- theta.mean[names(theta)]
  residual <- as.vector(theta - theta.mean)
  covmat <- (covmat * (i - 1) + (i - 1)/i * residual %*% t(residual))/i
  theta.mean <- theta.mean + residual/i
  return(list(covmat = covmat, theta.mean = theta.mean))
}

printNamedVector <- function (x, fmt = "%.2f", sep = " | ")
{
  paste(paste(names(x), sprintf(fmt, x), sep = " = "), collapse = sep)
}
