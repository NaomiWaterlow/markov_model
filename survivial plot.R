plot.msm <- function (x, from = NULL, to = NULL, range = NULL, covariates = "mean", 
                      legend.pos = NULL, xlab = "Time", ylab = "Fitted survival probability", 
                      lwd = 1, ...) 
{
  if (!inherits(x, "msm")) 
    stop("expected x to be a msm model")
  if (is.null(from)) 
    from <- transient.msm(x)
  else {
    if (!is.numeric(from)) 
      stop("from must be numeric")
    if (any(!(from %in% 1:x$qmodel$nstates))) 
      stop("from must be a vector of states in 1, ..., ", 
           x$qmodel$nstates)
  }
  if (is.null(to)) {
    if (length(absorbing.msm(x)) == 0) 
      stop("\"to\" not specified, and no absorbing state in the model")
    to <- max(absorbing.msm(x))
  }
  else {
    if (!is.numeric(to)) 
      stop("to must be numeric")
    if (!(to %in% absorbing.msm(x))) 
      stop("to must be an absorbing state")
  }
  if (is.null(range)) 
    rg <- range(model.extract(x$data$mf, "time"))
  else {
    if (!is.numeric(range) || length(range) != 2) 
      stop("range must be a numeric vector of two elements")
    rg <- range
  }
  timediff <- (rg[2] - rg[1])/50
  times <- seq(rg[1], rg[2], timediff)
  pr <- numeric()
  cols <- rainbow(length(from))
  for (t in times) pr <- c(pr, pmatrix.msm(x, t, times[1], 
                                           covariates)[from[1], to])

  df <- data.frame(times = times, "From state SS" = 1-pr)
  
  for (st in from[-1]) {
    pr <- numeric()
    for (t in times) pr <- c(pr, pmatrix.msm(x, t, times[1], 
                                             covariates)[st, to])
    df[,st+1] <- 1-pr
colnames(df)[st+1] <- paste0("From state ",names(from)[st])
  }
  colnames(df)[2] <- "From state SS"

   
df_m <- melt(df, id.vars= "times")

SURVIVAL <- ggplot(df_m, aes(x = times, y = value, colour = variable)) + geom_line() + 
  theme_linedraw() + labs(x = "Days from start", y = "Fitted survival probability", 
                          colour= "transition to RR", title = "B")  +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15))
return(SURVIVAL)

}
