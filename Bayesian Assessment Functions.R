#Model assessment functions
#calculate the WAIC
#From: https://github.com/stan-dev/stan/issues/473
# R script for Waic example.  Also needed in this directory:  data file hibbs.dat and Stan file lm_waic.stan

# Little function to calculate posterior variances from simulation
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function(stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              pointwise=pointwise, total=total, se=se))
}

#Function from ggmcmc, but not in the package - copied it from github, it calculates the percent of correct predictions
ggs_pcp <- function(D, outcome, threshold = "observed", bins = 30) {
  # Clarify the threshold
  if (threshold == "observed") {
    threshold <- length(which(outcome == 1)) / length(outcome)
  } else {
    if (!is.numeric(threshold) | (is.numeric(threshold) & (threshold < 0 | threshold > 1))) {
      stop("Error: threshold must be either 'observed' or a number between 0 and 1")
    }
    threshold <- threshold
  }
  # Calculate the percent correctly predicted
  S <- dplyr::inner_join(D, dplyr::data_frame(Observed=outcome, Parameter=unique(D$Parameter)), by="Parameter") %>%
    dplyr::mutate(Correct = ifelse( (value < threshold & Observed == 0) | # 0 predicted and 0 observed
                                      (value > threshold & Observed == 1), # 1 predicted and 1 observed
                                    TRUE, FALSE)) %>%
    dplyr::group_by(Iteration, Chain) %>%
    dplyr::summarize(PCP = length(which(Correct)) / dplyr::n())
  # Adjust binwidth
  pcp.bw <- calc_bin(S$PCP, bins = bins)
  names(pcp.bw)[names(pcp.bw)=="x"] <- "Percent_correctly_predicted"
  # Plot
  f <- ggplot(pcp.bw, aes(x=Percent_correctly_predicted, y = count, width = width)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_vline(xintercept = mean(pcp.bw$Percent_correctly_predicted), color="red")+
    annotate("text", label = paste(round(mean(pcp.bw$Percent_correctly_predicted), 2)), x = mean(pcp.bw$Percent_correctly_predicted)+.08, y = 1500, size = 8, colour = "red")+
    expand_limits(x = c(0, 1))+
    labs(x = "Percent Correctly Predicted", y = "Count")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(f)
}

