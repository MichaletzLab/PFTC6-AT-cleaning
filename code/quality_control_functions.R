# Helper functions for Assimilation-Temperature response quality control
# Plant Functional Traits Course 6 - Western Norway
# Last updated, 9 April 2025, Josef Garen



# Function to remove multimodal curves from dataset (i.e. curves with multiple "optima")
# Argument "data" is LI-6800 gas exchange dataframe
cut.multimodal = function(data) {

  data$no_modes = NA # Initialize holder for number of modes (optima)

  # Loop over curves
  for(i in unique(data$curveID)) {

    # Grab one curve
    d1 = subset(data, curveID == i)

    # Fit a cubic spline regression to the A-T data
    z = smooth.spline(x = d1$Tleaf, y = d1$A, nknots = 10)

    # Identify local maxima numerically by finding places where first derivative switches from positive to negative
    Sx = seq(min(d1$Tleaf)+1, max(d1$Tleaf)-1, 0.01)
    len = length(Sx)
    pred.spline <- predict(z, Sx)
    d0 <- data.frame(pred.spline)
    pred.prime <- predict(z, Sx, deriv=1)
    d0$local_max = c(pred.prime$y[1:len-1] > 0 & pred.prime$y[2:len] < 0, F)

    # Count local maxima
    Maxes = which(d0$local_max == T)
    data[data$curveID == i,]$no_modes = length(Maxes)
  }

  # Remove curves with more than one mode (or zero modes), then return
  data.cut = subset(data, no_modes == 1)
  return(data.cut)

}



# Function to remove data w/ poor r-squared or with fitted T_opt too close to Tmin or Tmax
# Argument "data" is LI-6800 gas exchange dataframe
cut.topt.out.of.bounds = function(data) {

  # Grab only columns of interest
  data2 = data %>% select(Tleaf, A, curveID)

  # Fit modified Sharpe-Schoolfield model to all curves (time consuming)
  results = fit_curves_all_with_Amax(data2, x = "Tleaf", y = "A")

  # For each curve, check that T_opt is within bounds, and that r_sq is above 0.5
  results$keep = NA
  for (i in unique(data$curveID)) {
    curdat = subset(data, curveID == i)
    Tmin = min(curdat$Tleaf)
    Tmax = max(curdat$Tleaf)
    j = which(results$curveID == i)
    Topt = results[j,]$topt
    r_sq = results[j,]$r_sq

    results[j,]$keep = Topt > Tmin+3 & Topt < Tmax-3 & r_sq > 0.5
  }

  # Remove curves that don't meet quality criteria
  results = subset(results, keep == T)

  # Return subsetted data
  data = subset(data, curveID %in% results$curveID)
  return(data)
}



# New formulation of Sharpe-Schoolfield with explicit Amax, see Garen & Michaletz (2024) for derivation
# Arguments:
# temp = temperature (C)
# r_max = maximum rate, i.e. rate at topt (arbitrary units)
# e = activation energy (eV)
# eh = deactivation energy (eV)
# topt = optimum temperature for rate (C)
ss_with_amax = function (temp, r_max, e, eh, topt) {
  k <- 8.62e-05
  boltzmann.term <- r_max * (1 + (e/(eh-e))) * exp(e/k * (1/(topt+273.15) - 1/(temp + 273.15)))
  inactivation.term <- 1/(1 + (e/(eh - e)) * exp(eh/k * (1/(topt + 273.15) - 1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}



# Function to fit many curves with modified Sharpe-Schoolfield model. Returns dataframe of best fit parameters
# Arguments:
# Data = Dataframe containing data to be fit with "x" and "y" columns, plus a curveID column
# x = name of "x" column in Data
# y = name of "y" column in Data
fit_curves_all_with_Amax <- function(Data, x = "Tleaf", y = "Photo") {

  set.seed(1) # Set RNG seed for consistency across multiple runs

  # Standardize names of x and y for use throughout
  if(x != "Tleaf") {
    Data$Tleaf = Data[[x]] # Rename X as Tleaf
    Data[[x]] = NULL      # Get rid of originals
  }
  if(y != "Photo") {
    Data$Photo = Data[[y]] # Rename Y as Photo
    Data[[y]] = NULL      # Get rid of originals
  }

  # Remove negative Photo
  Data = subset(Data, Photo > 0)

  Data = Data %>% arrange(by_group = curveID)
  results = c()

  ######################
  # Curve fitting loop #
  ######################
  message("Fitting curves; this may take some time.")
  n_iter <- length(unique(Data$curveID))
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar

  for (i in unique(Data$curveID)) {

    d_1 <- subset(Data, curveID == i)    # Grab only current curveID

    # Fit Schoolfield function using nls_multstart
    fit <- nls_multstart(Photo ~ ss_with_amax(temp = Tleaf, r_max, e, eh, topt),
                         data = d_1,
                         iter = 1000,
                         start_lower = c(r_max = 0, e = 0, eh = 0.2, topt = 0),
                         start_upper = c(r_max = 20, e = 2, eh = 5, topt = 50),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         lower = c(r_max = 0, e = 0, eh = 0, topt = 0))

    # If the fit was successful, extract parameters estimates and SE
    if(typeof(fit) != "NULL") {
      r_max <- coef(fit)["r_max"]
      r_max_SE <- summary(fit)$coefficients[,'Std. Error']['r_max']
      e <- coef(fit)["e"]
      e_SE <- summary(fit)$coefficients[,'Std. Error']['e']
      eh <- coef(fit)["eh"]
      eh_SE <- summary(fit)$coefficients[,'Std. Error']['eh']
      topt <- coef(fit)["topt"]
      topt_SE <- summary(fit)$coefficients[,'Std. Error']['topt']
      AIC <- AIC(fit)
      r_sq = 1-sum(resid(fit)^2)/sum((d_1$Photo-mean(d_1$Photo))^2)

      breadth = get_breadth(fit, level = 0.8)

    # Otherwise, likely a convergence failure. All set to NA
    } else {
      r_max = r_max_SE = e = e_SE = eh = eh_SE = topt = topt_SE = breadth = AIC = r_sq = NA
      d_1$failure_status = "Convergence failure"
    }

    # Increment progress bar
    pb$tick()

    # Build dataframe with results
    results = rbind(results, data.frame(d_1[1,] %>% select(-Tleaf, -Photo), r_max, r_max_SE, e, e_SE,
                                        eh, eh_SE, topt, topt_SE, breadth, AIC, r_sq))
  }

  # Return results
  return(results)
}
