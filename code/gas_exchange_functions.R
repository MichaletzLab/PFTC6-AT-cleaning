# Helper functions for reading and correcting LI-6800 gas exchange files
# Plant Functional Traits Course 6 - Western Norway
# Last updated, 9 April 2025, Josef Garen



# This function reads in a specified LI-6800 file (raw format, not excel)
# The function additionally reads in the coefficients used for calculating the
# boundary layer conductance (BLC); this is later used to recalculate BLC, to correct
# an error in certain versions of the LI-6800 firmware. This also extracts energy
# balance parameters from the header so that Tleaf can be recalculated from EB.
# Argument "filename" is the name of a raw LI-6800 file
read_6800_with_BLC <- function(filename) {

  raw_input = readLines(filename)  # Read in raw datafile

  # Grab the BLC parameters
  blc_data = data.frame(raw = raw_input[grep("blc_a", raw_input):grep("blc_Po", raw_input)]) %>%
    separate(raw, into = c(NA, "blc_params"), sep = ":") %>%
    separate(blc_params, into = c("name", "value"), sep = "\t")
  blc_data$value = as.numeric(blc_data$value)
  blc_data = blc_data %>% spread(name, value)

  # Grab energy balance parameters
  eb_data = data.frame(raw = raw_input[grep("LTConst", raw_input)]) %>%
    separate(raw, into = c(NA, "EBparams"), sep = ":") %>%
    separate(EBparams, into = c("name", "value"), sep = "\t")
  eb_data$value = as.numeric(eb_data$value)
  eb_data = eb_data %>% spread(name, value)

  data_start = grep("\\[Data\\]", raw_input) + 2 # Find where data starts (at "[Data]" line)
  data_end = grep("\\[Header\\]", raw_input) - 1 # Find where data ends (at "[Header]" line)
  data_end = c(data_end, length(raw_input))[-1]  # Add data end at end of file

  data_end = data_end[!(data_end < data_start[1])] # Find real data endings

  compiled_data = c()                                           # Initiate blank holder

  # Loop over each data set in the file
  for (i in 1:length(data_start)) {
    trimmed = raw_input[data_start[i]:data_end[i]]              # Grab data

    if(!is_empty(grep("Stability Definition", trimmed))) {
      trimmed = trimmed[-grep("Stability Definition", trimmed)] # Remove stability definitions
    }
    if(!is_empty(grep("^[^\t]*\t?[^\t]*$", trimmed))) {
      trimmed = trimmed[-grep("^[^\t]*\t?[^\t]*$", trimmed)]    # Remove comments
    }

    trimmed = trimmed[-2]                                       # Remove units line
    current_data = read.csv(text = trimmed, sep = "\t")         # Convert to dataframe
    compiled_data = rbind(compiled_data, current_data)          # Merge together into one
  }
  compiled_data$filename = filename                             # Add filename column

  compiled_data = bind_cols(compiled_data, blc_data, eb_data)   # Add BLC and EB data as columns instead of headers

  return(compiled_data)
}



# This function recalculates all (or most) "derived" gas exchange variables in LI-6800 data.
# Formulae are taken from LI-6800 manual and/or excel sheet outputs. Recalculation of boundary
# layer conductance and energy balance require coefficients extracted from the header (see previous function)
# Argument "licor" is a LI-6800 gas exchange dataframe, produce by the previous function
calc_licor6800 = function(licor) {

  x = licor$Fan_speed*licor$Pa / (1000*licor$blc_Po)
  y = pmax(pmin(licor$S, licor$blc_maxS) , licor$blc_minS)

  licor %>%
    mutate(
      gbw = ifelse(Geometry == "0: Broadleaf",
                   # Then apply broadleaf formula.
                   blc_a + blc_b*x + blc_c*x*y*y + blc_d*x*y + blc_e*x*x,
                   # Otherwise, check if needle geometry is specified
                   ifelse(Geometry == "1: Needle",
                          # Then apply needle value.
                          3,
                          # Otherwise, apply custom BLC value
                          Custom
                   )
      ),
      E = Flow * CorrFact * (H2O_s-H2O_r)/(100*S*(1000-CorrFact*H2O_s)),
      A = Flow * CorrFact * (CO2_r-CO2_s*(1000-CorrFact*H2O_r)/(1000-CorrFact*H2O_s))/(100*S),
      Ca = CO2_s - ifelse(CorrFact>1, A*S*100 , 0),
      Rabs = Qin * convert,
      VPcham = H2O_s * (Pa+`ΔPcham`)/1000,
      SVPcham = 0.61365 * exp(17.502*Tair/(240.97+Tair)),
      RHcham = VPcham/SVPcham*100,
      TleafEB = (Tair+(Rabs+2*0.95*0.0000000567*(((Tair+deltaTw)+273)^4 - (Tair+273)^4)-44100*E)/(1.84*29.3*gbw+8*0.95*0.0000000567*(Tair+273)^3)),
      TleafCnd = fT1*Tleaf + fT2*Tleaf2 + fTeb*TleafEB,
      SVPleaf = 0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd)),
      VPDleaf = (SVPleaf-H2O_s*(Pa+`ΔPcham`)/1000),
      LatHFlux = -E*44100,
      SenHFlux = 2*29.3*gbw*0.92*(Tair-TleafCnd),
      NetTherm = 2*0.95*0.0000000567*(((Tair+deltaTw)+273)^4-(TleafCnd+273)^4),
      EBSum = Rabs+NetTherm+LatHFlux+SenHFlux,
      gtw = E*(1000-(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)+H2O_s)/2)/(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)-H2O_s),
      gsw = 2 / ((1/gtw - 1/gbw) + sign(gtw)*sqrt((1/gtw - 1/gbw)^2 + 4*K/((K+1)^2)*(2/(gtw*gbw) - 1/(gbw^2)))),
      gtc = 1 / ((K+1)/(gsw/1.6)+1/(gbw/1.37)) + K/((K+1)/(gsw/1.6) + K/(gbw/1.37)),
      Ci = ((gtc-E/2)*Ca-A)/(gtc+E/2),
      Pci = Ci*(Pa+`ΔPcham`)/1000,
      Pca = (CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0)) * (Pa+`ΔPcham`)/1000
    )
}



# This function applies a variable match correction to gas exchange data acquired using
# the FAsTeR method; see Garen and Michaletz (2024), New Phytologist, for details.
# Argument "data" is a LI-6800 gas exchange dataframe
match_correct = function(data) {

  return_data = c()

  for (i in unique(data$curveID)) {

    # Grab each curve in our dataset and figure out how long it is
    a = subset(data, curveID == i)
    len = dim(a)[1]

    # Interpolate match offset corrections from original match value and test point, assuming match accumulates linearly with time
    correction_co2 = (a$time-a$time[1])*(a$co2_adj[len]-a$co2_adj[len-1])/(a$time[len]-a$time[1])
    correction_h2o = (a$time-a$time[1])*(a$h2o_adj[len]-a$h2o_adj[len-1])/(a$time[len]-a$time[1])

    # Add to data frame
    a$correction_co2 = correction_co2
    a$correction_h2o = correction_h2o

    # Change sample CO2 and H2O values using new correction
    b=a
    b$CO2_s = b$CO2_s + b$correction_co2
    b$H2O_s = b$H2O_s + b$correction_h2o

    # Recalculate gas exchange variables and drop test point
    b = b %>% calc_licor6800()
    b = b[1:len-1,]

    # Bind together
    return_data = bind_rows(return_data, b)

  }

  return(return_data)
}




# Function for computing assimilation rate using the dynamic method; see Garen and
# Michaletz (2024), and Saathof and Welles (2021) for complete details.
# Arguments:
# data is a dataframe with LI-6800 data
# dt1 is the calibration tuning constant from licor for REF
# dt2 is the calibration tuning constant from licor for SAMPLE
# aV is the effective volume tuning constant
noneq_correct_full = function(data, dt1, dt2, aV) {

  return_data = c()

  # Iterate over curves in dataset
  for (i in unique(data$curveID)) {

    # Grab first curve
    non.eq.AT = subset(data, curveID == i)

    len = dim(non.eq.AT)[1]

    co = non.eq.AT$CO2_s/1e6 #umol/mol
    ce = non.eq.AT$CO2_r/1e6 #umol/mol

    wo = non.eq.AT$H2O_s/1e3 #mmol/mol
    we = non.eq.AT$H2O_r/1e3 #mmol/mol

    # 1. Get Crd and Csd at each t
    Crd = ce/(1-we)*1e6
    Csd = co/(1-wo)*1e6

    # 2. Interpolate Crd at -dt1
    # Find the two points which straddle time-dt1 in time
    # If one of them matches exactly, use that one
    # Otherwise, interpolate between values

    time = non.eq.AT$time
    past_time = time-dt1

    Crd_dt1 = rep(NA, len)

    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        Crd_dt1[i] = NA
        next
      }

      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        Crd_dt1[i] = Crd[which(past_time[i] == time)]
        next
      }

      # Otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)

      time_diff = past_time[i]-time[last_obs]
      slope = (Crd[next_obs]-Crd[last_obs])/(time[next_obs]-time[last_obs])
      Crd_last = Crd[last_obs]

      Crd_dt1[i] = Crd_last + slope*time_diff

    }

    # 3. Compute dCsd/dt at each t
    dCsc_dt = (Csd[2:len]-Csd[1:len-1])/(time[2:len]-time[1:len-1])
    dCsc_dt = c(NA, dCsc_dt)

    dCsc_dt_dt2 = rep(NA, len)

    # 4. Interpolate dCsd/dt at -dt2
    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        dCsc_dt_dt2[i] = NA
        next
      }

      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        dCsc_dt_dt2[i] = dCsc_dt[which(past_time[i] == time)]
        next
      }

      # otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)

      time_diff = past_time[i]-time[last_obs]
      slope = (dCsc_dt[next_obs]-dCsc_dt[last_obs])/(time[next_obs]-time[last_obs])
      dCsc_dt_last = dCsc_dt[last_obs]

      dCsc_dt_dt2[i] = dCsc_dt_last + slope*time_diff

    }

    # 5. Finally, compute Adyn

    S = non.eq.AT$S # cm^2
    Tair = non.eq.AT$Tair
    Flow = non.eq.AT$Flow
    Pa = non.eq.AT$Pa
    H2OR = non.eq.AT$H2O_r

    Adyn = (Crd_dt1 - Csd - Pa*1000/(8.314*(Tair+273.15)) * aV/Flow * dCsc_dt_dt2) * Flow/(100*S) * (1000-H2OR)/1000

    non.eq.AT$A = Adyn

    return_data = bind_rows(return_data, non.eq.AT)
  }

  # Drop any with A = NA
  return_data = subset(return_data, !is.na(A))

  return(return_data)

}

