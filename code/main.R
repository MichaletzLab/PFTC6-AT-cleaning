# Download and clean assimilation-temperature response data
# Plant Functional Traits Course 6 - Western Norway
# Last updated, 9 April 2025, Josef Garen

#### 1. Header ####

# Load required libraries
library(osfr)
library(readxl)
library(tidyverse)
library(LeafArea)
library(stringr)
library(progress)
library(rTPC)
library(nls.multstart)

# Load helper functions
source("code/gas_exchange_functions.R")
source("code/quality_control_functions.R")

# Create directory structure
dir.create("data")



#### 2. Download and extract data ####

# Download and extract leaf scans
osf_leaf_scans = osf_retrieve_file("h8wjy")
osf_download(x = osf_leaf_scans, path = "data")
unzip("data/AT_leaf_scans.zip", exdir = "data")
file.remove("data/AT_leaf_scans.zip")

# Download and extract AT files
osf_AT_raw = osf_retrieve_file("ubf37")
osf_download(x = osf_AT_raw, path = "data")
unzip("data/AT_raw_data.zip", exdir = "data/AT_raw_data")
file.remove("data/AT_raw_data.zip")

# Download LI-6800 file metadata
osf_AT_meta = osf_retrieve_file("jqvr7")
osf_download(x = osf_AT_meta, path = "data")

# Clean up
rm(osf_leaf_scans, osf_AT_raw, osf_AT_meta)



#### 3. Load and amend metadata ####

# Load metadata; remove unneeded columns
AT_metadata <- read_excel("data/PFTC6_Norway_Leaf_traits_2022.xlsx", sheet = "Data")
AT_metadata = AT_metadata %>% subset(remark == "cut") %>%
  select(ID, siteID, taxon, curveID = individual_nr)

# Make corrections to metadata:
#Curve 1083 was measured on 30 Jul prior to the collection day at site 4. Change the siteID in sean.data to site 2 to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1083', "siteID"] <- "Hog"
#For curve 1104 both sites 2 and 4 were measured on 31 Jul. AT data is more reliable, especially as 'rest' envelope was labeled Lia. Change the siteID in sean.data to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1104', "siteID"] <- "Lia"
#For curve 1106 site 3 was not measured on 31 Jul. AT data is more reliable, especially as 'rest' envelope was labeled Lia. Change the siteID in sean.data to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1106', "siteID"] <- "Lia"
#For curve 1113 site 3 was not measured on 31 Jul. AT data is more reliable, especially as 'rest' envelope was labeled Lia. Change the siteID in sean.data to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1113', "siteID"] <- "Lia"
#For curve 1116 there was a species entry error. The scans Ecx0771 and DYJ1799 are both Alchemilla. Change taxon to reflect that
AT_metadata[AT_metadata$curveID=='1116', "taxon"] <- "Alchemilla alpina"
#For curve 1118 both sites 2 and 4 were measured on 31 Jul. AT data is more reliable, especially as 'cut' envelope was labeled Hog. Change the siteID in sean.data to site 2 to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1118', "siteID"] <- "Hog"
#For curve 1130 there was a species entry error. The scans HTW0662 and HUA5840 are both Alchemilla. Change HTw0662 to reflect that
AT_metadata[AT_metadata$curveID=='1130', "taxon"] <- "Alchemilla alpina"
#For curve 1152 AT data is more reliable, especially as 'cut' envelope was labeled Joa. Change the siteID in sean.data to reflect that and match the AT data
AT_metadata[AT_metadata$curveID=='1152', "siteID"] <- "Joa"



#### 4. Process leaf scans for area ####

# Find install location of imageJ
ij.location = find.ij()[1]
if(ij.location == "ImageJ not found") {
  stop("ImageJ not found. Please install ImageJ in the common install location or specify the install location of ImageJ")
}

# Run imageJ on leaf scans
leaf.area.1 = run.ij(path.imagej = ij.location,
                     set.directory = file.path(getwd(), "data", "AT_leaf_scans", "1"),
                     distance.pixel = 2480,
                     trim.pixel = 300,
                     low.size = 0.01,
                     prefix = "\\.")
leaf.area.2 = run.ij(path.imagej = ij.location,
                     set.directory = file.path(getwd(), "data", "AT_leaf_scans", "2"),
                     distance.pixel = 2480,
                     trim.pixel = 300,
                     low.size = 0.01,
                     prefix = "\\.")
leaf.area.3 = run.ij(path.imagej = ij.location,
                     set.directory = file.path(getwd(), "data", "AT_leaf_scans", "3"),
                     distance.pixel = 2480,
                     trim.pixel = 300,
                     low.size = 0.01,
                     prefix = "\\.")
leaf.area.4 = run.ij(path.imagej = ij.location,
                     set.directory = file.path(getwd(), "data", "AT_leaf_scans", "4"),
                     distance.pixel = 2480,
                     trim.pixel = 300,
                     low.size = 0.01,
                     prefix = "\\.")
leaf.area.5 = run.ij(path.imagej = ij.location,
                     set.directory = file.path(getwd(), "data", "AT_leaf_scans", "5"),
                     distance.pixel = 2480,
                     trim.pixel = 300,
                     low.size = 0.01,
                     prefix = "\\.")

# Bind all leaf areas together
leaf.area.all = bind_rows(leaf.area.1, leaf.area.2, leaf.area.3, leaf.area.4, leaf.area.5)

# Clean up
rm(leaf.area.1, leaf.area.2, leaf.area.3, leaf.area.4, leaf.area.5, ij.location)



#### 5. Load LI-6800 gas exchange data ####
# Load list of file paths for LI-6800 data
at.files = list.files("data/AT_raw_data/", full.names=T)

# Extract list of IDs from the file names
at.ids <- lapply(at.files, str_sub, -4, -1) # save the individual IDs, extracted from file names

# Load LI-6800 data from all files
file.x = list()
for(i in 1:length(at.files)){
  x = at.files[[i]] # take a single file name

  if(file.info(x)$size >2000){ # Remove files that are empty and would throw an error
    y <- read_6800_with_BLC(x) # Read raw LI-6800 file
    at.id <- at.ids[[i]]       # Take corresponding individual ID

    y$curveID <- rep(at.id, nrow(y)) # Add ID as a column to file
    y$obs <- as.integer(y$obs)
    y$time <- as.double(y$time)

    file.x[[i]] <- y # Add file (now a dataframe) to the list of files/dataframes
  }
}

file.x <- file.x[lengths(file.x) != 0]    # Remove empty placeholders
at.all <- do.call(bind_rows, file.x)       # Combine files into one dataframe
at.all <- subset(at.all, A >= -5)           # Remove any unreasonable A values
at.all$curveID = as.numeric(at.all$curveID) # Change data type

# Clean up
rm(file.x, at.ids, at.files, x, y, at.id, i)



#### 6. Merge gas exchange, leaf area, and metadata ####

# Join gas exchange data to metadata
at.all <- left_join(at.all, AT_metadata, by="curveID")

# Join gas exchange data to leaf area
leaf.area.all = leaf.area.all %>% rename(sampleID = sample)
at.all = at.all %>% rename(sampleID = ID.y)
at.all = left_join(at.all, leaf.area.all, by = "sampleID")

# Replace S (LI-6800 variable for area) with correct leaf area
at.all = at.all %>% mutate(S = total.leaf.area)

# Set any NAs to 2 (default small aperture size)
at.all[which(is.na(at.all$S)),]$S = 2

# Clean up
rm(leaf.area.all, AT_metadata)



#### 7. Perform post-measurement corrections ####
at.corr = match_correct(at.all) # Apply match offset correction

# Apply non-equilibrium correction:
# Here we  split this up based on flow rate; they are all either 250 or 300
at.corr.250 = subset(at.corr, Flow < 275)
at.corr.300 = subset(at.corr, Flow >= 275)

# Apply noneq calculation; Calibration coefficients were obtained using the Saathof and Welles method
at.corr.250 = noneq_correct_full(at.corr.250, dt1 = 5.53, dt2 = 3.52, aV = 66.44)
at.corr.300 = noneq_correct_full(at.corr.300, dt1 = 4.32, dt2 = 2.77, aV = 66.46)

# Rejoin datasets
at.all = bind_rows(at.corr.250, at.corr.300)

# Clean up
rm(at.corr, at.corr.250, at.corr.300)



#### 8. Quality checks on AT data ####
# Perform quality checks and remove low quality curves; this step is time consuming
at.clean = cut.multimodal(at.all)         # Remove curves with multiple "optima"
at.clean = cut.topt.out.of.bounds(at.clean) # Remove curves without a well-defined T_opt

# Clean up
rm(at.all)


#### 9. Output cleaned data ####
write.csv(at.clean, file = "data/AT_clean.csv", row.names = F)



#### END OF SCRIPT ####
