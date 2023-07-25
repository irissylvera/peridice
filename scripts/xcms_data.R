library(xcms)
library(tidyverse)

# generate raw data list - CHANGE FOR YOUR FILEPATH
raw_data <- list.files('HILICpos/', full.names = TRUE)
raw_data <- as.data.frame(raw_data) %>% 
  filter(!str_detect(raw_data, "DDA"))

# extract metadata from filenames
metadataframe <- raw_data %>%
  # Grab just the unique filenames and rename the column
  distinct(filename=`raw_data`) %>%
  # Create a new column with sample type information
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  # Create a new column with timepoint information
  mutate(timepoint=str_extract(filename, "27June|30June|14July|21July|26July|27July")) %>%
  # Extract treatment information from the filename
  # Longest one needs to go first!
  mutate(treatment=str_extract(filename, "Std|Poo|Blk|-C|ZF|ZL|ZH|RL|RH|LL|LH")) %>%
  mutate(treatment = str_remove(treatment, "-")) %>% 
  # Replace accidental "P" treaments from "Pooled" with NAs
  mutate(colid = treatment)

# xcms file generation
# msnexp MS1
register(BPPARAM = SerialParam(progressbar = TRUE))
msnexp <- readMSData(
  files = metadataframe$filename,
  pdata = new("NAnnotatedDataFrame", metadataframe),
  msLevel. = 1,
  mode = "onDisk"
)

# save so don't have to generate again
saveRDS(msnexp, file = "msnexp.rds")

# read in saved RDS object - don't need to run unless environment cleared
# msnexp <- readRDS("msnexp.rds")

# set prefilter QC
prefilter_versioned <- c(3, 1e6)

# chromatogram peak picking
register(BPPARAM = SnowParam(workers = 3, tasks = nrow(metadataframe), progressbar = TRUE))
cwp <- CentWaveParam(
  ppm = 5,
  peakwidth = c(20, 80),
  prefilter = prefilter_versioned,
  snthresh = 0,
  verboseColumns = TRUE,
  extendLengthMSW = TRUE,
  integrate = 2
)
# generate msnexp object with peak detection
msnexp_withpeaks <- findChromPeaks(msnexp, cwp)

# set retention time correction
register(BPPARAM = SerialParam(progressbar = TRUE))
obp <- ObiwarpParam(
  binSize = 0.1,
  centerSample = round(nrow(metadataframe)/2),
  response = 1,
  distFun = "cor_opt"
)

# generate msnexp object with retention time correction
msnexp_rtcor <- adjustRtime(msnexp_withpeaks, obp)

# generate msnexp object with grouping of metabolites
pdp <- PeakDensityParam(
  sampleGroups = metadataframe$colid,
  bw = 12,
  minFraction = 0.1,
  binSize = 0.001,
  minSamples = 2
)
msnexp_grouped <- groupChromPeaks(msnexp_rtcor, pdp)

# generate msnexp object with filled peaks
fpp <- FillChromPeaksParam(ppm = 2.5)
msnexp_filled <- fillChromPeaks(msnexp_grouped, fpp)

# save RDS files
saveRDS(msnexp_withpeaks, file = "msnexp_withpeaks.rds")
msnexp <- readRDS("msnexp.rds")
saveRDS(msnexp_rtcor, file = "msnexp_rtcor.rds")
saveRDS(msnexp_grouped, file = "msnexp_grouped.rds")
saveRDS(msnexp_filled, file = "msnexp_filled.rds")