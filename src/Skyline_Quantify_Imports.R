# Imports for the quantification step

Ingalls.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                              stringsAsFactors = FALSE, header = TRUE) %>%
  # filter(Column == Column.Type) %>%
  rename(Metabolite.Name = Compound_Name) %>%
  select(Metabolite.Name, Compound_Type, Concentration_uM, HILIC_Mix, Empirical_Formula) %>% #QE.RF.ratio
  filter(!is.na(Concentration_uM)) 
Ingalls.Standards$Metabolite.Name <- TrimWhitespace(Ingalls.Standards$Metabolite.Name)


# Import BMIS'd sample file ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))
BMISd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

# Import QC'd files and remove parameter data ------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = QC.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value, Retention.Time, SN.Flag, contains("Flag"))) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area.with.QC, Area, everything())

# Import Internal standards key ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_extras/", pattern = names.pattern))
filepath <- file.path("data_extras", paste(filename, ".csv", sep = ""))

IS.key <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  rename(FinalBMIS = Internal_Standards)

# Rename Ingalls standards to join with QC data
Ingalls.Standards.Join <- Ingalls.Standards %>% 
  rename(Precursor.Ion.Name = Metabolite.Name)

# Apply appropriate filters and isolate standards ---------------------------------------------------------------
Full.stds.data <- QCd.data %>%
  filter(Precursor.Ion.Name %in% Ingalls.Standards.Join$Precursor.Ion.Name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls.Standards.Join, by = "Precursor.Ion.Name") %>%
  select(Replicate.Name, Precursor.Ion.Name, Compound_Type, everything()) %>%
  unique()