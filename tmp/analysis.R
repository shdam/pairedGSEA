# nice /home/ctools/opt/R-4.0.5/bin/R
library(rhdf5)

### Define samples of interest
samples <- c("GSM1383738","GSM1383739","GSM1383740","GSM1383741","GSM1383742")
### Define file to read from
archsdb <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"

txCount <- loadArchs4(samples, archsdb)
