## 00_dataClean : Script for cleaning/organizing the spec data
# Ceili DeMarais

library(spectrolab)
library(dplyr)
library(tidyr)
library(signal)

## set where to find the stuff
data_base        <- "data_in/spec_data" 
data_out_raw     <- "data_work/spec_clean.csv"
data_out_proc    <- "data_work/spec_clean_processed.csv"
spectra_per_leaf <- 5
leaf_label       <- c("bot", "mid", "top")

## func to label leaf pos
get_leaf_pos <- function(spec_index, spectra_per_leaf, leaf_label) {
  leaf_num <- floor(spec_index / spectra_per_leaf) + 1
  if (leaf_num > length(leaf_label)) return(NA_character_)
  return(leaf_label[leaf_num])
}

## preprocessing function (SNV + detrend + SG smooth)
preprocess_spectra <- function(raw_spectra) {
  snv <- function(x) (x - mean(x)) / sd(x) #bright norn
  snv_spectra <- t(apply(raw_spectra, 1, snv))
  
  detrend_spectrum <- function(x) { #detrend
    lm_fit <- lm(x ~ seq_along(x))
    residuals(lm_fit)
  }
  detrended_spectra <- t(apply(snv_spectra, 1, detrend_spectrum))
  
  smoothed_spectra <- t(apply(detrended_spectra, 1, function(x) { #SG smooooth
    sgolayfilt(x, p = 2, n = 11, m = 0)
  }))
  
  processed_df <- as.data.frame(smoothed_spectra)
  colnames(processed_df) <- colnames(raw_spectra)
  return(processed_df)
}

## point to the time pt folders
tp_folders <- list.dirs(data_base, recursive = FALSE, full.names = TRUE)
tp_folders <- tp_folders[grepl("data_\\d+$", tp_folders)]

target_wvl <- seq(400, 2400, by = 1)
all_data   <- list()

for (tp_path in tp_folders) {
  tp_name <- basename(tp_path)
  
  spectra_raw <- read_spectra(tp_path, format = "asd")
  
  if (inherits(spectra_raw, "spectra")) {
    spectra <- spectrolab::resample(spectra_raw, new_bands = target_wvl, fwhm = 1)
  } else {
    spectra <- combine(lapply(spectra_raw, resample, new_bands = target_wvl, fwhm = 1))
  }
  
  dati <- as.data.frame(spectra) %>%
    mutate(
      raw_name  = sample_name,
      treatment = sub("^([A-Za-z]+).*", "\\1", sample_name),
      plantID   = sub("^[A-Za-z]+(\\d+)_.*", "\\1", sample_name),
      specNum   = as.integer(sub(".*_(\\d+)$", "\\1", sample_name)),
      leafPos   = sapply(specNum, get_leaf_pos,
                         spectra_per_leaf = spectra_per_leaf,
                         leaf_label = leaf_label),
      timePt    = tp_name
    )
  
  all_data[[tp_name]] <- dati
}

combined <- bind_rows(all_data)

## average 5 spectra per leaf of thee raw
averaged_raw <- combined %>%
  dplyr::filter(!is.na(leafPos)) %>%
  group_by(timePt, treatment, plantID, leafPos) %>%
  summarise(
    n_spectra = n(),
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )

## write raw
write.csv(averaged_raw, data_out_raw, row.names = FALSE)

## process averaged spectra
meta_cols <- c("timePt", "treatment", "plantID", "leafPos", "n_spectra")
wvl_cols  <- names(averaged_raw)[grepl("^[0-9]", names(averaged_raw))]

raw_matrix    <- as.matrix(averaged_raw[, wvl_cols])
processed_mat <- preprocess_spectra(raw_matrix)

averaged_proc <- cbind(
  averaged_raw[, meta_cols],
  processed_mat
)

## write processed
write.csv(averaged_proc, data_out_proc, row.names = FALSE)