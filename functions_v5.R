extract_abundances <- function(data,
                              normalized = 0) {
  print.noquote("Extracting data")

  # Separating the columns to extract from PQI data set
  iNorm <- which(colnames(data) == "Normalised.abundance")
  iRaw <- which(colnames(data) == "Raw.abundance")
  N <- as.numeric(iRaw - iNorm)
  fNorm <- iNorm + (N - 1)
  fRaw <- iRaw + (N - 1)

  if (normalized == 0) {
    subset_range <- as.vector(iRaw:fRaw)
    print.noquote("Raw data extracted")
  }
  else {
    subset_range <- as.vector(iNorm:fNorm)
    print.noquote("Normalized data extracted")
  }

  print.noquote("Restructuring")
  subset_data <- as.data.frame(data[, c(1:3, 5, subset_range)])
  
  print.noquote("Done")
  
  return(subset_data)
}
  

fill_groups <- function(data) {
  df <- as.data.frame(data)

  # Setting group's coordinates
  groups <- list(
  name = as.character(df[1, which(df[1, ] != "")]),
  from = as.numeric(c(which(df[1, ] != ""))),
  to = as.numeric(c((which(df[1, ] != "")[-1] - 1), ncol(df)))  
)
  print.noquote("Group Coordinates [OK]")

  # Filling group names on PQI data set
  for (g in 1:length(groups[["name"]])) {
    df[1, (groups$from[g]:groups$to[g])] <- groups$name[g]
  }
  print.noquote("Filling group names on PQI data set [OK]")
  
  print.noquote("Done")
  
  return(df)
}

extract_samples <- function(data) {
  # Isolating the "Samples" Data
  samples <- t(as.data.frame(data[, -c(2:4)]))
  samples[1, c(1, 2)] <- c("Groups", "Samples")
  colnames(samples) <- as.character(samples[1, ])
  rownames(samples) <- NULL
  samples <- as.data.frame(samples[-1, ])
  samples <- as.data.frame(samples[,c(2, 1, 3:ncol(samples))])

  samples[, 1] <- as.character(samples[, 1])
  samples[, 2] <- as.factor(samples[, 2])
  for (i in 3:ncol(samples)) {
    samples[, i] <- as.double(samples[, i])
  }
  print.noquote("Done")

  return(as.data.frame(samples))
}

rm_ghost_peaks <- function(samples_data, 
                             method = 1, # 0 = mean // 1 = median
                             referenceGroups = c("Blank")) { # 0 = mean // 1 = median
                              
  # Creates a `features` reference list
    featuresList <- as.vector(colnames(samples_data)[-c(1,2)])

  # Splits the data set into a list of data sets fragmented by groups.
    splitSamples <- split.data.frame(samples_data[,-c(1, 2)],
                                f = samples_data[, "Groups"]) # NEW FEATURE ? ref_group = "Groups"
    
    highestMedianGroup <- NULL
    highestMedianGroup <- vector(mode = "character", length(featuresList))
    names(highestMedianGroup) <- as.vector(featuresList)
  
  # Creating a temporary placeholder vector for median values and its group
  #   it is rather **unitary** than cyclic, thus, it creates dependency on the
  #   _features_, `f` via a loop. It mainly uses vectorized operations though.
    for (f in 1:length(featuresList)) {
      tmpMedianVector <- NULL
      tmpMedianVector <- unlist(
                           lapply(
                             sapply(splitSamples, 
                                    function(x) x[[as.character(featuresList[f])]]),
                             ifelse(method == 1, median, mean) # Everything works fine here, except that, for some reason, when the median is calculated, it returns 'Warning: number of items to replace is not a multiple of replacement length'. I couldn't dig into it right now. 
                             )
                           )
      
      # Assigning the Group which each respective variable has the highest median
      highestMedianGroup[
        which(
          names(highestMedianGroup) == as.character(featuresList[f])
          )
        ] <- as.character(
        names(
          which(
            tmpMedianVector == max(tmpMedianVector)
            )
          )
        )
    }
  
  print.noquote("Peaks detected")
  
  
    return(cbind.data.frame(
      samples_data[,c(1,2)], 
      samples_data[, colnames(samples_data) %in% names(which(highestMedianGroup != referenceGroups))]))
  
  print.noquote("Peaks cleaned")  
}


rm_missing_values <- function(samples_data,
                            featureValue = 0,
                            mvThreshold = 0.9) {

  toExclude <- vector(mode = "numeric", length = length(as.vector(colnames(samples_data)[-c(1,2)])))
  names(toExclude) <- as.vector(colnames(samples_data)[-c(1,2)])

  correctionFactor <- 0.5 # round up > .5 and down otherwise
  
  samplesThreshold <- floor(
    (nrow(samples_data) * mvThreshold) + correctionFactor
    )
  
  for (f in 1:length(toExclude)) {
    
  nMissing <- length(which(samples_data[, f] == featureValue))
  
  toExclude[f] <- ifelse(nMissing >= samplesThreshold, 1, 0)
  }
  
  return(
      cbind.data.frame(
        samples_data[, c(1,2)],
        samples_data[, colnames(samples_data) %in% names(toExclude[toExclude == 0])]
        )
    )
  
}

rsd_calc <- function(input_data,
                threshold_rsd = 0.3,
                reference_group = "QC") {
  
  # RSD threshold (0.15 - 0.30)
  
  subset_rsd <- as.data.frame(input_data[which(input_data$Groups == reference_group), ]) # Sub-setting necessary data
  
  rsd_stdev <- apply(subset_rsd[, -c(1, 2)], 2, sd)
  rsd_mean <- apply(subset_rsd[, -c(1, 2)], 2, mean)
  rsd_results <- rsd_stdev / rsd_mean # Calculating the actual RSD
  
  rsd_true <- ifelse(rsd_results < threshold_rsd, TRUE, FALSE)
  
  rsd_in <- input_data[, names(input_data) %in% names(rsd_true[rsd_true])]
  
  output_data <- cbind.data.frame(input_data[, c(1, 2)], rsd_in)
  
  return(output_data)
}



