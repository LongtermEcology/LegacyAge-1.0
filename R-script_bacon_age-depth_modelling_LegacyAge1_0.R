# Code for LegacyAge 1.0
# Alexander Karl Postl
# Alfred-Wegener-Institute Potsdam, Germany 2021

#-----Manual---------------------------------------------

# Step 1: Rows 24-29 -> Install the packages you do not have
# Step 2: Row 32 -> Insert <our folder for the output data
# Step 3: Row 35 -> Insert the Dataset_IDs of your interest

#-----Resultfolders-----------------------------------

# Ages.txt      ->  Chronology tables by Bacon
# Bacon.pdf     ->  Outputplot by Bacon
# Calibration   ->  Plots from Calibration
# ID.Subsets    ->  Summarized data of the ID
# Plot.png      ->  Plot to compare with other chronologies
# Plot.flipped  ->  the same plot but flipped
# Sites         ->  all data concerning the ID

#-----Begin Input-------------------------------------

# R packages that are needed. If you do not have one of the packages just delete the "#" one time
# install.packages("tidyverse")
# install.packages("stringr")
# install.packages("ggpubr")
# install.packages("rlist")
# install.packages("rbacon")
# install.packages("IntCal")

# Set the path to a folder that will contain the output
folder <- "../LegacyAge" 

# IDs that are requested | the the commented-out second row would select all IDs
parameter.ID <- c(47502,4320)
# parameter.ID<-parameter$Dataset_ID

# Download tables that are needed from pangaea
metadata <- read.csv2(
  "https://download.pangaea.de/reference/111158/attachments/Table-S1_chronological_control_points_metadata.csv",
  stringsAsFactors = FALSE,
  sep = "\t", dec = ".")

parameter <- read.csv2(
  "https://download.pangaea.de/reference/111160/attachments/Table-S3_bacon_parameter_settings.csv",
  stringsAsFactors = FALSE,
  sep = "\t", dec = ".")

AgeDepthPollen <- read.csv2(
  "https://download.pangaea.de/reference/111161/attachments/Table-S4_original_chronology_metadata_by_pollen_records.csv",
  stringsAsFactors = FALSE, sep = "\t", dec = ".")

#-----End Input | You do not need to insert anything beyond this line----------------------

# Creating Directories
dir.create(paste0(folder, "/Sites"))
dir.create(paste0(folder, "/Bacon.pdf"))
dir.create(paste0(folder, "/Ages.txt"))
dir.create(paste0(folder, "/ID.Subsets"))
dir.create(paste0(folder, "/Plot.png"))
dir.create(paste0(folder, "/Plot.flipped"))
dir.create(paste0(folder, "/Calibration"))

# Libraries
library(tidyverse)
library(stringr)
library(ggpubr)
library(rlist)
library(rbacon)
library(IntCal)

# Additional functions
`%notin%` <- Negate(`%in%`)
numextract <- function(string) {
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

# Define specific Values
delete.chrono <- c("Radiocarbon years BP", "No chronology", NA)
c14 <- c("Carbon-14", "14C", "C14", "14c", "c14", "", "Radiocarbon years BP")
plot.c14.correct <- "age.dating"
plot.AWI <- "AWI.age"

#-----Loop for Every ID-----------

# loop run for every ID
for (ID in parameter.ID) {

  # creating specific subsets for the ID
  parameter.subset <- parameter[which(parameter[, 1] == ID), ]
  AgeDepthPollen.subset <- AgeDepthPollen[which(AgeDepthPollen$Dataset_ID == ID), ]
  metadata.subset <- metadata[which(metadata$Dataset_ID == ID), ]
  if (is.na(parameter.subset$Reservoir)) {
    parameter.subset$Reservoir <- 0
  }
  if (is.na(parameter.subset$Waterline)) {
    parameter.subset$Waterline <- 0
  }
  calibration <- data.frame(as.numeric(metadata.subset$Depth..cm.), as.numeric(metadata.subset$Age_Uncalibrated..kyr.BP.) * 1000, (as.numeric(metadata.subset$Dating.Error_Older..kyr.) + as.numeric(metadata.subset$Dating.Error_Younger..kyr.)) * 500, metadata.subset$Dating_Method, as.numeric(metadata.subset$Age_Calibrated..kyr.BP.) * 1000, as.numeric(metadata.subset$Calibrated.dating_Error..kyr.) * 1000, stringsAsFactors = FALSE)
  names(calibration) <- c("depth", "age", "e.older", "age.type", "cal.age", "cal.age.se")
  cal.age <- cal.age.se <- cal.sp <- NULL

  # Choosing Calibration Curve

  if (AgeDepthPollen.subset$Latitude..DD.[1] >= 0) {
    cc <- 1
    cc1 <- min(IntCal::ccurve(1)[, 2] - IntCal::ccurve(1)[, 3])
    cc2 <- max(IntCal::ccurve(1)[, 2] + IntCal::ccurve(1)[, 3])
  }
  if (AgeDepthPollen.subset$Latitude..DD.[1] < 0) {
    cc <- 3
    cc1 <- min(IntCal::ccurve(3)[, 2] - IntCal::ccurve(3)[, 3])
    cc2 <- max(IntCal::ccurve(3)[, 2] + IntCal::ccurve(3)[, 3])
  }
  if (parameter.subset$Marine) {
    cc <- 2
    cc1 <- min(IntCal::ccurve(2)[, 2] - IntCal::ccurve(2)[, 3])
    cc2 <- max(IntCal::ccurve(2)[, 2] + IntCal::ccurve(2)[, 3])
  }


  # Defining the postbomb

  if (AgeDepthPollen.subset$Latitude..DD.[1] > 40) {
    postbomb <- 1
  }
  if (AgeDepthPollen.subset$Latitude..DD.[1] < 40) {
    postbomb <- 2
  }
  if (AgeDepthPollen.subset$Latitude..DD.[1] < 0) {
    postbomb <- 4
  }

  # Sperating uncalibrated C14 from the other datings
  if (is.na(parameter.subset$Waterline)) {
    parameter.subset$Waterline <- 0
  }
  for (j in (1:nrow(calibration))) {
    if (is.na(calibration[j, 3])) {
      if (calibration[j, 2] < 11460) {
        calibration[j, 3] <- 30
      } else {
        calibration[j, 3] <- 60
      }
    }
    if (calibration[j, 4] %in% c14 & (calibration[j, 2] + calibration[j, 3] - parameter.subset$Reservoir) > cc1 & (calibration[j, 2] - calibration[j, 3] - parameter.subset$Reservoir) < cc2) {
      calibration[j, 4] <- cc
    } else {
      calibration[j, 4] <- 0
    }
  }
  calibration.0 <- calibration[which(calibration$age.type == 0), ]
  calibration.0[, (5:6)] <- calibration.0[, (2:3)]
  calibration.1 <- calibration[which(calibration$age.type == cc), ]

  # Are there non calibrated C14 datings?
  if (dim(calibration.1)[1] >= 1) {

    # Calibration (Intcal)
    for (Intcal in 1:dim(calibration.1)[1]) {
      calibration.calc <- calibrate(
        age = calibration.1$age[Intcal],
        error = calibration.1$e.older[Intcal],
        cc = as.numeric(calibration.1$age.type[Intcal]),
        postbomb = FALSE,
        reservoir = parameter.subset$Reservoir,
        graph = TRUE,
        title = paste0(ID, "_", calibration.1$depth[Intcal]),
        mar = c(3.5, 3, 2, 1),
        mgp = c(1.7, 0.8, 0),
        bty = "n"
      )
      dev.print(pdf, paste0(folder, "/Calibration/", ID, "_", calibration.1$depth[Intcal], ".pdf"))

      # Standardization of the deviation & weighted mean age
      standart <- (1 / sum(calibration.calc[[1]][, 2])) * calibration.calc[[1]][, 2]
      new.memory <- sum(standart * calibration.calc[[1]][, 1])
      new.se <- sum(abs(standart * (calibration.calc[[1]][, 1] - new.memory)))
      spread <- data.frame(calibration.1$depth[Intcal], calibration.calc[[2]][1], calibration.calc[[2]][2], stringsAsFactors = FALSE)

      # setting values (weighted mean age, weighted standart deviation)
      cal.age <- c(cal.age, new.memory)
      cal.age.se <- c(cal.age.se, new.se)
      cal.sp <- rbind(cal.sp, spread)
    }

    # Finalize calibrated data
    calibration.1$cal.age <- cal.age
    calibration.1$cal.age.se <- cal.age.se
  }

  # Create calibration table
  calibration <- rbind(calibration.0, calibration.1)[order(calibration[, 1]), ]

  #-----Bacon modelling-----

  # Create csv for bacon
  data <- data.frame(sample_id = ID, age = calibration[, 2], error = calibration[, 3], depth = calibration[, 1], cc = calibration$age.type, stringsAsFactors = FALSE)
  zero <- data.frame(ID, max(-70, parameter.subset$Starting.age, na.rm = TRUE), 30, parameter.subset$Waterline, 0, stringsAsFactors = FALSE)
  colnames(zero) <- colnames(data)

  # Is in parameter csv a surface added or not
  if (parameter.subset$Add.surface) {
    data <- rbind(zero, data)
  }

  # Write csv for Bacon
  dir.create(paste0(folder, "/Sites/", ID))
  write.table(data[order(data$depth), ], file = paste0(folder, "/Sites/", ID, "/", ID, ".csv"), row.names = FALSE, sep = ",", dec = ".")

  # Calculate accumulation rate
  if (is.na(parameter.subset$Acc.mean)) {
    acc.rate <- (max(data[, 2]) - min(data[, 2])) / (max(data$depth) - min(data$depth))
  } else {
    acc.rate <- parameter.subset$Acc.mean
  }

  # Age-Depth-Modelling (Bacon)
  try(Bacon(as.character(ID),
    ask = FALSE,
    suggest = FALSE,
    thick = if (is.na(parameter.subset$Resolution.cm)) {
      (max(AgeDepthPollen.subset$Depth_Pollen, calibration[, 1]) - parameter.subset$Waterline) / parameter.subset$Resolution.section
    } else {
      parameter.subset$Resolution.cm
    },
    acc.shape = 1.5,
    acc.mean = acc.rate,
    mem.strength = 20,
    mem.mean = parameter.subset$Memory,
    d.min = min(parameter.subset$Waterline, AgeDepthPollen.subset$Depth_Pollen),
    d.max = max(AgeDepthPollen.subset$Depth_Pollen, calibration[, 1]),
    MinYr = parameter.subset$Starting.age,
    postbomb = postbomb,
    hiatus.depths = unique(c(as.numeric(parameter.subset[1, which(names(parameter.subset) == "Hiatus.1"):(which(names(parameter) == "Add.depth.1") - 1)]))),
    coredir = paste0(folder, "/Sites/"),
    close.connections = TRUE,
    delta.R = parameter.subset$Reservoir
  ), silent = TRUE)

  #-----Bacon review-----

  # if bacon was not able to create a model
  if (file.exists(paste0(folder, "/Sites/", ID, "/", list.files(path = paste0(folder, "/Sites/", ID), ".pdf"))) == FALSE) {
    parameter.ID <- parameter.ID[-1]
  } else {

    # Copy model files to new folders
    file.copy(from = paste0(folder, "/Sites/", ID, "/", list.files(path = paste0(folder, "/Sites/", ID), ".pdf")), to = paste0(folder, "/Bacon.pdf/", ID, ".pdf"), overwrite = TRUE)
    file.copy(from = paste0(folder, "/Sites/", ID, "/", list.files(path = paste0(folder, "/Sites/", ID), "ages.txt")), to = paste0(folder, "/Ages.txt/", ID, ".txt"), overwrite = TRUE)

    #-----Age Allocation-----

    # Define basic values for the age allocation
    model.AWI <- read.table(paste0(folder, "/Ages.txt/", ID, ".txt"), header = TRUE)
    model.rest <- round(model.AWI$depth[length(model.AWI$depth)] - floor(model.AWI$depth[length(model.AWI$depth)]), 4)
    new.memory <- data.frame(min = NA, max = NA, median = NA, mean = NA, stringsAsFactors = FALSE)

    # loop for every pollen entry
    for (g in 1:length(AgeDepthPollen.subset$Depth_Pollen)) {

      # Define basic values for the specific pollen entry
      pollen.num <- as.numeric(as.character(floor(AgeDepthPollen.subset$Depth_Pollen[g] - model.rest) + model.rest))
      pollen.rest <- AgeDepthPollen.subset$Depth_Pollen[g] - pollen.num

      # Case I
      if (AgeDepthPollen.subset$Depth_Pollen[g] < model.AWI$depth[1]) {
        lower.age <- model.AWI[1, c(2:5)]
        upper.age <- model.AWI[2, c(2:5)]
        m <- (upper.age - lower.age) / (model.AWI$depth[2] - model.AWI$depth[1])
        b <- upper.age - m * model.AWI$depth[2]
        new.memory[g, ] <- m * AgeDepthPollen.subset$Depth_Pollen[g] + b
      }

      # Case II
      if (AgeDepthPollen.subset$Depth_Pollen[g] >= model.AWI$depth[1] & AgeDepthPollen.subset$Depth_Pollen[g] < model.AWI$depth[length(model.AWI$depth)]) {
        lower.age <- model.AWI[model.AWI$depth == pollen.num, c(2:5)]
        upper.age <- model.AWI[model.AWI$depth == pollen.num + 1, c(2:5)]
        new.memory[g, ] <- lower.age * (1 - pollen.rest) + upper.age * pollen.rest
      }

      # Case III
      if (AgeDepthPollen.subset$Depth_Pollen[g] == model.AWI$depth[length(model.AWI$depth)]) {
        new.memory[g, ] <- model.AWI[length(model.AWI$depth), c(2:5)]
      }

      # Case IV
      if (AgeDepthPollen.subset$Depth_Pollen[g] > model.AWI$depth[length(model.AWI$depth)]) {
        lower.age <- model.AWI[length(model.AWI$depth) - 1, c(2:5)]
        upper.age <- model.AWI[length(model.AWI$depth), c(2:5)]
        m <- (upper.age - lower.age) / (model.AWI$depth[length(model.AWI$depth)] - model.AWI$depth[length(model.AWI$depth) - 1])
        b <- upper.age - m * model.AWI$depth[length(model.AWI$depth)]
        new.memory[g, ] <- m * AgeDepthPollen.subset$Depth_Pollen[g] + b
      }
    }

    #-----Preparing AgeDepthPollen table-----

    # create & write table
    AgeDepthPollen.subset <- cbind(AgeDepthPollen.subset[, 1:13], new.memory, AgeDepthPollen.subset[, 14:37])
    write.table(AgeDepthPollen.subset, paste0(folder, "/ID.Subsets/", ID, ".csv"), row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")


    #-----Plotting-----

    # Define Chronologies to plot (only calibrated chronologies)
    if (any(AgeDepthPollen.subset$Original_Age.type_1 %notin% delete.chrono)) {
      chronology1 <- data.frame(AgeDepthPollen.subset$Original_chronology_1, as.numeric(AgeDepthPollen.subset$Estimated_Age_1) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_1 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology1 <- NULL
    }
    if (any(AgeDepthPollen.subset$Original_Age.type_2 %notin% delete.chrono)) {
      chronology2 <- data.frame(AgeDepthPollen.subset$Original_chronology_2, as.numeric(AgeDepthPollen.subset$Estimated_Age_2) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_2 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology2 <- NULL
    }
    if (any(AgeDepthPollen.subset$Original_Age.type_3 %notin% delete.chrono)) {
      chronology3 <- data.frame(AgeDepthPollen.subset$Original_chronology_3, as.numeric(AgeDepthPollen.subset$Estimated_Age_3) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_3 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology3 <- NULL
    }
    if (any(AgeDepthPollen.subset$Original_Age.type_4 %notin% delete.chrono)) {
      chronology4 <- data.frame(AgeDepthPollen.subset$Original_chronology_4, as.numeric(AgeDepthPollen.subset$Estimated_Age_4) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_4 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology4 <- NULL
    }
    if (any(AgeDepthPollen.subset$Original_Age.type_5 %notin% delete.chrono)) {
      chronology5 <- data.frame(AgeDepthPollen.subset$Original_chronology_5, as.numeric(AgeDepthPollen.subset$Estimated_Age_5) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_5 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology5 <- NULL
    }
    if (any(AgeDepthPollen.subset$Original_Age.type_6 %notin% delete.chrono)) {
      chronology6 <- data.frame(AgeDepthPollen.subset$Original_chronology_6, as.numeric(AgeDepthPollen.subset$Estimated_Age_6) * 1000, as.numeric(AgeDepthPollen.subset$Estimated.age_Error_6 * 1000), stringsAsFactors = FALSE)
    } else {
      chronology6 <- NULL
    }
    chronology <- list(chronology1, chronology2, chronology3, chronology4, chronology5, chronology6)
    chronology <- list.clean(chronology, function(x) length(x) == 0L, TRUE)
    if (length(chronology) != 0) {
      chronology <- data.frame(AgeDepthPollen.subset$Depth_Pollen, chronology, stringsAsFactors = FALSE)
    } else {
      chronology <- data.frame(AgeDepthPollen.subset$Depth_Pollen, stringsAsFactors = FALSE)
    }

    # Collect datings to plot
    zero <- data.frame(zero$depth, zero$age, zero$error, "", zero$age, zero$error, stringsAsFactors = FALSE)
    colnames(zero) <- names(calibration)[1:6]
    if (parameter.subset$Add.surface) {
      calibration <- rbind(calibration, zero)
    }

    # Values for ggplot without a chronology
    if (dim(chronology)[2] >= 1) {
      y.plot.1 <- chronology[, 1]
      y.plot.2 <- chronology[, 1]
      y.plot.3 <- chronology[, 1]
      y.plot.4 <- chronology[, 1]
      y.plot.5 <- chronology[, 1]
      color.plot.1 <- "meanAgeBP"
      color.plot.2 <- "meanAgeBP"
      color.plot.3 <- "meanAgeBP"
      color.plot.4 <- "meanAgeBP"
      color.plot.5 <- "meanAgeBP"
      alpha.plot.1 <- 0
      alpha.plot.2 <- 0
      alpha.plot.3 <- 0
      alpha.plot.4 <- 0
      alpha.plot.5 <- 0
      naming <- list(plot.c14.correct, plot.AWI)
      color.names <- c("Calibrated.Age" = "red", "meanAgeBP" = "blue")
    }

    # With one chronology
    if (dim(chronology)[2] >= 4) {
      y.plot.1 <- chronology[, 3]
      color.plot.1 <- "Neotoma.1"
      alpha.plot.1 <- 0.8
      naming <- list(plot.c14.correct, plot.AWI, chronology[which(!is.na(chronology[, 2])), 2][1])
      color.names <- c(color.names, "Neotoma.1" = "green")
    }

    # With two chronologies
    if (dim(chronology)[2] >= 7) {
      y.plot.2 <- chronology[, 6]
      color.plot.2 <- "Neotoma.2"
      alpha.plot.2 <- 0.8
      naming <- list(plot.c14.correct, plot.AWI, chronology[which(!is.na(chronology[, 2])), 2][1], chronology[which(!is.na(chronology[, 5])), 5][1])
      color.names <- c(color.names, "Neotoma.2" = "yellow")
    }

    # With three chronologies
    if (dim(chronology)[2] >= 10) {
      y.plot.3 <- chronology[, 9]
      color.plot.3 <- "Neotoma.3"
      alpha.plot.3 <- 0.8
      naming <- list(plot.c14.correct, plot.AWI, chronology[which(!is.na(chronology[, 2])), 2][1], chronology[which(!is.na(chronology[, 5])), 5][1], chronology[which(!is.na(chronology[, 8])), 8][1])
      color.names <- c(color.names, "Neotoma.3" = "#FFCCFF")
    }

    # With four chronologies
    if (dim(chronology)[2] >= 13) {
      y.plot.4 <- chronology[, 12]
      color.plot.4 <- "Neotoma.4"
      alpha.plot.4 <- 0.8
      naming <- list(plot.c14.correct, plot.AWI, chronology[which(!is.na(chronology[, 2])), 2][1], chronology[which(!is.na(chronology[, 5])), 5][1], chronology[which(!is.na(chronology[, 8])), 8][1], chronology[which(!is.na(chronology[, 11])), 11][1])
      color.names <- c(color.names, "Neotoma.4" = "orange")
    }

    # With five chronologies
    if (dim(chronology)[2] >= 16) {
      y.plot.5 <- chronology[, 12]
      color.plot.5 <- "Neotoma.5"
      alpha.plot.5 <- 0.8
      naming <- list(plot.c14.correct, plot.AWI, chronology[which(!is.na(chronology[, 2])), 2][1], chronology[which(!is.na(chronology[, 5])), 5][1], chronology[which(!is.na(chronology[, 8])), 8][1], chronology[which(!is.na(chronology[, 11])), 11][1], chronology[which(!is.na(chronology[, 14])), 14][1])
      color.names <- c(color.names, "Neotoma.5" = "#FF00FF")
    }

    # ggplot
    output <- ggplot() +

      # curves
      geom_ribbon(data = model.AWI, aes(depth, ymin = min, ymax = max), alpha = 0.4) +
      geom_line(aes(x = chronology[, 1], y = y.plot.1, color = color.plot.1), alpha = alpha.plot.1) +
      geom_point(aes(x = chronology[, 1], y = y.plot.1, color = color.plot.1), alpha = min(alpha.plot.1, 0.6), shape = "|", size = 2) +
      geom_line(aes(x = chronology[, 1], y = y.plot.2, color = color.plot.2), alpha = alpha.plot.2) +
      geom_point(aes(x = chronology[, 1], y = y.plot.2, color = color.plot.2), alpha = min(alpha.plot.2, 0.6), shape = "|", size = 2) +
      geom_line(aes(x = chronology[, 1], y = y.plot.3, color = color.plot.3), alpha = alpha.plot.3) +
      geom_point(aes(x = chronology[, 1], y = y.plot.3, color = color.plot.3), alpha = min(alpha.plot.3, 0.6), shape = "|", size = 2) +
      geom_line(aes(x = chronology[, 1], y = y.plot.4, color = color.plot.4), alpha = alpha.plot.4) +
      geom_point(aes(x = chronology[, 1], y = y.plot.4, color = color.plot.4), alpha = min(alpha.plot.4, 0.6), shape = "|", size = 2) +
      geom_line(aes(x = chronology[, 1], y = y.plot.5, color = color.plot.5), alpha = alpha.plot.5) +
      geom_point(aes(x = chronology[, 1], y = y.plot.5, color = color.plot.5), alpha = min(alpha.plot.5, 0.6), shape = "|", size = 2) +
      geom_line(data = model.AWI, aes(x = depth, y = mean, color = "meanAgeBP")) +
      geom_point(data = AgeDepthPollen.subset, aes(x = Depth_Pollen, y = mean, color = "meanAgeBP"), alpha = 0.9, shape = "|", size = 2) +

      # datings
      geom_errorbar(data = calibration, aes(x = depth, ymin = cal.age - cal.age.se, ymax = cal.age + cal.age.se, color = "Calibrated.Age"), width = (max(calibration$depth) - min(calibration$depth)) * 0.01, alpha = 1) +
      geom_point(data = calibration, aes(x = depth, y = cal.age, color = "Calibrated.Age"), shape = 1, alpha = 1) +
      geom_linerange(aes(x = cal.sp[, 1], ymin = cal.sp[, 2], ymax = cal.sp[, 3], color = "Calibrated.Age")) +

      # lengend, theme and colour
      labs(title = paste0(AgeDepthPollen.subset$Site_Name[1], " | ID ", ID, " | long ", AgeDepthPollen.subset$Longitude..DD.[1], "? lat ", AgeDepthPollen.subset$Latitude..DD.[1], "? | elev ", AgeDepthPollen.subset$Elevation..m.a.s.l.., "\n", "acc.rate  ", round(acc.rate, digits = 2), " a/cm  ", round(1 / acc.rate, digits = 4), " cm/a | resolution ", round(if (is.na(parameter.subset$Resolution.cm)) {
        (max(AgeDepthPollen.subset$Depth_Pollen, calibration$depth) - parameter.subset$Waterline) / parameter.subset$Resolution.section
      } else {
        parameter.subset$Resolution.cm
      }, digits = 2), " cm | mem.mean ", parameter.subset$Memory, " | reservoir ", parameter.subset$Reservoir, " a  ")) +
      labs(x = "depth (cm)", y = "years BP") +
      scale_color_manual(name = "", values = color.names, labels = naming) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title = element_text(size = 10),
        legend.position = "bottom"
      ) +
      theme(
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 10)
      ) +
      theme(
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 15, b = 00, l = 0)),
        axis.text.y = element_text(size = 10)
      )

    # Write png
    output
    ggsave(paste0(folder, "/Sites/", ID, "/", ID, ".png"), width = 297, height = 210, units = "mm", dpi = 300)

    # Flipped version
    output.flipped <- output +
      coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
    output.flipped
    ggsave(paste0(folder, "/Sites/", ID, "/", ID, ".flipped.png"), width = 297, height = 210, units = "mm", dpi = 300)

    # Copy plots to new subfolders

    file.copy(from = paste0(folder, "/Sites/", ID, "/", ID, ".png"), to = paste0(folder, "/Plot.png/", ID, ".png"), overwrite = TRUE)
    file.copy(from = paste0(folder, "/Sites/", ID, "/", ID, ".flipped.png"), to = paste0(folder, "/Plot.flipped/", ID, ".flipped.png"), overwrite = TRUE)
  }

  #-----complete ID run-----


  # delete the ID from the list for the outer loop and go to the next ID
  parameter.ID <- parameter.ID[-1]
}
