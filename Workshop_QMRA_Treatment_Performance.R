graphics.off() 
rm(list=ls(all=TRUE))
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(tidyr)
library(vtable) 

###Quantification of pathogen reduction across a treatment barrier: empirical approach###

#Exemple adapted from the book: Haas (2023). Doing Engineering Calculations.

#Background: 
# The protozoan parasite Cryptosporidium can be waterborne and survives because it is excreted by humans and warm-blooded animals as resistant cysts.
# It can be found in wastewater. 
# At a certain wastewater treatment plant, the concentration of the organism was measured periodically over the course of a year, both in the raw wastewater and after secondary biological treatment.
# The reduction of oocysts through treatment is of particular interest in characterizing performance.

#Data import
CryptosporidiumData = read.csv2( file="Cryptosporidium.csv" )
# To be able to read the file "Cryptosporidium.csv", download it from Github and put it in a specific folder. 
# Then, in R studio, go to Session -> Set Working Directory -> Choose Directory, and select your specific folder.
# You should then be able to read the file by running line 21.

Count_Raw = as.numeric(CryptosporidiumData[, "Count_inflow"]) #oocyst
Volume_Raw = as.numeric(CryptosporidiumData[, "Volume_inflow"]) #L
Count_Secondary = as.numeric(CryptosporidiumData[, "Count_outflow"]) #oocyst
Volume_Secondary = as.numeric(CryptosporidiumData[, "Volume_outflow"]) #L

#Concentrations in raw wastewater and secondary effluent
Concentration_Raw <- Count_Raw/Volume_Raw #oocyst/L
Concentration_Secondary <- Count_Secondary/Volume_Secondary #oocyst/L
log_Concentration_Raw <-log(Concentration_Raw)
log_Concentration_Secondary <-log(Concentration_Secondary)

# Add a sample identifier 
CryptosporidiumData$SampleID <- 1:nrow(CryptosporidiumData)

# Add log10 concentrations to the dataset
CryptosporidiumData$log_Concentration_Raw <- log_Concentration_Raw
CryptosporidiumData$log_Concentration_Secondary <- log_Concentration_Secondary

# Convert data to long format for creating plot
CryptosporidiumData_long <- pivot_longer(CryptosporidiumData, 
                                         cols = c("log_Concentration_Raw", "log_Concentration_Secondary"),
                                         names_to = "Type",
                                         values_to = "Log_Concentration")

# Create the plot 
ggplot(CryptosporidiumData_long, aes(x = SampleID, y = Log_Concentration, color = Type, shape = Type)) +
  geom_point() + 
  labs(x = "Sample Order",
       y = "Log Concentration") 

# ------------------------

# Check for correlation 
rho <- cor.test (log_Concentration_Raw, log_Concentration_Secondary, method ="pearson")
print(rho)

# -> Influent and effluent data are not correlated, and can thus be assumed statistically independant.

# ------------------------

# Histogram of natural log Giardia concentrations
par(mfrow = c(1 ,2))
hist(log_Concentration_Raw)
hist (log_Concentration_Secondary)

# Summary statistics of natural log Giardia concentrations
logData <-data.frame(log_Concentration_Raw, log_Concentration_Secondary)
sumtable(logData, add.median = TRUE )

# ------------------------
# Test of Normality
print(shapiro.test(log_Concentration_Raw))
print( shapiro.test(log_Concentration_Secondary))

# -> Computing the ratio of two independant lognormally distributed random variables
#    would have been straightfoward as the ratio also follows a lognormal distribution.
#    However, log-transformed data are not normally distributed. 
#    To address this problem, we can use Monte Carlo simulations.

# ------------------------

# Monte Carlo Analysis for Log Removal
library(EnvStats) # For density plots
trials = 10000

# We construct two sets of 10,000 randome integers between <1 23>
indexRaw = sample.int(23, trials, replace = TRUE)
indexSecondary = sample.int(23, trials, replace = TRUE)

MCRaw <- log_Concentration_Raw[indexRaw]
MCSecondary <- log_Concentration_Secondary[indexSecondary]

NaturalLogsReduced <- MCRaw - MCSecondary

par(mfrow = c(1, 1))
epdfPlot(NaturalLogsReduced, xlab =" Natural Log Removal ",
         ylab =" Probability Density ",
         main =" Empirical pdf ",
         curve.fill = TRUE)

# Microbial reduction performance are commonly reported in log10
log10_Concentration_Raw <-log10(Concentration_Raw)
log10_Concentration_Secondary <-log10(Concentration_Secondary)

MCRaw <- log10_Concentration_Raw[indexRaw]
MCSecondary <- log10_Concentration_Secondary[indexSecondary]

Logs10Reduced <-MCRaw - MCSecondary

epdfPlot(Logs10Reduced, xlab ="Log10 Removal",
         ylab =" Probability Density ",
         main =" Empirical pdf ",
         curve.fill = TRUE)

# For risk assessment purposed, the arithemtic mean probability of passage 
# is of interest (also known as effective log10 removal)
Removal <- 10^(-Logs10Reduced)
Arithmetic_mean_removal <- mean(Removal)
Effective_LRV <- -log10(Arithmetic_mean_removal)
# Add vertical lines for Effective LRV
abline(v = Effective_LRV, col = "blue", lwd = 2, lty = 1)

# ------------------------
