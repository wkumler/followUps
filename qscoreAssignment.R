# Assign quality scores to each kind of peakpicked software by looking
# at the raw data and testing a correlation and estimating signal-to-noise

# Setup things ----
library(pbapply)

# Functions ----
grabSingleFileData <- function(filename, rt_lim=NULL){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  if(length(rt_lim)){ #If rt_lim isn't NULL, subset by rt_lim
    if(!length(rt_lim)==2){stop("rt_lim should be a vector of length 2")}
    if(!is.numeric(rt_lim)){stop("rt_lim must be numeric")}
    fullhd <- subset(fullhd, subset = fullhd$retentionTime>min(rt_lim)&
                       fullhd$retentionTime<max(rt_lim))
  }
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<3){
    return(0)
  }
  #Create an "ideal" peak of the same width
  perf_peak <- dnorm(seq(-3, 3, length.out = nrow(eic)))
  #Calculate the correlation between the perfect peak and the observed values
  peak_cor <- cor(perf_peak, eic$int)
  #Calculate the normalized residuals
  residuals <- eic$int/max(eic$int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(eic$int)-min(eic$int))/sd(norm_residuals*max(eic$int))
  #Return the quality score
  return(SNR*peak_cor^10*log10(max(eic$int)))
}
xcmsQscoreCalculator <- function(peak_row, peak_data, raw_data){
  #Extract the relevant EIC
  peak_row_data <- peak_data[peak_row, ]
  eic <- subset(raw_data, raw_data$rt>=peak_row_data$rtmin&raw_data$rt<=peak_row_data$rtmax&
                  raw_data$mz>=peak_row_data$mzmin&raw_data$mz<=peak_row_data$mzmax)
  return(qscoreCalculator(eic))
}
msdialQscoreCalculator <- function(peak_row, peak_data, raw_data, ppm){
  peak_row_data <- peak_data[peak_row, ]
  eic <- subset(raw_data, raw_data$rt>=peak_row_data$RT.left.min.*60&
                  raw_data$rt<=peak_row_data$RT.right..min.*60&
                  raw_data$mz>=peak_row_data$Precursor.m.z*(1-ppm/1000000)&
                  raw_data$mz<=peak_row_data$Precursor.m.z*(1+ppm/1000000))
  return(qscoreCalculator(eic))
}
mzmineQscoreCalculator <- function(peak_row, peak_data, raw_data, ppm){
  peak_row_data <- peak_data[peak_row, ]
  eic <- subset(raw_data, raw_data$rt>=peak_row_data$Peak.RT.start*60&
                  raw_data$rt<=peak_row_data$Peak.RT.end*60&
                  #raw_data$mz>=peak_row_data$Peak.m.z.min&
                  #raw_data$mz<=peak_row_data$Peak.m.z.max
                  raw_data$mz>=peak_row_data$row.m.z*(1-ppm/1000000)&
                  raw_data$mz<=peak_row_data$row.m.z*(1+ppm/1000000))
  return(qscoreCalculator(eic))
}

# Raw data extraction ----
raw_data <- grabSingleFileData(filename = "noisy_data.mzML")


# XCMS things ----
xcms_peaks <- read.csv("xcms_peaks.csv")
head(xcms_peaks)
plot(xcms_peaks$mz, xcms_peaks$sn)
xcms_qscores <- pbsapply(seq_len(nrow(xcms_peaks)), xcmsQscoreCalculator,
                       peak_data=xcms_peaks, raw_data=raw_data)
plot(xcms_peaks$mz, xcms_qscores)
plot(xcms_peaks$mz[xcms_peaks$sn<2000], xcms_peaks$sn[xcms_peaks$sn<2000])
plot(xcms_peaks$mz[xcms_qscores>1], xcms_qscores[xcms_qscores>1])

# MSDIAL things ----
msdial_peaks <- read.csv("msdial_peaks.csv")
head(msdial_peaks[-(ncol(msdial_peaks)-1)])
plot(msdial_peaks$Precursor.m.z, msdial_peaks$S.N)
msdial_qscores <- pbsapply(seq_len(nrow(msdial_peaks)), msdialQscoreCalculator,
                           peak_data=msdial_peaks, raw_data=raw_data, ppm=2.5)
plot(msdial_peaks$Precursor.m.z, msdial_qscores)
plot(msdial_peaks$Precursor.m.z[msdial_qscores>1], msdial_qscores[msdial_qscores>1])


# MzMine things ----
mzmine_peaks <- read.csv("mzmine_peaks.csv")
names(mzmine_peaks) <- gsub("noisy_data.mzML.", "", names(mzmine_peaks))
head(mzmine_peaks)
mzmine_qscores <- pbsapply(seq_len(nrow(mzmine_peaks)), mzmineQscoreCalculator,
                           peak_data=mzmine_peaks, raw_data=raw_data, ppm=2.5)
plot(mzmine_peaks$row.m.z, mzmine_qscores)
plot(mzmine_peaks$row.m.z[mzmine_qscores>1], mzmine_qscores[mzmine_qscores>1])



# Interpretation ----
# Plot the "good" quality scores for each
par(mfrow=c(3,1))
par(mar=c(2.1,4.1,0.1,0.1))
plot(xcms_peaks$mz[xcms_qscores>1], 
     xcms_qscores[xcms_qscores>1],
     xlim=c(100, 101), ylab="")
plot(msdial_peaks$Precursor.m.z[msdial_qscores>1], 
     msdial_qscores[msdial_qscores>1],
     xlim=c(100, 101), ylab="Quality score of 'good' peaks", xpd=NA)
plot(mzmine_peaks$row.m.z[mzmine_qscores>1],
     mzmine_qscores[mzmine_qscores>1],
     xlim=c(100, 101), ylab="")

# Plot them in log scale
plot(xcms_peaks$mz[xcms_qscores>1], 
     log(xcms_qscores[xcms_qscores>1]),
     xlim=c(100, 101), ylab="")
plot(msdial_peaks$Precursor.m.z[msdial_qscores>1], 
     log(msdial_qscores[msdial_qscores>1]),
     xlim=c(100, 101), ylab="Quality score of 'good' peaks", xpd=NA, xlab="")
plot(mzmine_peaks$row.m.z[mzmine_qscores>1],
     log(mzmine_qscores[mzmine_qscores>1]),
     xlim=c(100, 101), ylab="")
layout(1)
par(mar=c(4.1,5.1,2.1,2.1))

# Some heatmap things!
possible_peaks <- round(seq(100,101,0.001), digits = 3)
xcms_found <- possible_peaks%in%round(xcms_peaks$mz[xcms_qscores>1], digits = 3)
msdial_found <- possible_peaks%in%round(msdial_peaks$Precursor.m.z[msdial_qscores>1], digits = 3)
mzmine_found <- possible_peaks%in%round(mzmine_peaks$row.m.z[mzmine_qscores>1], digits = 3)
# Show a heatmap thing of all the picked peaks in each one
image(cbind(xcms_found, msdial_found, mzmine_found), axes=FALSE)
axis(1, at = pretty(possible_peaks)-100, labels = pretty(possible_peaks))
axis(2, at = c(0, 0.5, 1), labels = c("XCMS", "MS-DIAL", "MzMine"), las=1)
mtext("Mass (aka signal-to-noise)", side = 1, line = 2.5, at = 0.5)
# Show a barplot of which peaks were picked by all 3
barplot(table(cut(possible_peaks[which(xcms_found&msdial_found&mzmine_found)], 
                  breaks = seq(100,101,0.01))))


# Look at the sketchiest peaks they all liked
peakCheck <- function(mzr, df=raw_data, mars=TRUE){
  eic <- subset(df, df$mz==mzr&df$rt>400&df$rt<600)
  if(!mars){
    par(mar=c(0,0,0,0));ylabel="";yaxtype="n";xaxtype="n"
  } else {
    ylabel=NULL;yaxtype=NULL;xaxtype=NULL
  }
  plot(eic$rt, eic$int, type="l", lwd=2, ylab=ylabel, yaxt=yaxtype, xaxt=xaxtype)
  if(!mars){par(mar=c(4.1, 4.1, 2.1, 2.1))}
}
peakCheck(100.329)
