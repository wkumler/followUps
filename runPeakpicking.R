#Layout:
# Create a noisy data file
# XCMS things
# MSDIAL things
# MzMine things

#Strategy
# Create a noisy data .mzML file and place in processingThings/
# Run XCMS on it and output peaklist to main directory
# Convert .mzML to .abf in processingThings/ (MANUAL STEP)
# Run MSDIAL on it and output peaklist to main directory
# Run MzMine on it and output peaklist to main directory



# Create a noisy data file in clean processing directory ----
if(dir.exists("processingThings/")){unlink('processingThings', recursive=TRUE)}
dir.create("processingThings")
library(MSnbase)
writeSingleFile <- function(ms_data_frame, outfile, colnames=c("rt", "mz", "intensity")){
  ms_data_frame$file <- 1
  scan_list <- split(ms_data_frame, ms_data_frame$rt)
  #Turn each scan into its own spectrum
  scan_spectra <- lapply(seq_along(scan_list), function(x){
    new("Spectrum1", 
        rt=unique(scan_list[[x]]$rt), 
        mz=scan_list[[x]]$mz, 
        intensity=scan_list[[x]]$intensity, 
        fromFile=as.integer(unique(scan_list[[x]]$file)), 
        acquisitionNum=x,
        scanIndex=x,
        polarity=NA_integer_)
  })
  #Add the names
  names(scan_spectra) <- sapply(scan_list, function(x){
    scannum <- as.numeric(unique(x$rt))
    scannum <- paste0(paste0(numeric(4-nchar(scannum)), collapse = ""), scannum)
    paste0("F", unique(x$file), ".S", scannum)
  })
  #Convert to environment to make MSnExp object happy
  scan_env <- as.environment(scan_spectra)
  
  #Create new object and fill with info
  y <- new("MSnExp")
  y@processingData@files <- "Artificial1"
  y@assayData <- scan_env
  y@featureData <- AnnotatedDataFrame(
    data.frame(spectrum=1:length(scan_spectra),
               polarity=rep(1, length(y)),
               totIonCurrent=sapply(scan_spectra, function(x)sum(x@intensity), USE.NAMES = FALSE),
               basePeakIntensity=sapply(scan_spectra, function(x)max(x@intensity), USE.NAMES = FALSE),
               injectionTime=rnorm(length(y), mean = 50, sd = 5),
               filterString=rep("FTMS + p ESI Full ms [100.0000-101.0000]", length(y)),
               spectrumId=paste0("controllerType=0 controllerNumber=1 scan=", seq_along(y)),
               centroided=rep(TRUE, length(y)),
               scanWindowLowerLimit=rep(100, length(y)),
               scanWindowUpperLimit=rep(101, length(y)),
               row.names = names(scan_spectra), stringsAsFactors = FALSE))
  fData(y)$injectionTime <- 0
  #Write out to file
  # if(file.exists(outfile)){
  #   readline(prompt = "File already exists. Press Enter to continue, Esc to abort.")
  #   file.remove(outfile)
  # }
  file.remove(outfile)
  writeMSData(y, file = outfile, outformat="mzml", verbose=T)
}
varNoiseMaker <- function(noise, total_width=1000, leaveNaughts=TRUE){
  pints <- dnorm(x = seq(-3, 3, length.out=50))
  pints <- 1000000/max(pints)*pints
  noisevec <- rnorm(length(pints), mean = 0, sd=noise*100)
  pints <- pints+noisevec
  if(leaveNaughts){
    pints[pints<0] <- 0
  }
  buffer <- numeric((total_width-length(pints))/2)
  buffer <- runif((total_width-length(pints))/2)*10000
  pints <- c(buffer, pints, buffer)
  
  prt <- 1:length(pints)
  pmz <- rep(100, length(pints))+noise/10000
  full_df <- data.frame(intensity=pints, rt=prt, mz=pmz)
  if(!leaveNaughts){
    full_df <- subset(full_df, full_df$intensity>0)
  }
  return(full_df)
}
noise_df <- do.call(rbind, lapply(seq(0, 10000-10, 10), varNoiseMaker, leaveNaughts=TRUE))
writeSingleFile(noise_df, outfile = "noisy_data.mzML")
file.copy("noisy_data.mzML", "processingThings/noisy_data.mzML")
message("Remember to convert that to .abf, now!")


# XCMS things ----
library(xcms)
xcms_object <- readMSData(files = "processingThings/noisy_data.mzml", 
                          msLevel. = 1, mode = "onDisk")
cwp <- CentWaveParam(snthresh = 0, prefilter = c(0,0), ppm = 0.1, fitgauss = TRUE)
xcms_output <- findChromPeaks(xcms_object, param = cwp)
xcms_peaks <- as.data.frame(chromPeaks(xcms_output))
xcms_peaks <- xcms_peaks[order(xcms_peaks$mz),]
write.csv(xcms_peaks, file = "xcms_peaks.csv")



# MZmine things ----
mzmine_param_maker <- function(){
  writeLines('<?xml version="1.0" encoding="UTF-8"?><batch>
    <batchstep method="io.github.mzmine.modules.io.rawdataimport.RawDataImportModule">
        <parameter name="Raw data file names">
            <file>G:\\My Drive\\followUps\\processingThings\\noisy_data.mzML</file>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scans">
            <ms_level>1</ms_level>
        </parameter>
        <parameter name="Mass detector" selected="Wavelet transform">
            <module name="Centroid">
                <parameter name="Noise level"/>
            </module>
            <module name="Exact mass">
                <parameter name="Noise level"/>
            </module>
            <module name="Local maxima">
                <parameter name="Noise level"/>
            </module>
            <module name="Recursive threshold">
                <parameter name="Noise level"/>
                <parameter name="Min m/z peak width"/>
                <parameter name="Max m/z peak width"/>
            </module>
            <module name="Wavelet transform">
                <parameter name="Noise level">0.0</parameter>
                <parameter name="Scale level">20</parameter>
                <parameter name="Wavelet window size (%)">0.8</parameter>
            </module>
        </parameter>
        <parameter name="Mass list name">masses</parameter>
        <parameter name="Output netCDF filename (optional)" selected="false"/>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_ADAPchromatogrambuilder.ADAPChromatogramBuilderModule">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scans">
            <ms_level>1</ms_level>
        </parameter>
        <parameter name="Mass list">masses</parameter>
        <parameter name="Min group size in # of scans">5</parameter>
        <parameter name="Group intensity threshold">0.0</parameter>
        <parameter name="Min highest intensity">0.0</parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>1.0E-4</absolutetolerance>
            <ppmtolerance>0.1</ppmtolerance>
        </parameter>
        <parameter name="Suffix">chromatograms</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_gridmass.GridMassModule">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scans">
            <ms_level>1</ms_level>
        </parameter>
        <parameter name="Suffix">chromatograms</parameter>
        <parameter name="Minimum height">20.0</parameter>
        <parameter name="M/Z Tolerance">0.0001</parameter>
        <parameter name="Min-max width time (min)">
            <min>0.1</min>
            <max>3.0</max>
        </parameter>
        <parameter name="Smoothing time (min)">0.05</parameter>
        <parameter name="Smoothing m/z">0.05</parameter>
        <parameter name="False+: Intensity similarity ratio">0.5</parameter>
        <parameter name="False+: Ignore times">0-0</parameter>
        <parameter name="Debugging level">No debug</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.csvexport.CSVExportModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Filename">
            <current_file>G:\\My Drive\\followUps\\mzmine_peaks.csv</current_file>
        </parameter>
        <parameter name="Field separator">,</parameter>
        <parameter name="Export common elements">
            <item>Export row ID</item>
            <item>Export row m/z</item>
            <item>Export row retention time</item>
            <item>Export row identity (main ID)</item>
            <item>Export row identity (all IDs)</item>
            <item>Export row identity (main ID + details)</item>
            <item>Export row comment</item>
            <item>Export row number of detected peaks</item>
        </parameter>
        <parameter name="Export data file elements">
            <item>Peak status</item>
            <item>Peak name</item>
            <item>Peak m/z</item>
            <item>Peak RT</item>
            <item>Peak RT start</item>
            <item>Peak RT end</item>
            <item>Peak duration time</item>
            <item>Peak height</item>
            <item>Peak area</item>
            <item>Peak charge</item>
            <item>Peak # data points</item>
            <item>Peak FWHM</item>
            <item>Peak tailing factor</item>
            <item>Peak asymmetry factor</item>
            <item>Peak m/z min</item>
            <item>Peak m/z max</item>
        </parameter>
        <parameter name="Export quantitation results and other information">false</parameter>
        <parameter name="Identification separator">;</parameter>
        <parameter name="Filter rows">ALL</parameter>
    </batchstep>
</batch>
', con="G:/My Drive/followUps/processingThings/mzmine_params_auto.xml")
}
mzmine_param_maker()
file.create("G:/My Drive/followUps/mzmine_peaks.csv")
system(paste('"C:/Program Files/MZmine-2.53-Windows/startMZmine-Windows.bat"',
             '"G:/My Drive/followUps/processingThings/mzmine_params_auto.xml"'),
       show.output.on.console = TRUE)



# MSDIAL things ----
if(!file.exists("processingThings/noisy_data.abf")){stop("Remember to convert to abf!")}
file.remove("processingThings/noisy_data.mzML")
msdial_param_maker <- function(){
  writeLines('#Data type
MS1 data type: Centroid
MS2 data type: Centroid
Ion mode: Positive
DIA file: 

#Data collection parameters
Retention time begin: 0
Retention time end: 100
Mass range begin: 0
Mass range end: 2000

#Centroid parameters
MS1 tolerance for centroid: 0.01
MS2 tolerance for centroid: 0.05

#Peak detection parameters
Smoothing method: LinearWeightedMovingAverage
Smoothing level: 3
Minimum peak width: 5
Minimum peak height: 1
Mass slice width: 0.0001

#Deconvolution parameters
Sigma window value: 0.5
Amplitude cut off: 0

#Adduct list
Adduct list: [M+H]+

#MSP file and MS/MS identification setting
#MSP file: D:\\Msdial-ConsoleApp-Demo files\\Msdial-ConsoleApp-Demo files for DDA\\MSDIAL-LipidDB-VS23.msp
Retention time tolerance for identification: 4
Accurate ms1 tolerance for identification: 0.025
Accurate ms2 tolerance for identification: 0.25
Identification score cut off: 70

#Text file and post identification (retention time and accurate mass based) setting
#Text file: D:\\Msdial-ConsoleApp-Demo files\\Msdial-ConsoleApp-Demo files for DDA\\Lipid_Nega_IS_PostIdentification_vs1.txt
Retention time tolerance for post identification: 0.5
Accurate ms1 tolerance for post identification: 0.01
Post identification score cut off: 85

#Alignment parameters setting
Retention time tolerance for alignment: 0.05
MS1 tolerance for alignment: 0.025
Retention time factor for alignment: 0.5
MS1 factor for alignment: 0.5
Peak count filter: 0
QC at least filter: True
', con = "G:/My Drive/followUps/processingThings/msdial_params_auto.txt")
}
msdial_param_maker()
system('"C:/Program Files/MSDIAL ver 4.12/MsdialConsoleApp.exe" lcmsdda 
       -i "G:/My Drive/followUps/processingThings" 
       -o "G:/My Drive/followUps/processingThings" 
       -m "G:/My Drive/followUps/processingThings/msdial_params_auto.txt"')
write.csv(read.csv("processingThings/noisy_data.msdial", sep = "\t"), 
          file = "msdial_peaks.csv")


# Spacer ----




