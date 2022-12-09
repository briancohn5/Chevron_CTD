# CREATED DATE: 13 March 2018
# MODIFIED DATE: 19 May 2020
# 
# AUTHORS: Matt Robart mrobart@oxy.edu
# 
# this script was originally created by Matt Robart in the Matlab environment
# and was subsequently converted into R for VRG use
# 
# the May 2020 modification was dramatic so I changed the name of the script so 
# we can preserve the previous version. Now the averaging is all done in dplyr
# pipes, so theoretically it should be faster, and the code is now much easier to 
# read. The user must set the year, bin size, and whether the output will be in
# feet. The depth units are included in the final *.csv file, so pay attention!
# 
# This script is used to identify the downcast and clip CTD data files. It saves 
# the clipped data into a *.csv file that includes the date, time, and serial
# number in the filename. It also averages the data into user-defined depth
# bins and saves only the columns used for the Chevron reports (depth, temp (F),
# DO (ml/l), pH, and conductivity (mmhos/cm)).
# 
# This script uses files that have undergone initial processing by SBE Data
# Processing, then split into ascii files with a separate header file using the 
# module "ASCII out."
# 
# SBE Data Processing modules used for the Chevron project are:
# 1. Data Conversion
#       Scan Count
#       Julian Days
#       Pressure (db)
#       Temperature (ITS-90, degF)
#       Temperature (ITS-90, degC)
#       Conductivity (S/m)
#       Oxygen raw (V)
#       pH
#       Beam Transmission, WET Labs C-Star (%)
#       PAR/Logararithmic, Satlantic (umol phontons/m^2/sec)
# 2. Filter
#       LP filter A: 0.1s (Pressure only)
#       LP filter B: 0.5s (everything else)
# 3. Align CTD
#       T: 0.55s
#       DO: 1.5s
# 6. Derive
#       Depth (salt water, m)
#       Salinity, Practical (PSU)
#       Density (sigma-t, kg/m^3)
#       Oxygen, SBE 43 (ml/l)
# 15. ASCII Out
#       all variables are output
#       
# After the SBE Data Processing steps are completed, there should be a separate
# *.asc file for each cast. The filenames do not matter at all at this point.
# 
# Put all the *.asc files into a working folder for processing.
# 
# USER CAN ZOOM IN TO IDENTIFY THE END OF SURFACE SOAK AND END OF DOWNCAST 
# 
# 
# Using depth through time, the surface soak will be shallow depth at
# the beginning of the record.  The downcast begins when the depth
# begins to increase, either after the CTD was raised to the surface
# following the soak, or dropped directly after the soak.  In the
# second case the surface bin is lost unless the end user wants to try
# and pull it out of the upcast data, which gets messy.....
# anyway, the procedure here is to identify the scan numbers of the
# first and last good scans of the downcast only
# 
# on 12 Nov. 2019: added beam transmission and PAR


library(ggplot2)
library(dplyr)
library(plotly)

################ REQUIRES USER INPUT: CHANGE BEFORE RUNNING SCRIPT ############
# set working directory
setwd("~/VRG Files/R Code/CTD/data") 
# set serial number, year of sampling, and size of depth bins (in meters)
sn <- 7157
year <- 2022
binsize <- 1 # use 1-m bins for Chevron
feet <- FALSE # to convert depth to feet, set to TRUE
################ YOU SHOULDN'T HAVE TO CHANGE MUCH AFTER THIS POINT ###########


# turn the sampling year into a POSIXt date-time class
x <- paste(year,1,1,sep="-")
nt <- as.POSIXct(x,"%Y-%m-%d")

# create list of files in the working directory with *.asc filenames
f <- list.files(pattern="*.asc")

for (i in 1:length(f)) {
  dat <- read.table(f[i]) # read in all the data from each cast, one at a time
  # name all the columns of data (this will be helpful later)
  names(dat) <- c("Scan_Count", "Julian_Days", "Pressure", "temp_degF", 
                  "temp_degC", "Cond", "DO_raw", "pH", "light_transmission_25cm", "par", "depth",
                  "sal_PSU", "dens", "DO_ppt", "flag")
  #dat$Cond <- dat$Cond*10 # convert from S/m to mS/cm
  # plot depth (y-axis) vs scan count (x-axis):
  if(!is.null(dev.list())) dev.off() # I can't remember why I had to put this part in
  
  # make the plot of depth vs. scan count
  p <- ggplot(dat, aes(Scan_Count, depth)) +
    geom_point() +
    ggtitle(f[i])
  print(ggplotly(p)) # makes plot that you can zoom in on to find stop/start points
  
  # forces the user to choose the beginning and ending scans for the downcast
  print('Enter beginning index, then hit enter twice:'); start <- scan()
  print('Enter ending index, then hit enter twice:'); stop <- scan()
  
  # subset the downcast portion of the cast to be saved:
  downcast <- dat[start:stop,]
  
  # convert depth to feet if the chosen output is binned by feet
  if (feet == TRUE) {
    downcast$depth <- downcast$depth * 3.28084 # convert meters to feet
    depthunit <- "feet" # string to put into output file
  } else {
    depthunit <- "meters" # string to put into output file
    }

  # create depth bins and average data into them
  binned <- downcast %>% # choose the data frame to be binned
    mutate(bin_number = floor(depth/binsize)) %>% # create bins
    group_by(bin_number) %>% # group the data by bins
    add_count(bin_number) %>% # count how many scans are in each bin
    summarize(tempF = mean(temp_degF),
              DO_ppt = mean(DO_ppt),
              pH = mean(pH),
              cond = mean(Cond),
              tempC = mean(temp_degC),
              sal = mean(sal_PSU),
              light_transmission_25cm = mean(light_transmission_25cm),
              par = mean(par),
              count = max(n)) %>% # do all the averaging, and include the row count as a column
    mutate(depth = bin_number * binsize, depthunit = depthunit) %>% # convert bin number back to a real depth (helpful if it's not 1) and include a column showing the depth units
    select(bin_number, depthunit, depth, everything()) # reorder the columns
  
  dt <- dat$Julian_Days[1] # PDT time (use from mid March to early November, find the day/time of the cast (in decimal yearday)
  #dt <- dat$Julian_Days[1] - (1/24) # PST time (use from early November to mid March the following year), find the day/time of the cast (in decimal yearday)
  casttime <- nt + (86400*(dt-1)) # convert yearday to seconds and add the current year to find the actual date/time of the cast
  # 
  # # export the data as a *.csv
  setwd('../binned') # go to the parent directory and then into /binned
  
  # write the *.csv file 
  ####Change the PDT/PST identifier based on time of year###
  write.csv(binned, file=paste("CTD",as.character(casttime, format = "%Y%m%d__%H%M"),
                               "PDT_sn",sn,"_binned.csv",sep = ""))
  
  setwd('../data') # go to the parent directory, then back into /data, then repeat the loop for the next cast
}
