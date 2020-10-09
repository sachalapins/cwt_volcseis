# Example R script for producing Fig 1 of Lapins et al., 2020, JVGR, https://doi.org/10.1016/j.jvolgeores.2019.106728

# Packages required - if you don't have these packages installed, uncomment these first lines
#install.packages("viridis") # for scalogram colour palette
#install.packages("zoo") # for storing / subsetting time series
#install.packages("MASS") # for fractions axis labels

source("R/morlet_cwt.R") # Source function to compute CWT with Morlet wavelet
source("R/readSACbinary.R") # Source function to read in binary SAC file

# Read in example SAC file
sacfile <- "SAC/Kilauea_VT.SAC"
trace <- readSAC(sacfile)

# Check data read in okay:
plot(trace$x1, type='l', xlab="Time", ylab="Counts")

# Compute CWT
freq_high <- 50
freq_low <- 1/20
x_cwt <- morlet_cwt(x=coredata(trace$x1), dt=trace$delta, dj=1/64, lp=1/freq_high, up=1/freq_low, omega0=6, standardise_x = T, normalise_by_scale = T)

#### Plot CWT as per Fig 1 of JVGR paper ####
# Colour scheme for CWT plots
col <- rev(viridis::magma(100))

# Uncomment this line if you want to save plot as jpeg:
#jpeg(filename="Kilauea_VT.jpeg", units="in", width = 8.3, height = 3.9, res=600, bg="transparent")

# Start plotting. Set margins etcs:
par(omi=c(0,0,0,0), mgp=c(2,0.25,0), plt=c(0.12, 0.79, 0.15, 0.79))

# Plot scalogram with fourth-root scaling as per JVGR paper, and reverse to plot with increasing frequency (instead of increasing period) along y-axis
image(x=x_cwt$axis_x, y=x_cwt$axis_y, z=t(apply(sqrt(sqrt(x_cwt$power)), 2, rev)), ylim=c(max(x_cwt$axis_y), min(x_cwt$axis_y)), useRaster=T, col=col, xaxt="n", yaxt="n", xlab="", ylab="")

# Add COI
polygon(x_cwt$coi_x, x_cwt$coi_y, border = NA, col = rgb(1, 1, 1, 0.5))

# Add empty plot to do log y-axis and control tick marks:
par(new=T)
plot(x=x_cwt$axis_x, y=1/seq(from=min(x_cwt$periods), to=max(x_cwt$periods), length.out=length(x_cwt$axis_x)), log="y", type="n", xlim=c(0, max(x_cwt$axis_x)), ylim=c(freq_low, freq_high), xaxs="i", yaxs="i", xlab="", ylab="", cex.axis=0.8, yaxt="n", xaxt="n")

# Add x-axis:
axis(side=1, tck=-0.02, cex.axis=0.8)
title(xlab = "Time (s)", cex.lab = 1, line=1.25)

# Add y-axis:
tick_factor <- c(1,2,5) # Ticks at every tick_factor divided by some power of 10; e.g., 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, etc
myTicks <- c(1 / (rev(tick_factor) * 1000), 1 / (rev(tick_factor) * 100), 1 / (rev(tick_factor) * 10), 1 / rev(tick_factor), tick_factor, tick_factor * 10, tick_factor * 100)
myTicks <- unique(myTicks[myTicks >= freq_low & myTicks <= freq_high])
axis(side=2,at=myTicks,
     labels=ifelse(myTicks >= 1, sprintf("%.0f", myTicks),
                   ifelse(myTicks >= 0.1, sprintf("%0.1f", myTicks), 
                          ifelse(myTicks >= 0.01, sprintf("%0.2f", myTicks), 
                                 sprintf("%0.3f", myTicks)))), tck=0.02, cex.axis=0.8, las=1)
title(ylab = "Frequency (Hz)", cex.lab = 1, line = 1.75)

# Add box around scalogram to look more tidy:
box()

# Plot global wavelet spectrum
par(new=T, plt=c(0.79, 0.89, 0.15, 0.79)) # set plot margins
plot(x=rowMeans(x_cwt$power), y=x_cwt$freqs, ylim=c(freq_low, freq_high), log='y', type="l", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="") # Plot global wavelet spectrum
axis(side=4, at=1/myTicks, labels=MASS::fractions(myTicks), tck=0.02, cex.axis=0.8, las=1)
mtext("Period (s)", side=4, line=1.25)

# Plot signal above scalogram:
par(new=T, plt=c(0.12, 0.79, 0.79, 0.94)) # set margins
options(digits.secs = 2) # Give time stamp 2 decimal places
plot(x=x_cwt$axis_x, y=coredata(trace$x1), type="l", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", lwd=0.6) # plot signal
legend("topleft", paste(trace$kstnm, trace$kcmpnm, sep="."), bty="n", cex=0.6, pt.cex = 1) # Add station and component label to plot
time_stamp <- as.POSIXct(paste(trace$nzyear, " ", trace$nzjday, " ", trace$nzhour, ":", trace$nzmin, ":", trace$nzsec, ".", trace$nzmsec, sep=""), tz="utc", format="%Y %j %H:%M:%OS") # Create time stamp
legend("topright", strftime(time_stamp), bty="n", cex=0.6, pt.cex = 1) # Add time stamp to plot

# Uncomment this line to save plot as jpeg:
#dev.off()

