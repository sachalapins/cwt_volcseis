
morlet_cwt <- function(x, dt, dj, lp, up, omega0=6, standardise_x = T, normalise_by_scale = T) {

    ##################
    
    # Function to run continuous wavelet transform with complex Morlet wavelet
    # Code used in Lapins et al., 2020, JVGR (https://doi.org/10.1016/j.jvolgeores.2019.106728)
        
    # Function inputs:
    # x: input data (1D vector)
    # dt: delta (single numeric)
    # dj: 1 / number of 'voices' per 'octave' - this is effectively resolution, how many sub-scales between each integer wavelet scale
    #       By contrast, the DWT only uses integer wavelet scales, so dj parameter in CWT adds 'over-completeness' (not a bad thing!)
    # lp: lower period (highest frequency) of signal to compute CWT over
    # up: upper period (lowest frequency) of signal to compute CWT over
    # omega0: Morlet wavelet central frequency (see Eq 2.13 or Fig 3 in JVGR paper)
    # standardise_x: Standardise the input data before running CWT (e.g., to make spectra more comparable between very different event sizes)
    # normalise_by_scale: CWT has inherent bias towards longer period components - if TRUE, this bias is removed by dividing wavelet power by wavelet scale (Liu et al, 2007, https://doi.org/10.1175/2007JTECHO511.1)
    
    # Function outputs:
    # A list with:
    # scales: sequence of wavelet scales used
    # periods: corresponding sequence of periods in seconds (using wavelet's central frequency)
    # freqs: corresponding sequence of frequencies in Hz
    # wave: complex matrix holding CWT coefficients
    # power: matrix holding CWT power (scalogram)
    # axis_1: x-axis for plotting
    # axis_2: y-axis for plotting (R doesn't handle log-scaled images well)
    # coi_x: cone of influence x-coords for plotting
    # coi_y: cone of influence y-coords for plotting
    
    ###################
    
    #### Setup ####
    
    # Signal length and padding
    series_length <- length(x) # Series length
    next_dyadic <- 2^ceiling(log2(series_length)) # Next power of 2
    pad_length <- next_dyadic - series_length # Amount of zero-padding to reach next power of 2
    
    # Define characteristic frequency conversion (from wavelet scale to periods in sec)
    fourier_factor <- (2*pi)/omega0 # Use wavelet's central frequency to determine conversion between scale and frequency
    
    # There are other options for this conversion, using e.g. the peak frequency in wavelet's energy spectrum:
    #fourier_factor <- (4*pi)/(omega0+sqrt(2+omega0^2)) # See Table 1 of Torrence & Compo, 1998, https://doi.org/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2
    # or the centre of the passband of wavelet's spectrum (e.g., when using a wavelet with skewed spectrum)
    # but I find no noticeable difference between the above for Morlet wavelet
    # See Addison 2016, The Illustrated Wavelet Transform Handbook (2nd edition), Section 2.9 for more on this.
    
    # Compute wavelet scales and corresponding periods
    min_scale <- lp / fourier_factor # Convert lower period to minimum wavelet scale
    max_scale <- up / fourier_factor # Convert upper period to maximum wavelet scale
    
    min_scale_index <- ceiling(log2(min_scale) / dj) * dj
    max_scale_index <- floor(log2(max_scale) / dj) * dj
    
    scales <- 2^seq(from=min_scale_index, to=max_scale_index, by=dj) # Create sequence of wavelet scales
    scales_length <- length(scales) # No. of wavelet scales
    
    periods <- fourier_factor * scales # Sequence of periods (from sequence of wavelet scales)
    
    # Computation of the angular frequencies - See Eq 2.12 of JVGR paper
    M <- next_dyadic
    omega_k_pt1 <- 0:(M/2) # For k <= M/2
    omega_k_pt1 <- (omega_k_pt1 * 2 * pi)/(M * dt) # top part of Eq 2.12
    omega_k_pt2 <- ((M/2)+1):(M-1) # For k > M/2 (up to M-1)
    omega_k_pt2 <- -(omega_k_pt2 * 2 * pi) / (M * dt) # bottom part of Eq 2.12
    omega_k <- c(omega_k_pt1, omega_k_pt2) # Combine two parts
    
    #### Morlet Wavelet Transform ####
    
    # Standardise x and pad with zeros
    x <- x - mean(x) # Substract mean (reduces edge effect)
    if (standardise_x) {
        x <- x / sd(x) # Standardise
    }
    xpad <- c(x, rep(0, pad_length)) # Pad x with zeros
    
    # Compute FFT of xpad - for Eq 2.11 of JVGR paper
    fft_xpad <- fft(xpad)
    
    # Prepare complex matrix to accommodate CWT coefficients
    wave <- matrix(0, nrow=scales_length, ncol=next_dyadic) # Matrix of zeros
    wave <- wave + (1i * wave) # Matrix of complex zeroes
    
    # Compute CWT for each wavelet scale through inverse FFT (Eq 2.11 of JVGR paper)
    for (i in 1:scales_length) {
        my_scale <- scales[i]
        norm_unit_energy <- sqrt(2 * pi * my_scale / dt) # Normalisation to unit energy (Eq 6 of Torrence & Compo, 1998)
        norm_constant <- pi^(-1/4) # Normalisation constant for Morlet wavelet (e.g., Table 1 of Torrence & Compo, 1998)
        
        expnt1 <- -((my_scale * omega_k - omega0)^2 / 2) # Exponent for Fourier transform of Morlet wavelet, e.g. Table 1 of T&C 1998
        expnt2 <- -(((my_scale * omega_k)^2 + omega0^2) / 2) # Correction for non-zero mean (becomes important for values of omega0 < 5), not in T&C 1998
        daughter <- norm_unit_energy * norm_constant * (exp(expnt1) - exp(expnt2)) # Full Fourier transform of Morlet mother wavelet ("daughter" wavelet)
        daughter <- daughter * (as.numeric(omega_k > 0)) # Heaviside step function from Table 1 of T&C 1998
        
        # Compute CWT via inverse FFT of product of FT of signal and FT of wavelet (Eq 2.11 of JVGR paper)
        wave[i,] <- fft(fft_xpad * daughter, inverse=TRUE) / M # R's fft function doesn't normalise so we do
    }
    
    # Remove zero-padded part
    wave <- wave[,1:series_length]
    
    # Compute wavelet power (scalogram) - e.g., Fig 1 of JVGR paper
    power <- Mod(wave)^2
    if (normalise_by_scale) {
        power <- power / matrix(rep(scales, series_length), nrow=scales_length) # Divide power by wavelet scale
        wave <- wave / matrix(rep(sqrt(scales), series_length), nrow=scales_length) # Divide wavelet transform by sqrt of wavelet scale
    }
    
    axis_1 = seq(from = dt, by = dt, length.out = series_length) # x-axis for plotting CWT
    axis_2 = log2(periods) # y-axis for plotting CWT (log axes don't seem to work well with images in R)
    
    # Compute cone of influence (see Fig 1 of JVGR paper) - reworked from old WaveletCo pacakge (https://cran.r-project.org/package=WaveletCo)
    coi = fourier_factor * sqrt(2) * dt * c(1e-5, 1:((series_length + 1) / 2 - 1), rev((1:(series_length / 2 - 1))), 1e-5) # Calculate COI
    # For plotting COI ploygon:
    coi_x = c(axis_1[c(1, 1)] - dt * 0.5, axis_1, axis_1[c(series_length, series_length)] + dt * 0.5)
    logyint = axis_2[2] - axis_2[1]
    yl = c(log2(periods[scales_length]) + 0.5*logyint, log2(periods[1]) - 0.5*logyint)
    yr = rev(yl)
    coi_y = c(yl, log2(coi), yr)
    
    # Create list of outputs and return
    output <- list(scales=scales, periods=periods, freqs=1/periods, 
                   wave=wave, power=power,
                   axis_x=axis_1, axis_y=axis_2,
                   coi_x=coi_x, coi_y=coi_y)
    
    return(output)
    
}