
########## CARGA DE LOS DATOS ##########

setwd('D:/AGN/data/final_plots') # directorio de trabajo
sdss_xmatch <- read.csv('final_xmatch.csv', header = TRUE) # Archivo con el cross-match entre SDSS y los demás catálogos
len <- length(sdss_xmatch$X) # longitud de los datos

SFR <- sdss_xmatch$sfr_tot_p50_1 # SFR (en logarítmico)
M <- sdss_xmatch$lgm_tot_p50 # Masa estelar de la galaxia (en logarítmico)

L_Soft <- sdss_xmatch$Soft_luminosity_.erg.s. # luminosidad en rayos X (baja energía, 0.2-2 keV)
L_Soft <- L_Soft*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_Soft <- log10(L_Soft) # pasamos a escala logarítmica
Infinites_soft <- which(is.infinite(L_Soft) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_Soft[Infinites_soft] <- NA # y los eliminamos

L_Hard <- sdss_xmatch$Hard_luminosity_.erg.s. # luminosidad en rayos X (alta energía, 2-12 keV) 
L_Hard <- L_Hard*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_Hard <- log10(L_Hard) # pasamos a escala logarítmica
Infinites_hard <- which(is.infinite(L_Hard) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_Hard[Infinites_hard] <- NA # y los eliminamos


L_Radio <- sdss_xmatch$Radio_luminosity..erg.s.Hz. # luminosidad en Radio
L_Radio <- L_Radio*(1.4*10^9) # a erg/s
L_Radio <- L_Radio*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_Radio <- log10(L_Radio) # pasamos a escala logarítmica
Infinites_Radio <- which(is.infinite(L_Radio) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_Radio[Infinites_Radio] <- NA # y los eliminamos

L_IR <- sdss_xmatch$IR_luminosity..erg.s.Hz. # luminosidad en IR (banda W4) 
L_IR <- L_IR*((3*10^8)/(22.1*10^-6)) # a erg/s
L_IR <- L_IR*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_IR <- log10(L_IR) # pasamos a escala logarítmica
Infinites_IR4 <- which(is.infinite(L_IR) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_IR[Infinites_IR4] <- NA # y los eliminamos

L_IR1 <- sdss_xmatch$IR1_luminosity_.erg.s.Hz. # luminosidad en IR (banda W1) 
L_IR1 <- L_IR1*((3*10^8)/(3.4*10^-6)) # a erg/s
L_IR1 <- L_IR1*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_IR1 <- log10(L_IR1) # pasamos a escala logarítmica
Infinites_IR1 <- which(is.infinite(L_IR1) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_IR1[Infinites_IR1] <- NA # y los eliminamos

L_IR2 <- sdss_xmatch$IR2_luminosity_.erg.s.Hz. # luminosidad en IR (banda W2) 
L_IR2 <- L_IR2*((3*10^8)/(4.6*10^-6)) # a erg/s
L_IR2 <- L_IR2*(100/70)^2 # en su momento tomamos la cte de hubble naive (100 km/s/Mpc) y con esto pasamos a 70 km/s/Mpc
L_IR2 <- log10(L_IR2) # pasamos a escala logarítmica
Infinites_IR2 <- which(is.infinite(L_IR2) == TRUE) # identificamos los posibles valores infinitos en la luminosidad
L_IR2[Infinites_IR2] <- NA # y los eliminamos

W1 <- sdss_xmatch$W1mag
W2 <- sdss_xmatch$W2mag
color_IR <- W1 - W2
color_W2_W1 <- L_IR2 - L_IR1
q_24 <- log10(sdss_xmatch$IR_flux..erg.s.cm.2.Hz./(sdss_xmatch$S1_4*(10^-26))) # el 10^-26 es el factor de conversión para mJy
color_Radio_IR <- L_Radio - L_IR
color_X <- L_Hard - L_Soft
HR <- (sdss_xmatch$Hard_flux/7 - sdss_xmatch$Soft_flux/1)/(sdss_xmatch$Hard_flux/7 + sdss_xmatch$Soft_flux/1)
# como tenemos flujo y no número de cuentas cogemos la expresión tradicional para el hardness ratio
# y dividimos cada flujo por la energía media de esa banda (7 keV en el caso Hard, 1 keV en el caso Soft)


############ LIMPIEZA DE LOS DATOS Y FILTROS EN SEÑAL-RUIDO #############

for (i in 1:len){ 
  if (SFR[i] == -9999){ 
    SFR[i] <- NA 
  }
  if (M[i] == -9999){
    M[i] <- NA
  }
}

sigma_SoftX <- sdss_xmatch$Soft_flux_error # Noise
SN_SoftX <- sdss_xmatch$Soft_flux/sigma_SoftX # Signal-to-Noise ratio

sigma_HardX <- sdss_xmatch$Hard_flux_error # Noise
SN_HardX <- sdss_xmatch$Hard_flux/sigma_HardX # Signal-to-Noise ratio

sigma_Radio <- sdss_xmatch$e_S1_4 # Noise
SN_Radio <- sdss_xmatch$S1_4/sigma_Radio # Signal-to-Noise ratio

sigma_IR4 <- sdss_xmatch$IR_flux_error..erg.s.cm.2.Hz. # Noise
SN_IR4 <- sdss_xmatch$IR_flux..erg.s.cm.2.Hz./sigma_IR4 # Signal-to-Noise ratio

sigma_IR1 <- sdss_xmatch$IR1_flux_error_.erg.s.cm2.Hz. # Noise
SN_IR1 <- sdss_xmatch$IR1_flux_.erg.s.cm2.Hz./sigma_IR1 # Signal-to-Noise ratio

sigma_IR2 <- sdss_xmatch$IR2_flux_error_.erg.s.cm2.Hz. # Noise
SN_IR2 <- sdss_xmatch$IR2_flux_.erg.s.cm2.Hz./sigma_IR2 # Signal-to-Noise ratio

# Filtro en flujo en rayos X blandos:
for (i in 1:len){ 
  if ((is.na(SN_SoftX[i]) == FALSE) & (SN_SoftX[i] < 2)){ 
    L_Soft[i] <- NA
  }
}

# Filtro en flujo en rayos X duros:
for (i in 1:len){ 
  if ((is.na(SN_HardX[i]) == FALSE) & (SN_HardX[i] < 2)){ 
    L_Hard[i] <- NA
  }
}

# Filtro en flujo en Radio:
for (i in 1:len){ 
  if ((is.na(SN_Radio[i]) == FALSE) & (SN_Radio[i] < 2)){ 
    L_Radio[i] <- NA
  }
}

# Filtro en flujo en IR (banda W4):
for (i in 1:len){ 
  if ((is.na(SN_IR4[i]) == FALSE) & (SN_IR4[i] < 2)){ 
    L_IR[i] <- NA
  }
}

# Filtro en flujo en IR (banda W1):
for (i in 1:len){ 
  if ((is.na(SN_IR1[i]) == FALSE) & (SN_IR1[i] < 2)){ 
    L_IR1[i] <- NA
  }
}

# Filtro en flujo en IR (banda W2):
for (i in 1:len){ 
  if ((is.na(SN_IR2[i]) == FALSE) & (SN_IR2[i] < 2)){ 
    L_IR2[i] <- NA
  }
}


# Filtro en flujo en colorW2_W1:

for (i in 1:len){ 
  if (((is.na(SN_IR1[i]) == FALSE) & (SN_IR1[i] < 2)) | ((is.na(SN_IR2[i]) == FALSE) & (SN_IR2[i] < 2))){ 
    color_W2_W1[i] <- NA
    color_IR[i] <- NA
  }
}

# Filtro en flujo en color_Radio_IR:

for (i in 1:len){ 
  if (((is.na(SN_IR4[i]) == FALSE) & (SN_IR4[i] < 2)) | ((is.na(SN_Radio[i]) == FALSE) & (SN_Radio[i] < 2))){ 
    color_Radio_IR[i] <- NA
    q_24[i] <- NA
  }
}

# Filtro en flujo en color_X:

for (i in 1:len){ 
  if (((is.na(SN_SoftX[i]) == FALSE) & (SN_SoftX[i] < 2)) | ((is.na(SN_HardX[i]) == FALSE) & (SN_HardX[i] < 2))){ 
    color_X[i] <- NA
  }
}


############# GRIDS Y FUNCIONES PARA EL CÁLCULO DE LOS HISTOGRAMAS 2D #############

SFR_y <- seq(min(na.omit(SFR)), max(na.omit(SFR)), (max(na.omit(SFR)) - min(na.omit(SFR)))/249)
SFR_grid <- matrix(SFR_y, nrow = 250, ncol = 250, byrow = TRUE)
d_SFR = SFR_y[2] - SFR_y[1]

M_X <- seq(min(na.omit(M)), max(na.omit(M)), (max(na.omit(M)) - min(na.omit(M)))/249)
M_grid <- matrix(M_X, nrow = 250, ncol = 250, byrow = FALSE)
d_M = M_X[2] - M_X[1]
grid_len <- length(SFR_grid)

densidad <- function(SFR_i, M_j, SFR_p, M_p, Delta_SFR_p, Delta_M_p){
  dangerous_terms <- c(SFR_p, M_p, Delta_SFR_p, Delta_M_p)
  term_1 = (M_j - M_p)/Delta_M_p
  term_2 = (SFR_i - SFR_p)/Delta_SFR_p
  term = -term_1^2 - term_2^2
  term = term/2
  term = exp(term)
  term = term/(2*pi*Delta_SFR_p*Delta_M_p)
  if (length(which(is.na(dangerous_terms) == TRUE)) > 0){
    term = matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  }
  return(term)
} 


densidad_IR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_Radio <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_Soft <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_Hard <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_W1 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_W2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_colorIR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_color_W2_W1 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_q24 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_color_Radio_IR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_colorX <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_IR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Radio <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Soft <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Hard <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_W1 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_W2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_colorIR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_color_W2_W1 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_q24 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_color_Radio_IR <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_colorX <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_IR_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Radio_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Soft_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_Hard_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_W1_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_W2_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_colorIR_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_color_W2_W1_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_q24_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_color_Radio_IR_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
densidad_pesada_colorX_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
suavizado_uniforme <- 0.3

for (p in 1:len){
  densidad_p = densidad(SFR_i = SFR_grid, M_j = M_grid, SFR_p = SFR[p], M_p = M[p],
                          Delta_SFR_p = suavizado_uniforme, Delta_M_p = suavizado_uniforme)
  L_p = L_IR[p]
  if (is.na(L_p) == FALSE){
    densidad_IR = densidad_IR + densidad_p
    densidad_pesada_IR = densidad_pesada_IR + L_p*densidad_p
    densidad_pesada_IR_2 = densidad_pesada_IR_2 + L_p^2*densidad_p
  }
  L_p = L_Radio[p]
  if (is.na(L_p) == FALSE){
    densidad_Radio = densidad_Radio + densidad_p
    densidad_pesada_Radio = densidad_pesada_Radio + L_p*densidad_p
    densidad_pesada_Radio_2 = densidad_pesada_Radio_2 + L_p^2*densidad_p
  }
  L_p = L_Soft[p]
  if (is.na(L_p) == FALSE){
    densidad_Soft = densidad_Soft + densidad_p
    densidad_pesada_Soft = densidad_pesada_Soft + L_p*densidad_p
    densidad_pesada_Soft_2 = densidad_pesada_Soft_2 + L_p^2*densidad_p
  }
  L_p = L_Hard[p]
  if (is.na(L_p) == FALSE){
    densidad_Hard = densidad_Hard + densidad_p
    densidad_pesada_Hard = densidad_pesada_Hard + L_p*densidad_p
    densidad_pesada_Hard_2 = densidad_pesada_Hard_2 + L_p^2*densidad_p
  }
  L_p = L_IR1[p]
  if (is.na(L_p) == FALSE){
    densidad_W1 = densidad_W1 + densidad_p
    densidad_pesada_W1 = densidad_pesada_W1 + L_p*densidad_p
    densidad_pesada_W1_2 = densidad_pesada_W1_2 + L_p^2*densidad_p
  }
  L_p = L_IR2[p]
  if (is.na(L_p) == FALSE){
    densidad_W2 = densidad_W2 + densidad_p
    densidad_pesada_W2 = densidad_pesada_W2 + L_p*densidad_p
    densidad_pesada_W2_2 = densidad_pesada_W2_2 + L_p^2*densidad_p
  }
  L_p = color_IR[p]
  if (is.na(L_p) == FALSE){
    densidad_colorIR = densidad_colorIR + densidad_p
    densidad_pesada_colorIR = densidad_pesada_colorIR + L_p*densidad_p
    densidad_pesada_colorIR_2 = densidad_pesada_colorIR_2 + L_p^2*densidad_p
  }
  L_p = color_W2_W1[p]
  if (is.na(L_p) == FALSE){
    densidad_color_W2_W1 = densidad_color_W2_W1 + densidad_p
    densidad_pesada_color_W2_W1 = densidad_pesada_color_W2_W1 + L_p*densidad_p
    densidad_pesada_color_W2_W1_2 = densidad_pesada_color_W2_W1_2 + L_p^2*densidad_p
  }
  L_p = q_24[p]
  if (is.na(L_p) == FALSE){
    densidad_q24 = densidad_q24 + densidad_p
    densidad_pesada_q24 = densidad_pesada_q24 + L_p*densidad_p
    densidad_pesada_q24_2 = densidad_pesada_q24_2 + L_p^2*densidad_p
  }
  L_p = color_Radio_IR[p]
  if (is.na(L_p) == FALSE){
    densidad_color_Radio_IR = densidad_color_Radio_IR + densidad_p
    densidad_pesada_color_Radio_IR = densidad_pesada_color_Radio_IR + L_p*densidad_p
    densidad_pesada_color_Radio_IR_2 = densidad_pesada_color_Radio_IR_2 + L_p^2*densidad_p
  }
  L_p = color_X[p]
  if (is.na(L_p) == FALSE){
    densidad_colorX = densidad_colorX + densidad_p
    densidad_pesada_colorX = densidad_pesada_colorX + L_p*densidad_p
    densidad_pesada_colorX_2 = densidad_pesada_colorX_2 + L_p^2*densidad_p
  }
}

###################### GETTING ISOCONTOUR LEVELS FOR THE MASS ########################

library(pracma)
library(colorRamps)
library(fields)
library(scales)

get_isocontour_levels = function(Density, dx, dy){
  
  mass_histogram = Density*dx*dy
  sorted_flat = sort.int(as.vector(Density), index.return = TRUE)$ix
  density_sorted = as.vector(Density)[sorted_flat]
  cumulative_mass = cumsum(as.vector(mass_histogram)[sorted_flat])
  plot(density_sorted, cumulative_mass, type = "l", col = "blue",
       xlab = "density", ylab = "cummulative mass")
  
  fraction_sorted = cumulative_mass/cumulative_mass[length(cumulative_mass)]
  fraction_interp = approxfun(x = density_sorted, y = fraction_sorted)
  fraction = matrix(fraction_interp(Density), nrow = nrow(Density), ncol = ncol(Density))
  
  percentil_10 = quantile(na.omit(fraction), probs = 0.10)
  percentil_50 = quantile(na.omit(fraction), probs = 0.50)
  percentil_90 = quantile(na.omit(fraction), probs = 0.90)
  
  image.plot(M_X, SFR_y, Density, col = matlab.like(1000), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
             xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density')
  contour(M_X, SFR_y, fraction, levels = percentil_10, add = TRUE, lty = 2)
  contour(M_X, SFR_y, fraction, levels = percentil_50, add = TRUE, lty = 3)
  contour(M_X, SFR_y, fraction, levels = percentil_90, add = TRUE, lty = 1)

  image.plot(M_X, SFR_y, fraction, col = matlab.like(1000), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
             xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density')
  contour(M_X, SFR_y, fraction, levels = percentil_10, add = TRUE, lty = 2)
  contour(M_X, SFR_y, fraction, levels = percentil_50, add = TRUE, lty = 3)
  contour(M_X, SFR_y, fraction, levels = percentil_90, add = TRUE, lty = 1)
  
  return(fraction)
}

get_IR1 = get_isocontour_levels(densidad_W1, d_M, d_SFR)
get_IR2 = get_isocontour_levels(densidad_W2, d_M, d_SFR)
get_IR4 = get_isocontour_levels(densidad_IR, d_M, d_SFR)
get_Radio = get_isocontour_levels(densidad_Radio, d_M, d_SFR)
get_Soft = get_isocontour_levels(densidad_Soft, d_M, d_SFR)
get_Hard = get_isocontour_levels(densidad_Hard, d_M, d_SFR)

get_color_w2_w1 = get_isocontour_levels(densidad_color_W2_W1, d_M, d_SFR)
get_color_Radio_IR = get_isocontour_levels(densidad_color_Radio_IR, d_M, d_SFR)
get_color_X = get_isocontour_levels(densidad_colorX, d_M, d_SFR)



############ CÁLCULO DE LoS <L> y <COLOR> Y DE SUS CORRESPONDIENTES DISPERSIONES #############

L_mean_IR = densidad_pesada_IR/densidad_IR
L_mean_Radio = densidad_pesada_Radio/densidad_Radio
L_mean_Soft = densidad_pesada_Soft/densidad_Soft
L_mean_Hard = densidad_pesada_Hard/densidad_Hard
L_mean_IR1 = densidad_pesada_W1/densidad_W1
L_mean_IR2 = densidad_pesada_W2/densidad_W2
colorIR_mean <- densidad_pesada_colorIR/densidad_colorIR
color_W2_W1_mean <- densidad_pesada_color_W2_W1/densidad_color_W2_W1
q24_mean <- densidad_pesada_q24/densidad_q24
color_Radio_IR_mean <- densidad_pesada_color_Radio_IR/densidad_color_Radio_IR
colorX_mean <- densidad_pesada_colorX/densidad_colorX

L_mean_squared_IR = densidad_pesada_IR_2/densidad_IR
L_mean_squared_Radio = densidad_pesada_Radio_2/densidad_Radio
L_mean_squared_Soft = densidad_pesada_Soft_2/densidad_Soft
L_mean_squared_Hard = densidad_pesada_Hard_2/densidad_Hard
L_mean_squared_IR1 <- densidad_pesada_W1_2/densidad_W1
L_mean_squared_IR2 <- densidad_pesada_W2_2/densidad_W2
colorIR_mean_squared <- densidad_pesada_colorIR_2/densidad_colorIR
color_W2_W1_mean_squared <- densidad_pesada_color_W2_W1_2/densidad_color_W2_W1
q24_mean_squared <- densidad_pesada_q24_2/densidad_q24
color_Radio_IR_mean_squared <- densidad_pesada_color_Radio_IR_2/densidad_color_Radio_IR
colorX_mean_squared <- densidad_pesada_colorX_2/densidad_colorX

Dispersion_IR <- sqrt(L_mean_squared_IR - L_mean_IR^2)
Dispersion_Radio <- sqrt(L_mean_squared_Radio - L_mean_Radio^2)
Dispersion_Soft <- sqrt(L_mean_squared_Soft - L_mean_Soft^2)
Dispersion_Hard <- sqrt(L_mean_squared_Hard - L_mean_Hard^2)
Dispersion_IR1 <- sqrt(L_mean_squared_IR1 - L_mean_IR1^2)
Dispersion_IR2 <- sqrt(L_mean_squared_IR2 - L_mean_IR2^2)
Dispersion_colorIR <- sqrt(colorIR_mean_squared - colorIR_mean^2)
Dispersion_color_W2_W1 <- sqrt(color_W2_W1_mean_squared - color_W2_W1_mean^2)
Dispersion_q24 <- sqrt(q24_mean_squared - q24_mean^2)
Dispersion_color_Radio_IR <- sqrt(color_Radio_IR_mean_squared - color_Radio_IR_mean^2)
Dispersion_colorX <- sqrt(colorX_mean_squared - colorX_mean^2)



############# HOMOGENEIZAR LOS RANGOS DINÁMICOS DE LOS HISTOGRAMAS 2D ################

#  Para que se vean mejor las figuras, acotamos el rango dinámico de todas las luminosidades a sus respectivos percentiles 15 y 85:

L_mean_IR_truncated <- L_mean_IR
L_mean_IR_truncated[L_mean_IR_truncated > quantile(L_mean_IR, probs = 0.85)] <- quantile(L_mean_IR, probs = 0.85)
L_mean_IR_truncated[L_mean_IR_truncated < quantile(L_mean_IR, probs = 0.15)] <- quantile(L_mean_IR, probs = 0.15)

L_mean_IR1_truncated <- L_mean_IR1
L_mean_IR1_truncated[L_mean_IR1_truncated > quantile(L_mean_IR1, probs = 0.85)] <- quantile(L_mean_IR1, probs = 0.85)
L_mean_IR1_truncated[L_mean_IR1_truncated < quantile(L_mean_IR1, probs = 0.15)] <- quantile(L_mean_IR1, probs = 0.15)

L_mean_IR2_truncated <- L_mean_IR2
L_mean_IR2_truncated[L_mean_IR2_truncated > quantile(L_mean_IR2, probs = 0.85)] <- quantile(L_mean_IR2, probs = 0.85)
L_mean_IR2_truncated[L_mean_IR2_truncated < quantile(L_mean_IR2, probs = 0.15)] <- quantile(L_mean_IR2, probs = 0.15)

L_mean_Radio_truncated <- L_mean_Radio
L_mean_Radio_truncated[L_mean_Radio_truncated > quantile(L_mean_Radio, probs = 0.85)] <- quantile(L_mean_Radio, probs = 0.85)
L_mean_Radio_truncated[L_mean_Radio_truncated < quantile(L_mean_Radio, probs = 0.15)] <- quantile(L_mean_Radio, probs = 0.15)

L_mean_Soft_truncated <- L_mean_Soft
L_mean_Soft_truncated[L_mean_Soft_truncated > quantile(L_mean_Soft, probs = 0.85)] <- quantile(L_mean_Soft, probs = 0.85)
L_mean_Soft_truncated[L_mean_Soft_truncated < quantile(L_mean_Soft, probs = 0.15)] <- quantile(L_mean_Soft, probs = 0.15)

L_mean_Hard_truncated <- L_mean_Hard
L_mean_Hard_truncated[L_mean_Hard_truncated > quantile(L_mean_Hard, probs = 0.85)] <- quantile(L_mean_Hard, probs = 0.85)
L_mean_Hard_truncated[L_mean_Hard_truncated < quantile(L_mean_Hard, probs = 0.15)] <- quantile(L_mean_Hard, probs = 0.15)


# Ídem para los colores:

colorIR_truncated <- colorIR_mean
colorIR_truncated[colorIR_truncated > quantile(colorIR_mean, probs = 0.85)] <- quantile(colorIR_mean, probs = 0.85)
colorIR_truncated[colorIR_truncated < quantile(colorIR_mean, probs = 0.15)] <- quantile(colorIR_mean, probs = 0.15)

color_W2_W1_truncated <- color_W2_W1_mean
color_W2_W1_truncated[color_W2_W1_truncated > quantile(color_W2_W1_mean, probs = 0.85)] <- quantile(color_W2_W1_mean, probs = 0.85)
color_W2_W1_truncated[color_W2_W1_truncated < quantile(color_W2_W1_mean, probs = 0.15)] <- quantile(color_W2_W1_mean, probs = 0.15)

color_Radio_IR_truncated <- color_Radio_IR_mean
color_Radio_IR_truncated[color_Radio_IR_truncated > quantile(color_Radio_IR_mean, probs = 0.85)] <- quantile(color_Radio_IR_mean, probs = 0.85)
color_Radio_IR_truncated[color_Radio_IR_truncated < quantile(color_Radio_IR_mean, probs = 0.15)] <- quantile(color_Radio_IR_mean, probs = 0.15)

colorX_truncated <- colorX_mean
colorX_truncated[colorX_truncated > quantile(colorX_mean, probs = 0.85)] <- quantile(colorX_mean, probs = 0.85)
colorX_truncated[colorX_truncated < quantile(colorX_mean, probs = 0.15)] <- quantile(colorX_mean, probs = 0.15)


# Para las dispersiones acotamos el rango dinámico al intervalo (0 - rango intercuartílico de la magnitud considerada):

Dispersion_IR_truncated <- Dispersion_IR
limit = (quantile(na.omit(L_IR), probs = 0.85) - quantile(na.omit(L_IR), probs = 0.15))
Dispersion_IR_truncated[Dispersion_IR_truncated > limit] = limit
Dispersion_IR_truncated[Dispersion_IR_truncated < 0] = 0

Dispersion_Radio_truncated <- Dispersion_Radio
limit = (quantile(na.omit(L_Radio), probs = 0.85) - quantile(na.omit(L_Radio), probs = 0.15))
Dispersion_Radio_truncated[Dispersion_Radio_truncated > limit] = limit
Dispersion_Radio_truncated[Dispersion_Radio_truncated < 0] = 0

Dispersion_Soft_truncated <- Dispersion_Soft
limit = (quantile(na.omit(L_Soft), probs = 0.85) - quantile(na.omit(L_Soft), probs = 0.15))
Dispersion_Soft_truncated[Dispersion_Soft_truncated > limit] = limit
Dispersion_Soft_truncated[Dispersion_Soft_truncated < 0] = 0

Dispersion_Hard_truncated <- Dispersion_Hard
limit = (quantile(na.omit(L_Hard), probs = 0.85) - quantile(na.omit(L_Hard), probs = 0.15))
Dispersion_Hard_truncated[Dispersion_Hard_truncated > limit] = limit
Dispersion_Hard_truncated[Dispersion_Hard_truncated < 0] = 0

Dispersion_IR1_truncated <- Dispersion_IR1
limit = (quantile(na.omit(L_IR1), probs = 0.85) - quantile(na.omit(L_IR1), probs = 0.15))
Dispersion_IR1_truncated[Dispersion_IR1_truncated > limit] = limit
Dispersion_IR1_truncated[Dispersion_IR1_truncated < 0] = 0

Dispersion_IR2_truncated <- Dispersion_IR2
limit = (quantile(na.omit(L_IR2), probs = 0.85) - quantile(na.omit(L_IR2), probs = 0.15))
Dispersion_IR2_truncated[Dispersion_IR2_truncated > limit] = limit
Dispersion_IR2_truncated[Dispersion_IR2_truncated < 0] = 0




Dispersion_color_W2_W1_truncated <- Dispersion_color_W2_W1
limit = (quantile(na.omit(color_W2_W1), probs = 0.85) - quantile(na.omit(color_W2_W1), probs = 0.15))
Dispersion_color_W2_W1_truncated[Dispersion_color_W2_W1_truncated > limit] = limit
Dispersion_color_W2_W1_truncated[Dispersion_color_W2_W1_truncated < 0] = 0

Dispersion_color_Radio_IR_truncated <- Dispersion_color_Radio_IR
limit = (quantile(na.omit(color_Radio_IR), probs = 0.85) - quantile(na.omit(color_Radio_IR), probs = 0.15))
Dispersion_color_Radio_IR_truncated[Dispersion_color_Radio_IR_truncated > limit] = limit
Dispersion_color_Radio_IR_truncated[Dispersion_color_Radio_IR_truncated < 0] = 0

Dispersion_colorX_truncated <- Dispersion_colorX
limit = (quantile(na.omit(color_X), probs = 0.85) - quantile(na.omit(color_X), probs = 0.15))
Dispersion_colorX_truncated[Dispersion_colorX_truncated > limit] = limit
Dispersion_colorX_truncated[Dispersion_colorX_truncated < 0] = 0



############ HISTOGRAMAS 2D DE LOS NUEVOS <L> Y <COLOR> Y DE SUS RESPECTIVAS DISPERSIONES ################


image.plot(M_X, SFR_y, L_mean_IR_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W4 band')
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, L_mean_IR1_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W1 band')
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, L_mean_IR2_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W2 band')
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, L_mean_Radio_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Radio')
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, L_mean_Soft_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Soft X-Rays')
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, L_mean_Hard_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Hard X-Rays')
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, color_W2_W1_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <IR color>')
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, color_Radio_IR_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <Radio-IR(W4) color>')
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, colorX_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <X-Rays color>')
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .95), add = TRUE, lty = 1)



image.plot(M_X, SFR_y, Dispersion_IR_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W4 band)')
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR4, levels = quantile(get_IR4, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_Radio_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio')
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Radio, levels = quantile(get_Radio, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_Soft_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Soft X-Rays')
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Soft, levels = quantile(get_Soft, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_Hard_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Hard X-Rays')
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_Hard, levels = quantile(get_Hard, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_IR1_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W1 band)')
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR1, levels = quantile(get_IR1, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_IR2_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W2 band)')
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_IR2, levels = quantile(get_IR2, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_color_W2_W1_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR color')
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_w2_w1, levels = quantile(get_color_w2_w1, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_color_Radio_IR_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio-IR(W4) color')
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_Radio_IR, levels = quantile(get_color_Radio_IR, probs = .95), add = TRUE, lty = 1)

image.plot(M_X, SFR_y, Dispersion_colorX_truncated, col = matlab.like(15), ylim = c(-3, 1.5), xlim = c(8.5, 11.5),
           xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion X-Rays color')
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .50), add = TRUE, lty = 3)
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .80), add = TRUE, lty = 2)
contour(M_X, SFR_y, get_color_X, levels = quantile(get_color_X, probs = .95), add = TRUE, lty = 1)



############ CÁLCULO DE LOS LUMINOSIDADES Y COLORES MEDIOS Y DE SUS EXCESOS ############

L_media_soft <- rep(0, len)
L_media_hard <- rep(0, len)
L_media_radio <- rep(0, len)
L_media_ir <- rep(0, len)
L_media_ir1 <- rep(0, len)
L_media_ir2 <- rep(0, len)
color_W2_W1_media <- rep(0, len)
color_Radio_IR_media <- rep(0, len)
colorX_media <- rep(0, len)

binsize_SFR <- (SFR_grid[1,2] - SFR_grid[1,1])/2
binsize_M <- (M_grid[2,1] - M_grid[1,1])/2


L_mean_vector <- function(L, L_media, L_mean){
  for (i in 1:len){
    if ((is.na(L[i]) == FALSE) & (is.na(SFR[i]) == FALSE) & (is.na(M[i]) == FALSE)){
      fila = 0
      columna = 0
      for (j in 1:length(M_X)){
        if (abs(M[i] - M_grid[j,1]) < binsize_M){
          fila = j
          break
        }
      }
      for (k in 1:length(SFR_y)){
        if (abs(SFR[i] - SFR_grid[1,k]) < binsize_SFR){
          columna = k
          break
        }
      }
      if ((fila != 0) & (columna != 0)){
        L_media[i] <- L_mean[fila, columna]
      }
      else{
        L_media[i] <- NA
      }
    }
    else{
      L_media[i] <- NA
    }
  }
  return(L_media)
}


L_media_soft <- L_mean_vector(L = L_Soft, L_media = L_media_soft, L_mean = L_mean_Soft)
L_media_hard <- L_mean_vector(L = L_Hard, L_media = L_media_hard, L_mean = L_mean_Hard)
L_media_radio <- L_mean_vector(L = L_Radio, L_media = L_media_radio, L_mean = L_mean_Radio)
L_media_ir <- L_mean_vector(L = L_IR, L_media = L_media_ir, L_mean = L_mean_IR)
L_media_ir1 <- L_mean_vector(L = L_IR1, L_media = L_media_ir1, L_mean = L_mean_IR1)
L_media_ir2 <- L_mean_vector(L = L_IR2, L_media = L_media_ir2, L_mean = L_mean_IR2)
color_W2_W1_media <- L_mean_vector(L = color_W2_W1, L_media = color_W2_W1_media, L_mean = color_W2_W1_mean)
color_Radio_IR_media <- L_mean_vector(L = color_Radio_IR, L_media = color_Radio_IR_media, L_mean = color_Radio_IR_mean)
colorX_media <- L_mean_vector(L = color_X, L_media = colorX_media, L_mean = colorX_mean)

L_reducida_Soft <- L_Soft - L_media_soft
L_reducida_Hard <- L_Hard - L_media_hard
L_reducida_Radio <- L_Radio - L_media_radio
L_reducida_IR <- L_IR - L_media_ir
L_reducida_IR1 <- L_IR1 - L_media_ir1
L_reducida_IR2 <- L_IR2 - L_media_ir2
color_W2_W1_reducido <- color_W2_W1 - color_W2_W1_media
color_Radio_IR_reducido <- color_Radio_IR - color_Radio_IR_media
color_X_reducido <- color_X - colorX_media



############## CÁLCULO DE LOS PARÁMETROS DE CONVERSIÓN ENTRE NUESTROS COLORES Y LOS DE LA LITERATURA ###############

W1_W2_norm <- log10((171*3.4*10^-6)/(309*4.6*10^-6)) # cociente entre flujos cero de WISE
nuradio_nu24 <- log10((1.4*10^9)/((3*10^8)/(22.1*10^-6))) # cociente de frecuencias (1.4 GHz/22.1 micras)

################################## ACTIVENESS PLOTS #####################################

sigmas = function(lum){
  
  x_ordenado = na.omit(sort.int(lum, na.last = TRUE, index.return = TRUE)$x)
  ranking = c(1:length(x_ordenado))
  ranking_interp = approxfun(x_ordenado, ranking)
  N = length(x_ordenado)
  z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
  N_greater = N - ranking_interp(z)
  N_normal = rep(0, length(N_greater))
  pivot_point = median(x_ordenado)
  # pivot_point = moda_color_W2_W1_red
  for (i in 1:length(z)){
    if (z[i] <= pivot_point){
      N_symmetric = 2*ranking_interp(pivot_point) - ranking_interp(z[i])
    }
    else {
      N_symmetric = ranking_interp(pivot_point - z[i])
    }
    N_normal[i] = min(N_symmetric, N_greater[i])
  }
  
  N_active = N_greater - N_normal
  f_active = N_active/N_greater
  z_active = z[min(which(f_active >= 0.5))]
  n_sigmas = (z_active - quantile(na.omit(lum), probs = 0.50))/
    (quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16))
  return(as.numeric(n_sigmas))
}


sigma_minima = min(sigmas(L_reducida_IR2), sigmas(L_reducida_Radio), sigmas(L_reducida_Hard),
                   sigmas(color_W2_W1_reducido), sigmas(color_Radio_IR_reducido), sigmas(color_X_reducido))

sigma_maxima = max(sigmas(L_reducida_IR2), sigmas(L_reducida_Radio), sigmas(L_reducida_Hard),
                   sigmas(color_W2_W1_reducido), sigmas(color_Radio_IR_reducido), sigmas(color_X_reducido))


activeness = function(lum, xlabel, colour){
  
  x_ordenado = na.omit(sort.int(lum, na.last = TRUE, index.return = TRUE)$x)
  ranking = c(1:length(x_ordenado))
  ranking_interp = approxfun(x_ordenado, ranking)
  N = length(x_ordenado)
  z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
  N_greater = N - ranking_interp(z)
  N_normal = rep(0, length(N_greater))
  pivot_point = median(x_ordenado)
  # pivot_point = moda_color_W2_W1_red
  for (i in 1:length(z)){
    if (z[i] <= pivot_point){
      N_symmetric = 2*ranking_interp(pivot_point) - ranking_interp(z[i])
    }
    else {
      N_symmetric = ranking_interp(pivot_point - z[i])
    }
    N_normal[i] = min(N_symmetric, N_greater[i])
  }
  
  N_active = N_greater - N_normal
  f_active = N_active/N_greater
  z_active = z[min(which(f_active >= 0.5))]
  minimo = sigma_minima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
    quantile(na.omit(lum), probs = 0.50)
  maximo = sigma_maxima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
    quantile(na.omit(lum), probs = 0.50)
  
  library(ggplot2)
  df = data.frame(z, log10(N_normal))
  df2 = data.frame(z, log10(N_greater))
  print(ggplot(df, aes(x = df$z)) +
          ylim(0, max(na.omit(log10(N_greater)))) +
          geom_line(aes(y = df$log10.N_normal.), color = colour, linetype = 'dashed') +
          geom_line(aes(y = df2$log10.N_greater.), color = 'black') +
          geom_rect(aes(xmin = minimo, xmax = maximo, ymin = -Inf, ymax = Inf), color = 'grey', alpha = 0.005) +
          geom_vline(xintercept = z_active) +
          theme_bw() + 
          xlab(xlabel) + 
          ylab('log (N > x)'))
  return(z_active)
}

activeness(L_reducida_IR2, 'Luminosity excess (IR W2 band)', 'red')
activeness(L_reducida_Radio, 'Luminosity excess (Radio)', 'green')
activeness(L_reducida_Hard, 'Luminosity excess (Hard X-Rays)', 'blue')

activeness(color_W2_W1_reducido, 'IR color excess', 'red')
activeness(color_Radio_IR_reducido, 'Radio-IR color excess', 'green')
activeness(color_X_reducido, 'X-Rays color excess', 'blue')

############## CÁLCULO DE LOS THRESHOLDS PARA CADA LUMINOSIDAD Y COLOR ################

threshold_colorIR_lit = 0.4*0.8 + W1_W2_norm
threshold_color_Radio_IR_lit = 0.23 + nuradio_nu24
threshold_colorX_lit = log10(7)

threshold_ir2_excess = as.numeric(activeness(L_reducida_IR2, 'Luminosity excess (IR W2 band)', 'red'))
threshold_colorIR_excess = as.numeric(activeness(color_W2_W1_reducido, 'IR color excess', 'red'))
threshold_radio_excess = as.numeric(activeness(L_reducida_Radio, 'Luminosity excess (Radio)', 'green'))
threshold_color_Radio_IR_excess = as.numeric(activeness(color_Radio_IR_reducido, 'Radio-IR color excess', 'green'))
threshold_hard_excess = as.numeric(activeness(L_reducida_Hard, 'Luminosity excess (Hard X-Rays)', 'blue'))
threshold_colorX_excess = as.numeric(activeness(color_X_reducido, 'X-Rays color excess', 'blue'))


################################### FINAL PLOTS ########################################

# Función que dados los excesos de color y luminosidad de una galaxia le asigna
# un color de acuerdo con nuestra clasificación:
# Rojo: active galaxy. Verde: luminosity extreme. Amarillo: color extreme. Azul: normal galaxy
# get_galaxy_type hace esto para un objeto dados sus excesos de luminosidad y color

get_galaxy_type = function(L_i, color_i, threshold_L, threshold_color, alpha_value){
  galaxy_type = alpha("blue", alpha_value/10)
  if ((is.na(L_i) == FALSE) & (L_i > threshold_L)){
    galaxy_type = alpha("green", alpha_value)
    if ((is.na(color_i) == FALSE) & (color_i > threshold_color)){
      galaxy_type = alpha("red4", alpha_value)
    }
  }
  else if ((is.na(color_i) == FALSE) & (color_i > threshold_color)){
    galaxy_type = alpha("gold", alpha_value)
  }
  return (galaxy_type)
}

# galaxy_classification lo generaliza y recorre todos los objetos para cada banda

galaxy_classification = function(L, color, threshold_L, threshold_color, alpha_value){
  N = length(L)
  galaxy_class = as.vector(rep(alpha("black", 1), N))
  for (i in 1:N){
    galaxy_class[i] = get_galaxy_type(L[i], color[i], threshold_L, threshold_color, alpha_value)
  }
  return(galaxy_class)
}

# Pintamos primero los excesos: exceso de color vs exceso de luminosidad:

galaxy_class_IR = galaxy_classification(L_reducida_IR, color_W2_W1_reducido,
                                           threshold_ir2_excess, threshold_colorIR_excess, 0.7)

plot(L_reducida_IR2, color_W2_W1_reducido, col = galaxy_class_IR, pch = 16, ylim = c(-0.5, 1),
     xlab = "IR (W2 band) luminosity excess", ylab = "IR color excess")
abline(h = threshold_colorIR_excess, lty = 2, lwd = 2)
abline(v = threshold_ir2_excess, lty = 4, lwd = 2)

# Pintamos ahora color frente a luminosidad y comparamos con la literatura:

par(mar = c(5, 4, 4, 4) + 0.1) # abrimos una nueva ventana para la figura
plot(L_IR2, color_W2_W1, col = galaxy_class_IR, pch = 16, ylim = c(-0.75, 0.6),
     xlab = "", ylab = "")
abline(h = threshold_colorIR_lit, lty = 2, lwd = 2) # comparación con la literatura
axis(2, col = "black", las = 0)  # (las = 1 makes horizontal labels). 2 es por el eje: es el de ordenadas de la izquierda
mtext("IR color", side = 2, line = 2.5) # side = 2 sitúa la etiqueta en el eje de ordenadas de la izquierda
par(new = TRUE) # nuevo overplot para pintar el eje de ordenadas de la derecha para la comparación con la literatura
plot(L_IR2, color_IR, pch = 16, col = alpha('blue', 0.1), # pintamos por encima una figura vacía (type = 'n')
     main = '', 
     xlab = '', ylab = '', type = 'n', axes = FALSE, 
     ylim = c((-0.75 - W1_W2_norm)/0.4, (0.6 - W1_W2_norm)/0.4)) # conversión del ylim de la figura principal a ylim de W1-W2
mtext("W1-W2", side = 4, col = "black", line = 2) # texto en el eje de ordenadas de la derecha
axis(4, col = "black", col.axis = "black", las = 0) # 4 es por el eje: es el de ordenadas de la derecha
axis(1, pretty(range(L_IR2))) # el eje 1 (abcisas de abajo) cubre todo el rango de la luminosidad que queremos pintar
mtext("log[L (IR W2 band)]", side = 1, col = "black", line = 2.5) # ponemos la etiqueta a este eje



galaxy_class_Radio = galaxy_classification(L_reducida_Radio, color_Radio_IR_reducido,
                                           threshold_radio_excess, threshold_color_Radio_IR_excess, 0.3)

plot(L_reducida_Radio, color_Radio_IR_reducido, col = galaxy_class_Radio, pch = 16, xlim = c(-1.5, 3), ylim = c(-2, 3),
     xlab = "Radio luminosity excess", ylab = "Radio-IR color excess")
abline(h = threshold_color_Radio_IR_excess, lty = 2, lwd = 2)
abline(v = threshold_radio_excess, lty = 4, lwd = 2)


par(mar = c(5, 4, 4, 4) + 0.1)
plot(L_Radio, color_Radio_IR, col = galaxy_class_Radio, pch = 16,
     xlab = "", ylab = "")
abline(h = threshold_color_Radio_IR_lit, lty = 2, lwd = 2)
abline(v = 40, lty = 4, lwd = 2)
axis(2, col = "black", las = 0) 
mtext("Radio-IR color", side = 2, line = 2.5)
par(new = TRUE)
plot(L_Radio, q_24, pch = 16, col = alpha('blue', 0.1),
     main = '', 
     xlab = '', ylab = '', type = 'n', axes = FALSE,
     ylim = rev(range(na.omit(q_24)))) # el rango de q24 va al revés que color_Radio_IR (de mayor a menor) 
mtext("q24", side = 4, col = "black", line = 2) 
axis(4, col = "black", col.axis = "black", las = 0)
axis(1, pretty(range(L_Radio)))
mtext("log[L (Radio)]", side = 1, col = "black", line = 2.5)



galaxy_class_X = galaxy_classification(L_reducida_Hard, color_X_reducido,
                                           threshold_hard_excess, threshold_colorX_excess, 0.7)

plot(L_reducida_Hard, color_X_reducido, col = galaxy_class_X, pch = 16, ylim = c(-1, 2),
     xlab = "Hard X-Rays luminosity excess", ylab = "X-Rays color excess")
abline(h = threshold_colorX_excess, lty = 2, lwd = 2)
abline(v = threshold_hard_excess, lty = 4, lwd = 2)


# par(mar=c(5, 4, 4, 4) + 0.1)
# range_color_X = seq(-3,3,1) # elegimos un rango de valores para el X-Rays color
# range_HR = round(((1/7)*10^range_color_X - 1)/((1/7)*10^range_color_X + 1), 2) # lo pasamos a HR despejando de la ecuación
# abline(h = log10(7), lty = 2, lwd = 2) # HR = 0 ---> color_X = log10(7)
# abline(a = 42, b = -1, lty = 4, lwd = 2) # pintamos también el criterio de selección en luminosidad Hard = 10^42 erg/s
# axis(2, col = "black", las = 0)
# mtext("X-Rays color", side = 2, line = 2.5)
# par(new = TRUE)
# plot(L_Hard, color_X, pch = 16, col = alpha('blue', 0.1),
#      main = '', 
#      xlab = '', ylab = '', type = 'n', axes = FALSE)
# axis(4, labels = range_HR, col = "black", col.axis = "black", las = 0, at = -3:3)
# mtext("HR", side = 4, col = "black", line = 2) 
# axis(1, pretty(range(L_Hard)))
# mtext("log[L (Hard X-Rays)]", side=1, col="black", line = 2.5)


######################################  Sparse Matrix  #########################################


Sparse = array(0, dim = c(5, 5, 5))

X = rep(0, len)


for (i in 1:len){
  if ((is.na(color_X[i]) == TRUE) | (is.na(color_X_reducido[i]) == TRUE) | (is.na(L_reducida_Hard[i]) == TRUE)){
    X[i] = 5 #undetected (grey)
  }
  else {
    if (L_reducida_Hard[i] > threshold_hard_excess){
      if ((color_X_reducido[i] > threshold_colorX_excess)){
        X[i] = 1 #active (red)
      }
      else {
        X[i] = 2 #luminous (green)
      }
    }
    else {
      if ((color_X_reducido[i] > threshold_colorX_excess)){
        X[i] = 3 #color (yellow)
      }
      else {
        X[i] = 4 #normal (blue)
      }
    }
  }
}



IR = rep(0, len)


for (i in 1:len){
  if ((is.na(color_W2_W1[i]) == TRUE) | (is.na(color_W2_W1_reducido[i]) == TRUE) | (is.na(L_reducida_IR2[i]) == TRUE)){
    IR[i] = 5 #undetected (grey)
  }
  else {
    if (L_reducida_IR2[i] > threshold_ir2_excess){
      if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
        IR[i] = 1 #active (red)
      }
      else {
        IR[i] = 2 #luminous (green)
      }
    }
    else {
      if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
        IR[i] = 3 #color (yellow)
      }
      else {
        IR[i] = 4 #normal (blue)
      }
    }
  }
}



R = rep(0, len)


for (i in 1:len){
  if ((is.na(color_Radio_IR[i]) == TRUE) | (is.na(color_Radio_IR_reducido[i]) == TRUE) | (is.na(L_reducida_Radio[i]) == TRUE)){
    R[i] = 5 #undetected (grey)
  }
  else {
    if (L_reducida_Radio[i] > threshold_radio_excess){
      if ((color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
        R[i] = 1 #active (red)
      }
      else {
        R[i] = 2 #luminous (green)
      }
    }
    else {
      if ((color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
        R[i] = 3 #color (yellow)
      }
      else {
        R[i] = 4 #normal (blue)
      }
    }
  }
}


for (i in 1:len){
  Sparse[X[i], IR[i], R[i]] = Sparse[X[i], IR[i], R[i]] + 1
}

