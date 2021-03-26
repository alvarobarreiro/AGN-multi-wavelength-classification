
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
HR <- (2*sdss_xmatch$Hard_flux - 5*sdss_xmatch$Soft_flux)/(2*sdss_xmatch$Hard_flux + 5*sdss_xmatch$Soft_flux)


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

SFR_X <- seq(min(na.omit(SFR)), max(na.omit(SFR)), (max(na.omit(SFR)) - min(na.omit(SFR)))/249)
SFR_grid <- matrix(SFR_X, nrow = 250, ncol = 250, byrow = FALSE)

M_y <- seq(min(na.omit(M)), max(na.omit(M)), (max(na.omit(M)) - min(na.omit(M)))/249)
M_grid <- matrix(M_y, nrow = 250, ncol = 250, byrow = TRUE)
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



############# HISTOGRAMAS 2D PARA LA DENSIDAD Y LA DENSIDAD PESADA #############

library(colorRamps)
library(viridis)
library(fields)
library(scales)

image.plot(SFR_X, M_y, densidad_IR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density IR (W4 band)')

image.plot(SFR_X, M_y, densidad_Radio, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density Radio')

image.plot(SFR_X, M_y, densidad_Soft, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density Soft X-Rays')

image.plot(SFR_X, M_y, densidad_Hard, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density Hard X-Rays')

image.plot(SFR_X, M_y, densidad_W1, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density IR (W1 band)')

image.plot(SFR_X, M_y, densidad_W2, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density IR (W2 band)')

image.plot(SFR_X, M_y, densidad_colorIR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density W1-W2')

image.plot(SFR_X, M_y, densidad_color_W2_W1, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density IR color')

image.plot(SFR_X, M_y, densidad_q24, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density q24')

image.plot(SFR_X, M_y, densidad_color_Radio_IR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density Radio - IR (W4) color')

image.plot(SFR_X, M_y, densidad_colorX, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Density X-Rays color')




image.plot(SFR_X, M_y, densidad_pesada_IR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5), 
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: IR (W4 band) weighted density')

image.plot(SFR_X, M_y, densidad_pesada_Radio, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5), 
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Radio weighted density')

image.plot(SFR_X, M_y, densidad_pesada_Soft, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5), 
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Soft X-Rays weighted density')

image.plot(SFR_X, M_y, densidad_pesada_Hard, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5), 
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Hard X-Rays weighted density')

image.plot(SFR_X, M_y, densidad_pesada_W1, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: IR (W1 band) weighted density')

image.plot(SFR_X, M_y, densidad_pesada_W2, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: IR (W2 band) weighted density')

image.plot(SFR_X, M_y, densidad_pesada_colorIR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: W1-W2 weighted density')

image.plot(SFR_X, M_y, densidad_pesada_color_W2_W1, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: IR color weighted density')

image.plot(SFR_X, M_y, densidad_pesada_q24, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: q24 weighted density')

image.plot(SFR_X, M_y, densidad_pesada_color_Radio_IR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Radio - IR(W4) color weighted density')

image.plot(SFR_X, M_y, densidad_pesada_colorX, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: X-Rays color weighted density')


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


############ HISTOGRAMAS 2D PARA LOS <L> Y <COLOR> #############

image.plot(SFR_X, M_y, L_mean_IR, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR (W4 band)')

image.plot(SFR_X, M_y, L_mean_Radio, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Radio')

image.plot(SFR_X, M_y, L_mean_Soft, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Soft X-Rays')

image.plot(SFR_X, M_y, L_mean_Hard, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Hard X-Rays')

image.plot(SFR_X, M_y, L_mean_IR1, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR (W1 band)')

image.plot(SFR_X, M_y, L_mean_IR2, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR (W2 band)')

image.plot(SFR_X, M_y, colorIR_mean, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <W1-W2>')

image.plot(SFR_X, M_y, color_W2_W1_mean, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <IR color>')

image.plot(SFR_X, M_y, q24_mean, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <q24>')

image.plot(SFR_X, M_y, color_Radio_IR_mean, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <Radio-IR(W4) color>')

image.plot(SFR_X, M_y, colorX_mean, col = viridis(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <X-Rays color>')


############## HISTOGRAMAS 2D PARA LAS DISPERSIONES DE <L> Y <COLOR> ##############

image.plot(SFR_X, M_y, Dispersion_IR, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W4 band)')

image.plot(SFR_X, M_y, Dispersion_Radio, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio')

image.plot(SFR_X, M_y, Dispersion_Soft, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Soft X-Rays')

image.plot(SFR_X, M_y, Dispersion_Hard, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
      ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Hard X-Rays')

image.plot(SFR_X, M_y, Dispersion_IR1, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W1 band)')

image.plot(SFR_X, M_y, Dispersion_IR2, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W2 band)')

image.plot(SFR_X, M_y, Dispersion_colorIR, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion W1-W2')

image.plot(SFR_X, M_y, Dispersion_color_W2_W1, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR color')

image.plot(SFR_X, M_y, Dispersion_q24, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion q24')

image.plot(SFR_X, M_y, Dispersion_color_Radio_IR, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio-IR(W4) color')

image.plot(SFR_X, M_y, Dispersion_colorX, col = matlab.like(1000), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion X-Rays color')


############# HOMOGENEIZAR LOS RANGOS DINÁMICOS DE LOS HISTOGRAMAS 2D ################

# Acotamos el rango dinámico de todas las luminosidades al intervalo (37.5 - 43.5):

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


# Los colores van más a su bola, así que aplicamos un criterio de percentiles para cada uno de ellos:

colorIR_truncated <- colorIR_mean
colorIR_truncated[colorIR_truncated > quantile(colorIR_mean, probs = 0.85)] <- quantile(colorIR_mean, probs = 0.85)
colorIR_truncated[colorIR_truncated < quantile(colorIR_mean, probs = 0.15)] <- quantile(colorIR_mean, probs = 0.15)

q24_truncated <- q24_mean
q24_truncated[q24_truncated > quantile(q24_mean, probs = 0.85)] <- quantile(q24_mean, probs = 0.85)
q24_truncated[q24_truncated < quantile(q24_mean, probs = 0.15)] <- quantile(q24_mean, probs = 0.15)

color_W2_W1_truncated <- color_W2_W1_mean
color_W2_W1_truncated[color_W2_W1_truncated > quantile(color_W2_W1_mean, probs = 0.85)] <- quantile(color_W2_W1_mean, probs = 0.85)
color_W2_W1_truncated[color_W2_W1_truncated < quantile(color_W2_W1_mean, probs = 0.15)] <- quantile(color_W2_W1_mean, probs = 0.15)

color_Radio_IR_truncated <- color_Radio_IR_mean
color_Radio_IR_truncated[color_Radio_IR_truncated > quantile(color_Radio_IR_mean, probs = 0.85)] <- quantile(color_Radio_IR_mean, probs = 0.85)
color_Radio_IR_truncated[color_Radio_IR_truncated < quantile(color_Radio_IR_mean, probs = 0.15)] <- quantile(color_Radio_IR_mean, probs = 0.15)


# colorX_truncated <- colorX_mean
# colorX_truncated[colorX_truncated > quantile(colorX_mean, probs = 0.85)] <- quantile(colorX_mean, probs = 0.85)
# colorX_truncated[colorX_truncated < quantile(colorX_mean, probs = 0.15)] <- quantile(colorX_mean, probs = 0.15)

colorX_truncated <- colorX_mean
colorX_truncated[colorX_truncated > quantile(na.omit(color_X), probs = 0.85)] <- quantile(na.omit(color_X), probs = 0.85)
colorX_truncated[colorX_truncated < quantile(na.omit(color_X), probs = 0.15)] <- quantile(na.omit(color_X), probs = 0.15)

# Para las dispersiones acotamos el rango dinámico al intervalo (0.1 - 1):

Dispersion_IR_truncated <- Dispersion_IR
Dispersion_IR_truncated[Dispersion_IR_truncated > 
                          (quantile(Dispersion_IR, probs = 0.85) - quantile(Dispersion_IR, probs = 0.15))/2] = 
  (quantile(Dispersion_IR, probs = 0.85) - quantile(Dispersion_IR, probs = 0.15))/2
Dispersion_IR_truncated[Dispersion_IR_truncated < 0] = 0

Dispersion_Radio_truncated <- Dispersion_Radio
Dispersion_Radio_truncated[Dispersion_Radio_truncated > 
                             (quantile(Dispersion_Radio, probs = 0.85) - quantile(Dispersion_Radio, probs = 0.15))/2] =
  (quantile(Dispersion_Radio, probs = 0.85) - quantile(Dispersion_Radio, probs = 0.15))/2
Dispersion_Radio_truncated[Dispersion_Radio_truncated < 0] = 0

Dispersion_Soft_truncated <- Dispersion_Soft
Dispersion_Soft_truncated[Dispersion_Soft_truncated > 
                            (quantile(Dispersion_Soft, probs = 0.85) - quantile(Dispersion_Soft, probs = 0.15))/2] =
  (quantile(Dispersion_Soft, probs = 0.85) - quantile(Dispersion_Soft, probs = 0.15))/2
Dispersion_Soft_truncated[Dispersion_Soft_truncated < 0] = 0

Dispersion_Hard_truncated <- Dispersion_Hard
Dispersion_Hard_truncated[Dispersion_Hard_truncated > 
                            (quantile(Dispersion_Hard, probs = 0.85) - quantile(Dispersion_Hard, probs = 0.15))/2] =
  (quantile(Dispersion_Hard, probs = 0.85) - quantile(Dispersion_Hard, probs = 0.15))/2
Dispersion_Hard_truncated[Dispersion_Hard_truncated < 0] = 0

Dispersion_IR1_truncated <- Dispersion_IR1
Dispersion_IR1_truncated[Dispersion_IR1_truncated > 
                           (quantile(Dispersion_IR1, probs = 0.85) - quantile(Dispersion_IR1, probs = 0.15))/2] =
  (quantile(Dispersion_IR1, probs = 0.85) - quantile(Dispersion_IR1, probs = 0.15))/2
Dispersion_IR1_truncated[Dispersion_IR1_truncated < 0] <- 0

Dispersion_IR2_truncated <- Dispersion_IR2
Dispersion_IR2_truncated[Dispersion_IR2_truncated > 
                           (quantile(Dispersion_IR2, probs = 0.85) - quantile(Dispersion_IR2, probs = 0.15))/2] =
  (quantile(Dispersion_IR2, probs = 0.85) - quantile(Dispersion_IR2, probs = 0.15))/2
Dispersion_IR2_truncated[Dispersion_IR2_truncated < 0] <- 0

# Dispersion_colorIR_truncated <- Dispersion_colorIR
# Dispersion_colorIR_truncated[Dispersion_colorIR_truncated > 
#                                (quantile(Dispersion_color_W2_W1, probs = 0.85) - quantile(Dispersion_color_W2_W1, probs = 0.15))/2] =
#   quantile(Dispersion_color_W2_W1, probs = 0.85)
# Dispersion_colorIR_truncated[Dispersion_colorIR_truncated < 0] = 0
# 
# Dispersion_q24_truncated <- Dispersion_q24
# Dispersion_q24_truncated[Dispersion_q24_truncated > 
#                            (quantile(Dispersion_q24, probs = 0.85) - quantile(Dispersion_q24, probs = 0.15))/2] =
#   (quantile(Dispersion_q24, probs = 0.85) - quantile(Dispersion_q24, probs = 0.15))/2
# Dispersion_q24_truncated[Dispersion_q24_truncated < 0] <- 0

Dispersion_color_W2_W1_truncated <- Dispersion_color_W2_W1
Dispersion_color_W2_W1_truncated[Dispersion_color_W2_W1_truncated >
                                   (quantile(Dispersion_, probs = 0.85) - quantile(Dispersion_Radio, probs = 0.15))/2] =
  quantile(Dispersion_color_W2_W1, probs = 0.85)
Dispersion_color_W2_W1_truncated[Dispersion_color_W2_W1_truncated < 0] = 0

Dispersion_color_Radio_IR_truncated <- Dispersion_color_Radio_IR
Dispersion_color_Radio_IR_truncated[Dispersion_color_Radio_IR_truncated > 
                                      (quantile(Dispersion_color_Radio_IR, probs = 0.85) - quantile(Dispersion_color_Radio_IR, probs = 0.15))/2] = 
  (quantile(Dispersion_color_Radio_IR, probs = 0.85) - quantile(Dispersion_color_Radio_IR, probs = 0.15))/2
Dispersion_color_Radio_IR_truncated[Dispersion_color_Radio_IR_truncated < 0] <- 0

Dispersion_colorX_truncated <- Dispersion_colorX
Dispersion_colorX_truncated[Dispersion_colorX_truncated > 
                              (quantile(Dispersion_colorX, probs = 0.85) - quantile(Dispersion_colorX, probs = 0.15))/2] =
  (quantile(Dispersion_colorX, probs = 0.85) - quantile(Dispersion_colorX, probs = 0.15))/2
Dispersion_colorX_truncated[Dispersion_colorX_truncated < 0] <- 0



############ HISTOGRAMAS 2D DE LOS NUEVOS <L> Y <COLOR> Y DE SUS RESPECTIVAS DISPERSIONES ################



image.plot(SFR_X, M_y, L_mean_IR, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W4 band')
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, L_mean_IR1_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W1 band')
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, L_mean_IR2_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> IR W2 band')
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, L_mean_Radio_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Radio')
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, L_mean_Soft_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Soft X-Rays')
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, L_mean_Hard_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Hard X-Rays')
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .95), add = TRUE, lty = 2)

# image.plot(SFR_X, M_y, colorIR_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
#            ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <W1-W2>')
# 
# image.plot(SFR_X, M_y, q24_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
#            ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <q24>')

image.plot(SFR_X, M_y, color_W2_W1_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <IR color>')
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, color_Radio_IR_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <Radio-IR(W4) color>')
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, colorX_truncated, col = matlab.like(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <X-Rays color>')
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .95), add = TRUE, lty = 2)



image.plot(SFR_X, M_y, Dispersion_IR_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W4 band)')
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_IR, levels = quantile(densidad_IR, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_Radio_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio')
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Radio, levels = quantile(densidad_Radio, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_Soft_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Soft X-Rays')
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Soft, levels = quantile(densidad_Soft, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_Hard_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Hard X-Rays')
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_Hard, levels = quantile(densidad_Hard, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_IR1_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W1 band)')
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W1, levels = quantile(densidad_W1, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_IR2_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR (W2 band)')
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_W2, levels = quantile(densidad_W2, probs = .95), add = TRUE, lty = 2)

# image.plot(SFR_X, M_y, Dispersion_colorIR_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
#            ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion W1-W2')
# 
# image.plot(SFR_X, M_y, Dispersion_q24_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
#            ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion q24')

image.plot(SFR_X, M_y, Dispersion_color_W2_W1_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion IR color')
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_W2_W1, levels = quantile(densidad_color_W2_W1, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_color_Radio_IR_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion Radio-IR(W4) color')
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_color_Radio_IR, levels = quantile(densidad_color_Radio_IR, probs = .95), add = TRUE, lty = 2)

image.plot(SFR_X, M_y, Dispersion_colorX_truncated, col = viridis(15), xlim = c(-3, 1.5), ylim = c(8.5, 11.5),
           ylab = 'log[M* (solar masses)]', xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: Dispersion X-Rays color')
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .50), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .75), add = TRUE, lty = 2)
contour(SFR_X, M_y, densidad_colorX, levels = quantile(densidad_colorX, probs = .95), add = TRUE, lty = 2)

############ CÁLCULO DE LOS LUMINOSIDADES Y COLORES MEDIOS Y DE SUS EXCESOS ############

L_media_soft <- rep(0, len)
L_media_hard <- rep(0, len)
L_media_radio <- rep(0, len)

L_media_ir <- rep(0, len)
L_media_ir1 <- rep(0, len)
L_media_ir2 <- rep(0, len)
colorIR_media <- rep(0, len)
color_W2_W1_media <- rep(0, len)
q24_media <- rep(0, len)
color_Radio_IR_media <- rep(0, len)
colorX_media <- rep(0, len)

binsize_SFR <- (SFR_grid[2,1] - SFR_grid[1,1])/2
binsize_M <- (M_grid[1,2] - M_grid[1,1])/2


L_mean_vector <- function(L, L_media, L_mean){
  for (i in 1:len){
    if ((is.na(L[i]) == FALSE) & (is.na(SFR[i]) == FALSE) & (is.na(M[i]) == FALSE)){
      fila = 0
      columna = 0
      for (j in 1:length(SFR_X)){
        if (abs(SFR[i] - SFR_grid[j,1]) < binsize_SFR){
          fila = j
          break
        }
      }
      for (k in 1:length(M_y)){
        if (abs(M[i] - M_grid[1,k]) < binsize_M){
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
colorIR_media <- L_mean_vector(L = color_IR, L_media = colorIR_media, L_mean = colorIR_mean)
color_W2_W1_media <- L_mean_vector(L = color_W2_W1, L_media = color_W2_W1_media, L_mean = color_W2_W1_mean)
q24_media <- L_mean_vector(L = q_24, L_media = q24_media, L_mean = q24_mean)
color_Radio_IR_media <- L_mean_vector(L = color_Radio_IR, L_media = color_Radio_IR_media, L_mean = color_Radio_IR_mean)
colorX_media <- L_mean_vector(L = color_X, L_media = colorX_media, L_mean = colorX_mean)

L_reducida_Soft <- L_Soft - L_media_soft
L_reducida_Hard <- L_Hard - L_media_hard
L_reducida_Radio <- L_Radio - L_media_radio
L_reducida_IR <- L_IR - L_media_ir
L_reducida_IR1 <- L_IR1 - L_media_ir1
L_reducida_IR2 <- L_IR2 - L_media_ir2
color_IR_reducido <- color_IR - colorIR_media
color_W2_W1_reducido <- color_W2_W1 - color_W2_W1_media
q24_reducido <- q_24 - q24_media
color_Radio_IR_reducido <- color_Radio_IR - color_Radio_IR_media
color_X_reducido <- color_X - colorX_media

#################################   CONTOUR PLOTS 18-05-20   ###############################

# image.plot(SFR_X, M_y, L_mean_Hard_truncated, col = matlab.like(1000), ylab = 'log[M* (solar masses)]',
#            xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Hard X-Rays', xlim = c(-3, 2), ylim = c(8.5, 11.5))
# contour(SFR_X, M_y, L_mean_Hard, levels = 42, lwd = 2, add = TRUE, lty = 1)
# contour(SFR_X, M_y, L_mean_Hard + Dispersion_Hard, levels = 42, lwd = 2, add = TRUE, lty = 2)
# contour(SFR_X, M_y, L_mean_Hard + 2*Dispersion_Hard, levels = 42, lwd = 2, add = TRUE, lty = 3)
# 
# image.plot(SFR_X, M_y, colorIR_mean, col = matlab.like(1000), ylab = 'log[M* (solar masses)]', 
#            xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <W1-W2>', xlim = c(-3, 2), ylim = c(8.5, 11.5))
# contour(SFR_X, M_y, colorIR_mean, levels = 0.8, lwd = 2, add = TRUE, lty = 1)
# contour(SFR_X, M_y, colorIR_mean + Dispersion_colorIR, levels = 0.8, lwd = 2, add = TRUE, lty = 2)
# contour(SFR_X, M_y, colorIR_mean + 2*Dispersion_colorIR, levels = 0.8, lwd = 2, add = TRUE, lty = 3)
# 
# 
# image.plot(SFR_X, M_y, q24_truncated, col = matlab.like(1000), ylab = 'log[M* (solar masses)]',
#            xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <q24>', xlim = c(-3, 2), ylim = c(8.5, 11.5))
# contour(SFR_X, M_y, q24_mean, levels = -0.23, lwd = 2, add = TRUE, lty = 1)
# contour(SFR_X, M_y, q24_mean - Dispersion_q24, levels = -0.23, lwd = 2, add = TRUE, lty = 2)
# contour(SFR_X, M_y, q24_mean - 2*Dispersion_q24, levels = -0.23, lwd = 2, add = TRUE, lty = 3)
# 
# 
# image.plot(SFR_X, M_y, colorX_truncated, col = matlab.like(1000), ylab = 'log[M* (solar masses)]',
#            xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <X-Rays color>', xlim = c(-3, 2), ylim = c(8.5, 11.5))
# contour(SFR_X, M_y, colorX_mean, levels = log10(5/2), lwd = 2, add = TRUE, lty = 1)
# contour(SFR_X, M_y, colorX_mean + Dispersion_colorX, levels = log10(5/2), lwd = 2, add = TRUE, lty = 2)
# contour(SFR_X, M_y, colorX_mean + 2*Dispersion_colorX, levels = log10(5/2), lwd = 2, add = TRUE, lty = 3)
# 
# 
# image.plot(SFR_X, M_y, L_mean_Radio_truncated, col = matlab.like(1000), ylab = 'log[M* (solar masses)]',
#            xlab = 'log[SFR (solar masses/yr)]', main = 'Histogram 2D: <L> Radio', xlim = c(-3, 2), ylim = c(8.5, 11.5))
# contour(SFR_X, M_y, L_mean_Radio, levels = 40, lwd = 2, add = TRUE, lty = 1)
# contour(SFR_X, M_y, L_mean_Radio + Dispersion_Radio, levels = 40, lwd = 2, add = TRUE, lty = 2)
# contour(SFR_X, M_y, L_mean_Radio + 2*Dispersion_Radio, levels = 40, lwd = 2, add = TRUE, lty = 3)


############################### HISTOGRAMAS 01-06-20 ######################################

h_soft <- hist(L_Soft, col = 'skyblue', breaks = (max(na.omit(L_Soft)) - min(na.omit(L_Soft)))/(0.3*sd(na.omit(L_Soft))),
               freq = FALSE, main = 'Histogram of log(L) (Soft X-Rays)')

h_hard <- hist(L_Hard, col = 'darkblue', breaks = (max(na.omit(L_Hard)) - min(na.omit(L_Hard)))/(0.3*sd(na.omit(L_Hard))),
               freq = FALSE, main = 'Histogram of log(L) (Hard X-Rays)')

h_radio <- hist(L_Radio, col = 'green', breaks = (max(na.omit(L_Radio)) - min(na.omit(L_Radio)))/(0.3*sd(na.omit(L_Radio))),
                freq = FALSE, main = 'Histogram of log(L) (Radio)')

h_ir <- hist(L_IR, col = 'red', breaks = (max(na.omit(L_IR)) - min(na.omit(L_IR)))/(0.3*sd(na.omit(L_IR))),
             freq = FALSE, main = 'Histogram of log(L) (IR W4 band)')

h_ir1 <- hist(L_IR1, col = 'yellow', breaks = (max(na.omit(L_IR1)) - min(na.omit(L_IR1)))/(0.3*sd(na.omit(L_IR1))),
              freq = FALSE, main = 'Histogram of log(L) (IR W1 band)')

h_ir2 <- hist(L_IR2, col = 'orange', breaks = (max(na.omit(L_IR2)) - min(na.omit(L_IR2)))/(0.3*sd(na.omit(L_IR2))),
              freq = FALSE, main = 'Histogram of log(L) (IR W2 band)')

h_w1_w2 <- hist(color_IR, col = 'darkred', breaks = (max(na.omit(color_IR)) - min(na.omit(color_IR)))/(0.3*sd(na.omit(color_IR))),
                freq = FALSE, main = 'Histogram of (W1-W2)')
# abline(v = moda_color_IR + 2*sigma_izquierda_colorIR, lty = 1, lwd = 2)
# abline(v = moda_color_IR - 2*sigma_izquierda_colorIR, lty = 1, lwd = 2)
# abline(v = 0.8, lty = 2, lwd = 2)

h_color_W2_W1 <- hist(color_W2_W1, col = 'darkred', breaks = (max(na.omit(color_W2_W1)) - min(na.omit(color_W2_W1)))/(0.3*sd(na.omit(color_W2_W1))),
                      freq = FALSE, main = 'Histogram of IR color')
# abline(v = 0.4*0.8 + W1_W2_norm, lty = 2, lwd = 2)
# abline(v = moda_color_W2_W1 - 2*sigma_izquierda_color_W2_W1, lty = 1, lwd = 2)
# abline(v = moda_color_W2_W1 + 2*sigma_izquierda_color_W2_W1, lty = 1, lwd = 2)

h_q24 <- hist(q_24, col = 'darkgreen', breaks = (max(na.omit(q_24)) - min(na.omit(q_24)))/(0.3*sd(na.omit(q_24))), 
              freq = FALSE, main = 'Histogram of q24')
# abline(v = moda_q24 - 2*sigma_izquierda_q24, lty = 1, lwd = 2)
# abline(v = moda_q24 + 2*sigma_izquierda_q24, lty = 1, lwd = 2)
# abline(v = -0.23, lty = 2, lwd = 2)

h_color_Radio_IR <- hist(color_Radio_IR, col = 'darkgreen', breaks = (max(na.omit(color_Radio_IR)) - min(na.omit(color_Radio_IR)))/(0.3*sd(na.omit(color_Radio_IR))), 
                         freq = FALSE, main = 'Histogram of Radio-IR(W4) color')

h_colorX <- hist(color_X, col = 'purple', breaks = (max(na.omit(color_X)) - min(na.omit(color_X)))/(0.3*sd(na.omit(color_X))),
                 freq = FALSE, main = 'Histogram of colorX', xlim = c(-2, 3.5))
# abline(v = moda_colorX + 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
# abline(v = moda_colorX - 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
# abline(v = log10(5/2), lty = 2, lwd = 2)

h_w1_w2_red <- hist(color_IR_reducido, col = 'darkred',
                    breaks = (max(na.omit(color_IR_reducido)) - min(na.omit(color_IR_reducido)))/(0.3*sd(na.omit(color_IR_reducido))),
                    freq = FALSE, main = expression(paste('Histogram of ', Delta, "(W1-W2)")), xlab = expression(paste(Delta, "(W1-W2)")))

h_color_W2_W1_red <- hist(color_W2_W1_reducido, col = 'darkred', breaks = (max(na.omit(color_W2_W1_reducido)) - min(na.omit(color_W2_W1_reducido)))/(0.3*sd(na.omit(color_W2_W1_reducido))),
                          freq = FALSE, main = expression(paste('Histogram of ', Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])")), 
                          xlab = expression(paste(Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])")))

h_q24_red <- hist(q24_reducido, col = 'darkgreen', freq = FALSE, breaks = (max(na.omit(q24_reducido)) - min(na.omit(q24_reducido)))/(0.3*sd(na.omit(q24_reducido))),
                  main = expression(paste('Histogram of ', Delta, "q24")),
                  xlab = expression(paste(Delta, "q24")))

h_color_Radio_IR_red <- hist(color_Radio_IR_reducido, col = 'darkgreen', freq = FALSE, breaks = (max(na.omit(color_Radio_IR_reducido)) - min(na.omit(color_Radio_IR_reducido)))/(0.3*sd(na.omit(color_Radio_IR_reducido))),
                             main = expression(paste('Histogram of ', Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(IR W4)])")),
                             xlab = expression(paste(Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(IR W4)])")))

h_hard_red <- hist(L_reducida_Hard, col = 'darkblue', freq = FALSE, breaks = (max(na.omit(L_reducida_Hard)) - min(na.omit(L_reducida_Hard)))/(0.3*sd(na.omit(L_reducida_Hard))),
                   main = expression(paste('Histogram of ', Delta, "(log[L(Hard X-Rays)])")),
                   xlab = expression(paste(Delta, "(log[L(Hard X-Rays)])")))

h_colorX_red <- hist(color_X_reducido, col = 'darkblue', freq = FALSE, breaks = (max(na.omit(color_X_reducido)) - min(na.omit(color_X_reducido)))/(0.3*sd(na.omit(color_X_reducido))),
                     main = expression(paste('Histogram of ', Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])")),
                     xlab = expression(paste(Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])")))

h_ir2_excess <- hist(L_reducida_IR2, col = 'darkred', breaks = (max(na.omit(L_reducida_IR2)) - min(na.omit(L_reducida_IR2)))/(0.3*sd(na.omit(L_reducida_IR2))),
                     freq = FALSE, main = expression(paste('Histogram of ', Delta, "(log[L(IR W2 band)])")))

h_radio_excess <- hist(L_reducida_Radio, col = 'darkgreen', breaks = (max(na.omit(L_reducida_Radio)) - min(na.omit(L_reducida_Radio)))/(0.3*sd(na.omit(L_reducida_Radio))),
                     freq = FALSE, main = expression(paste('Histogram of ', Delta, "(log[L(Radio)])")))

h_hard_excess <- hist(L_reducida_Hard, col = 'darkblue', breaks = (max(na.omit(L_reducida_Hard)) - min(na.omit(L_reducida_Hard)))/(0.3*sd(na.omit(L_reducida_Hard))),
                       freq = FALSE, main = expression(paste('Histogram of ', Delta, "(log[L(Hard X-Rays)])")))



################ CÁLCULO DE LAS MODAS Y DE LAS SIGMAS UTILIZANDO LOS HISTOGRAMAS ##################

moda_ir2_excess <- h_ir2_excess$mids[which.max(h_ir2_excess$counts)]
S_ir2_excess_izq <- as.numeric()
ir2_excess_h <- na.omit(L_reducida_IR2)
for (i in 1:length(ir2_excess_h)){
  if (ir2_excess_h[i] <= moda_ir2_excess){
    S_ir2_excess_izq <- c(S_ir2_excess_izq, ir2_excess_h[i])
  }
}
sigma_izquierda_ir2_excess <- sqrt(sum((S_ir2_excess_izq - moda_ir2_excess)^2)/(length(S_ir2_excess_izq) - 1))

moda_radio_excess <- h_radio_excess$mids[which.max(h_radio_excess$counts)]
S_radio_excess_izq <- as.numeric()
radio_excess_h <- na.omit(L_reducida_Radio)
for (i in 1:length(radio_excess_h)){
  if (radio_excess_h[i] <= moda_radio_excess){
    S_radio_excess_izq <- c(S_radio_excess_izq, radio_excess_h[i])
  }
}
sigma_izquierda_radio_excess <- sqrt(sum((S_radio_excess_izq - moda_radio_excess)^2)/(length(S_radio_excess_izq) - 1))


moda_hard_excess <- h_hard_excess$mids[which.max(h_hard_excess$counts)]
S_hard_excess_izq <- as.numeric()
hard_excess_h <- na.omit(L_reducida_Hard)
for (i in 1:length(hard_excess_h)){
  if (hard_excess_h[i] <= moda_hard_excess){
    S_hard_excess_izq <- c(S_hard_excess_izq, hard_excess_h[i])
  }
}
sigma_izquierda_hard_excess <- sqrt(sum((S_hard_excess_izq - moda_hard_excess)^2)/(length(S_hard_excess_izq) - 1))


moda_color_W2_W1 <- h_color_W2_W1$mids[which.max(h_color_W2_W1$counts)]
S_color_W2_W1_izq <- as.numeric() # vector con los valores menores que la moda
color_W2_W1_h <- na.omit(color_W2_W1) # nos quedamos sólo con los datos para los cuales tenemos un valor numérico
for (i in 1:length(color_W2_W1_h)){
  if (color_W2_W1_h[i] <= moda_color_W2_W1){
    S_color_W2_W1_izq <- c(S_color_W2_W1_izq, color_W2_W1_h[i]) # metemos los menores que la moda en una lista
  }
}
sigma_izquierda_color_W2_W1 <- sqrt(sum((S_color_W2_W1_izq - moda_color_W2_W1)^2)/(length(S_color_W2_W1_izq) - 1))
# y calculamos la desviación típica para los elementos de dicha lista


moda_color_Radio_IR <- h_color_Radio_IR$mids[which.max(h_color_Radio_IR$counts)]
S_color_Radio_IR_izq <- as.numeric()
color_Radio_IR_h <- na.omit(color_Radio_IR)
for (i in 1:length(color_Radio_IR_h)){
  if (color_Radio_IR_h[i] <= moda_color_Radio_IR){ # aquí lo hacemos al revés porque estamos interesados en los valores
    # mayores que la moda que son los que en principio se identifican con galaxias no activas
    S_color_Radio_IR_izq <- c(S_color_Radio_IR_izq, color_Radio_IR_h[i])
  }
}
sigma_izquierda_color_Radio_IR <- sqrt(sum((S_color_Radio_IR_izq - moda_color_Radio_IR)^2)/(length(S_color_Radio_IR_izq) - 1))

moda_colorX <- h_colorX$mids[which.max(h_colorX$counts)]
S_colorX_izq <- as.numeric()
color_X_h <- na.omit(color_X)
for (i in 1:length(color_X_h)){
  if (color_X_h[i] <= moda_colorX){
    S_colorX_izq <- c(S_colorX_izq, color_X_h[i])
  }
}
sigma_izquierda_colorX <- sqrt(sum((S_colorX_izq - moda_colorX)^2)/(length(S_colorX_izq) - 1))

moda_hard <- h_hard$mids[which.max(h_hard$counts)]
S_hard_izq <- as.numeric()
L_Hard_h <- na.omit(L_Hard)
for (i in 1:length(L_Hard_h)){
  if (L_Hard_h[i] <= moda_hard){
    S_hard_izq <- c(S_hard_izq, L_Hard_h[i])
  }
}
sigma_izquierda_hard <- sqrt(sum((S_hard_izq - moda_hard)^2)/(length(S_hard_izq) - 1))
# moda_hard + sigma_izquierda_hard = 42.5
# moda_hard + 2*sigma_izquierda_hard = 43.5


moda_color_W2_W1_red <- h_color_W2_W1_red$mids[which.max(h_color_W2_W1_red$counts)]
S_color_W2_W1_red_izq <- as.numeric()
color_W2_W1_h_red <- na.omit(color_W2_W1_reducido)
for (i in 1:length(color_W2_W1_h_red)){
  if (color_W2_W1_h_red[i] <= moda_color_W2_W1_red){
    S_color_W2_W1_red_izq <- c(S_color_W2_W1_red_izq, color_W2_W1_h_red[i])
  }
}
sigma_izquierda_color_W2_W1_red <- sqrt(sum((S_color_W2_W1_red_izq - moda_color_W2_W1_red)^2)/(length(S_color_W2_W1_red_izq) - 1))


moda_color_Radio_IR_red <- h_color_Radio_IR_red$mids[which.max(h_color_Radio_IR_red$counts)]
S_color_Radio_IR_red_izq <- as.numeric()
color_Radio_IR_h_red <- na.omit(color_Radio_IR_reducido)
for (i in 1:length(color_Radio_IR_h_red)){
  if (color_Radio_IR_h_red[i] <= moda_color_Radio_IR_red){ # aquí lo hacemos al revés porque estamos interesados en los valores
    # mayores que la moda que son los que en principio se identifican con galaxias no activas
    S_color_Radio_IR_red_izq <- c(S_color_Radio_IR_red_izq, color_Radio_IR_h_red[i])
  }
}
sigma_izquierda_color_Radio_IR_red <- sqrt(sum((S_color_Radio_IR_red_izq - moda_color_Radio_IR_red)^2)/(length(S_color_Radio_IR_red_izq) - 1))


moda_colorX_red <- h_colorX_red$mids[which.max(h_colorX_red$counts)]
S_colorX_red_izq <- as.numeric()
color_X_h_red <- na.omit(color_X_reducido)
for (i in 1:length(color_X_h_red)){
  if (color_X_h_red[i] <= moda_colorX_red){
    S_colorX_red_izq <- c(S_colorX_red_izq, color_X_h_red[i])
  }
}
sigma_izquierda_colorX_red <- sqrt(sum((S_colorX_red_izq - moda_colorX_red)^2)/(length(S_colorX_red_izq) - 1))
# moda_colorX + 2*sigma_izquierda_colorX = 1

moda_hard_red <- h_hard_red$mids[which.max(h_hard_red$counts)]
S_hard_red_izq <- as.numeric()
L_Hard_h_red <- na.omit(L_reducida_Hard)
for (i in 1:length(L_Hard_h_red)){
  if (L_Hard_h_red[i] <= moda_hard_red){
    S_hard_red_izq <- c(S_hard_red_izq, L_Hard_h_red[i])
  }
}
sigma_izquierda_hard_red <- sqrt(sum((S_hard_red_izq - moda_hard_red)^2)/(length(S_hard_red_izq) - 1))
# moda_hard + sigma_izquierda_hard = 42.5
# moda_hard + 2*sigma_izquierda_hard = 43.5



############## CÁLCULO DE LOS PARÁMETROS DE CONVERSIÓN ENTRE NUESTROS COLORES Y LOS DE LA LITERATURA ###############

W1_W2_norm <- log10((171*3.4*10^-6)/(309*4.6*10^-6)) # cociente entre flujos cero de WISE
nuradio_nu24 <- log10((1.4*10^9)/((3*10^8)/(22.1*10^-6))) # cociente de frecuencias (1.4 GHz/22.1 micras)


############## CÁLCULO DE LOS THRESHOLDS PARA CADA LUMINOSIDAD Y COLOR ################

threshold_colorIR_ours = moda_color_W2_W1 + 2*sigma_izquierda_color_W2_W1
threshold_colorIR_lit = 0.4*0.5 + W1_W2_norm
threshold_color_Radio_IR_ours = moda_color_Radio_IR + 2*sigma_izquierda_color_Radio_IR
threshold_color_Radio_IR_lit = 0.23 + nuradio_nu24
threshold_colorX_ours = moda_colorX + 2*sigma_izquierda_colorX
threshold_colorX_lit = log10(5/2)

threshold_colorIR_excess = moda_color_W2_W1_red + 2*sigma_izquierda_color_W2_W1_red
threshold_color_Radio_IR_excess = moda_color_Radio_IR_red + 2*sigma_izquierda_color_Radio_IR_red
threshold_colorX_excess = moda_colorX_red + 2*sigma_izquierda_colorX_red

threshold_ir2_excess = moda_ir2_excess + 2*sigma_izquierda_ir2_excess
threshold_radio_excess = moda_radio_excess + 2*sigma_izquierda_radio_excess
threshold_hard_excess = moda_hard_excess + 2*sigma_izquierda_hard_excess

threshold_ir2_excess_new = 0.5635 # estos valores salen del código excesses_threshold
threshold_colorIR_excess_new = 0.0732
threshold_radio_excess_new = 0.9076
threshold_color_Radio_IR_excess_new = 0.8074
threshold_hard_excess_new = 1.6467
threshold_colorX_excess_new = 0.5654

IR_count_red <- 0
IR_count_green <- 0
IR_count_yellow <- 0
IR_count_blue <- 0

Radio_count_red <- 0
Radio_count_green <- 0
Radio_count_yellow <- 0
Radio_count_blue <- 0

X_count_red <- 0
X_count_green <- 0
X_count_yellow <- 0
X_count_blue <- 0


######################################### COLOR PLOTS ##############################################

par(mar=c(5, 4, 4, 4) + 0.1)

plot(L_IR1, color_W2_W1, pch = 16, col = alpha('blue', 0.05),
     main = 'log[L (IR W2)] - log[L (IR W1)] vs log[L (IR W1)]',
     xlab = '', ylab = '', type = 'n', ylim = c(-0.5,0.5), xlim = c(38, 45))
for (i in 1:len){
  if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess)){
    if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
      points(L_IR1[i], color_W2_W1[i], pch = 16, col = alpha('red4', 0.3))
      IR_count_red <- IR_count_red + 1
    }
    else {
      points(L_IR1[i], color_W2_W1[i], pch = 16, col = alpha('green', 0.3))
      IR_count_green <- IR_count_green + 1
    }
  }
  else if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess)){
    if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
      points(L_IR1[i], color_W2_W1[i], pch = 16, col = alpha('gold', 0.07))
      IR_count_yellow <- IR_count_yellow + 1
    }
    else {
      points(L_IR1[i], color_W2_W1[i], pch = 16, col = alpha('blue', 0.02))
      IR_count_blue <- IR_count_blue + 1
    }
  }
}
abline(h = 0.4*0.8 + W1_W2_norm, lty = 2, lwd = 2)
abline(h = 0.4*0.5 + W1_W2_norm, lty = 4, lwd = 2)
# abline(h = (moda_color_W2_W1 + 2*sigma_izquierda_color_W2_W1), lty = 1, lwd = 2)
# abline(h = (moda_color_W2_W1 - 2*sigma_izquierda_color_W2_W1), lty = 1, lwd = 2)
axis(2, col="black", las = 0)  ## las=1 makes horizontal labels
mtext("log[L (IR W2)] - log[L (IR W1)]",side=2,line=2.5)
par(new=TRUE)
plot(L_IR1, color_IR, pch = 16, col = alpha('blue', 0.1),
     main = '', 
     xlab = '', ylab = '', type = 'n', axes = FALSE, ylim = c((-0.5 - W1_W2_norm)/0.4, (0.5 - W1_W2_norm)/0.4))
mtext("W1-W2",side=4,col="black",line=2) 
axis(4, col="black", col.axis="black", las = 0)
axis(1,pretty(range(L_IR1)))
mtext("log[L (IR W1)]",side=1,col="black",line=2.5)




par(mar=c(5, 4, 4, 4) + 0.1)

plot(L_IR, color_Radio_IR, pch = 16, col = alpha('blue', 0.1),
     main = 'log[L (Radio)] - log[L (IR W4)] vs log[L (IR W4)]', 
     xlab = '', ylab = '', type = 'n')
for (i in 1:len){
  if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess)){
    if ((color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
      points(L_IR[i], color_Radio_IR[i], pch = 16, col = alpha('red4', 0.1))
      Radio_count_red <- Radio_count_red + 1
    }
    else {
      points(L_IR[i], color_Radio_IR[i], pch = 16, col = alpha('green', 0.1))
      Radio_count_green <- Radio_count_green + 1
    }
  }
  else if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess)){
    if ((color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
      points(L_IR[i], color_Radio_IR[i], pch = 16, col = alpha('gold', 0.1))
      Radio_count_yellow <- Radio_count_yellow + 1
    }
    else {
      points(L_IR[i], color_Radio_IR[i], pch = 16, col = alpha('blue', 0.1))
      Radio_count_blue <- Radio_count_blue + 1
    }
  }
}
abline(h = 0.23 + nuradio_nu24, lty = 2, lwd = 2)
# abline(h = (moda_color_Radio_IR - 2*sigma_izquierda_color_Radio_IR), lty = 1, lwd = 2)
# abline(h = (moda_color_Radio_IR + 2*sigma_izquierda_color_Radio_IR), lty = 1, lwd = 2)
abline(a = 40, b = -1, lty = 4, lwd = 2)
axis(2, col="black", las = 0)  ## las=1 makes horizontal labels
mtext("log[L (Radio)] - log[L (IR W4)]",side=2,line=2.5)
par(new=TRUE)
plot(L_IR, q_24, pch = 16, col = alpha('blue', 0.1),
     main = '', 
     xlab = '', ylab = '', type = 'n', axes = FALSE, ylim = rev(range(na.omit(q_24))))
mtext("q24",side=4,col="black",line=2) 
axis(4, col="black", col.axis="black", las = 0)
axis(1,pretty(range(L_IR)))
mtext("log[L (IR W4)]",side=1,col="black",line = 2.5)



par(mar=c(5, 4, 4, 4) + 0.1)

plot(L_Soft, color_X, pch = 16, col = alpha('blue', 0.3),
     main = 'log[L (Hard X-Rays)] - log[L (Soft X-Rays)] vs log[L (Soft X-Rays)]',
     xlab = '', ylab = '', type = 'n')
for (i in 1:len){
  if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess)){
    if ((color_X_reducido[i] > threshold_colorX_excess)){
      points(L_Soft[i], color_X[i], pch = 16, col = alpha('red4', 0.3))
      X_count_red <- X_count_red + 1
    }
    else {
      points(L_Soft[i], color_X[i], pch = 16, col = alpha('green', 0.3))
      X_count_green <- X_count_green + 1
    }
  }
  else if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess)){
    if ((color_X_reducido[i] > threshold_colorX_excess)){
      points(L_Soft[i], color_X[i], pch = 16, col = alpha('gold', 0.3))
      X_count_yellow <- X_count_yellow + 1
    }
    else {
      points(L_Soft[i], color_X[i], pch = 16, col = alpha('blue', 0.3))
      X_count_blue <- X_count_blue + 1
    }
  }
}
range_color_X = seq(-3,3,1)
range_HR = round(((2/5)*10^range_color_X - 1)/((2/5)*10^range_color_X + 1), 2)
abline(h = log10(5/2), lty = 2, lwd = 2)
# abline(h = moda_colorX + 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
# abline(h = moda_colorX - 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
abline(a = 42, b = -1, lty = 4, lwd = 2)
axis(2, col="black", las = 0)  ## las=1 makes horizontal labels
mtext("log[L (Hard X-Rays)] - log[L (Soft X-Rays)]",side=2,line=2.5)
par(new=TRUE)
plot(L_Soft, color_X, pch = 16, col = alpha('blue', 0.1),
     main = '', 
     xlab = '', ylab = '', type = 'n', axes = FALSE)
axis(4, labels = range_HR, col="black", col.axis="black", las = 0, at = -3:3)
mtext("HR",side=4,col="black",line=2) 
axis(1,pretty(range(L_Soft)))
mtext("log[L (Soft X-Rays)]",side=1,col="black",line=2.5)


################################### FINAL PLOTS ########################################

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

galaxy_classification = function(L, color, threshold_L, threshold_color, alpha_value){
  N = length(L)
  galaxy_class = as.vector(rep(alpha("black", 1), N))
  for (i in 1:N){
    galaxy_class[i] = get_galaxy_type(L[i], color[i], threshold_L, threshold_color, alpha_value)
  }
  return(galaxy_class)
}

galaxy_class_IR = galaxy_classification(L_reducida_IR, color_W2_W1_reducido,
                                           threshold_ir2_excess_new, threshold_colorIR_excess_new, 0.7)

plot(L_reducida_IR, color_W2_W1_reducido, col = galaxy_class_IR, pch = 16, ylim = c(-0.5, 1),
     xlab = "IR (W2 band) luminosity excess", ylab = "IR color excess")
abline(h = threshold_colorIR_excess_new, lty = 2, lwd = 2)
abline(v = threshold_ir2_excess_new, lty = 4, lwd = 2)

plot(L_IR2, color_W2_W1, col = galaxy_class_IR, pch = 16, ylim = c(-0.75, 0.6),
     xlab = "log[L (IR W2 band)]", ylab = "IR color")
abline(h = threshold_colorIR_lit, lty = 2, lwd = 2)


galaxy_class_Radio = galaxy_classification(L_reducida_Radio, color_Radio_IR_reducido,
                                           threshold_radio_excess_new, threshold_color_Radio_IR_excess_new, 0.3)

plot(L_reducida_Radio, color_Radio_IR_reducido, col = galaxy_class_Radio, pch = 16, xlim = c(-1.5, 3), ylim = c(-2, 3),
     xlab = "Radio luminosity excess", ylab = "Radio-IR color excess")
abline(h = threshold_color_Radio_IR_excess_new, lty = 2, lwd = 2)
abline(v = threshold_radio_excess_new, lty = 4, lwd = 2)

plot(L_Radio, color_Radio_IR, col = galaxy_class_Radio, pch = 16, ylim = c(-6, -0.5),
     xlab = "log[L (Radio)]", ylab = "Radio-IR color")
abline(h = threshold_color_Radio_IR_lit, lty = 2, lwd = 2)
abline(v = 40, lty = 4, lwd = 2)


galaxy_class_X = galaxy_classification(L_reducida_Hard, color_X_reducido,
                                           threshold_hard_excess_new, threshold_colorX_excess_new, 0.7)

plot(L_reducida_Hard, color_X_reducido, col = galaxy_class_X, pch = 16, ylim = c(-1, 2),
     xlab = "Hard X-Rays luminosity excess", ylab = "X-Rays color excess")
abline(h = threshold_colorX_excess_new, lty = 2, lwd = 2)
abline(v = threshold_hard_excess_new, lty = 4, lwd = 2)

plot(L_Hard, color_X, col = galaxy_class_X, pch = 16, ylim = c(-1, 2.5),
     xlab = "log[L (Hard X-Rays)]", ylab = "X-Rays color")
# abline(h = threshold_colorX_lit, lty = 2, lwd = 2)
abline(v = 42, lty = 4, lwd = 2)



###################################### COLOR EXCESS PLOTS ########################################

plot(L_reducida_IR1, L_reducida_IR2 - L_reducida_IR1, pch = 16, col = alpha('blue', 0.05), xlim = c(-2, 3), ylim = c(-0.3, 0.5),
    # main = expression(paste(Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])", " vs ", Delta, "(log[L(IR W1)])")),
     xlab = expression(paste(Delta, "(log[L(IR W1)])")), 
     ylab = expression(paste(Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])")), type = 'n')
for (i in 1:len){
  if ((is.na(L_reducida_IR2[i] - L_reducida_IR1[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess_new)){
    if (((L_reducida_IR2[i] - L_reducida_IR1[i]) > threshold_colorIR_excess_new)){
      points(L_reducida_IR1[i], L_reducida_IR2[i] - L_reducida_IR1[i], pch = 16, col = alpha('red4', 0.3))
    }
    else {
      points(L_reducida_IR1[i], L_reducida_IR2[i] - L_reducida_IR1[i], pch = 16, col = alpha('green', 0.3))
    }
  }
  else if ((is.na((L_reducida_IR2[i] - L_reducida_IR1[i])) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess_new)){
    if (((L_reducida_IR2[i] - L_reducida_IR1[i]) > threshold_colorIR_excess_new)){
      points(L_reducida_IR1[i], L_reducida_IR2[i] - L_reducida_IR1[i], pch = 16, col = alpha('gold', 0.07))
    }
    else {
      points(L_reducida_IR1[i], L_reducida_IR2[i] - L_reducida_IR1[i], pch = 16, col = alpha('blue', 0.02))
    }
  }
}
# abline(h = (moda_color_W2_W1_red + 2*sigma_izquierda_color_W2_W1_red), lty = 1, lwd = 2)
# abline(h = (moda_color_W2_W1_red - 2*sigma_izquierda_color_W2_W1_red), lty = 1, lwd = 2)
abline(h = threshold_colorIR_excess_new, lty = 1, lwd = 2)
abline(a = threshold_ir2_excess_new, b = -1, lty = 4, lwd = 2)



plot(L_reducida_IR, L_reducida_Radio - L_reducida_IR, pch = 16, col = alpha('blue', 0.1),
     # main = expression(paste(Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(IR W4)])", " vs ", Delta, "(log[L(IR W4)])")),
     xlab = expression(paste(Delta, "(log[L(W4)])")), 
     ylab = expression(paste(Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(W4)])")), type = 'n')
for (i in 1:len){
  if ((is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess_new)){
    if (((color_Radio_IR_reducido[i]) > threshold_color_Radio_IR_excess_new)){
      points(L_reducida_IR[i], L_reducida_Radio[i] - L_reducida_IR[i], pch = 16, col = alpha('red4', 0.1))
    }
    else {
      points(L_reducida_IR[i], L_reducida_Radio[i] - L_reducida_IR[i], pch = 16, col = alpha('green', 0.1))
    }
  }
  else if ((is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess_new)){
    if (((color_Radio_IR_reducido[i]) > threshold_color_Radio_IR_excess_new)){
      points(L_reducida_IR[i], L_reducida_Radio[i] - L_reducida_IR[i], pch = 16, col = alpha('gold', 0.1))
    }
    else {
      points(L_reducida_IR[i], L_reducida_Radio[i] - L_reducida_IR[i], pch = 16, col = alpha('blue', 0.1))
    }
  }
}
# abline(h = (moda_color_Radio_IR_red - 2*sigma_izquierda_color_Radio_IR_red), lty = 1, lwd = 2)
# abline(h = (moda_color_Radio_IR_red + 2*sigma_izquierda_color_Radio_IR_red), lty = 1, lwd = 2)
abline(h = threshold_color_Radio_IR_excess_new, lty = 2, lwd = 2)
abline(a = threshold_radio_excess_new, b = -1, lty = 4, lwd = 2)



plot(L_reducida_Soft, L_reducida_Hard - L_reducida_Soft, pch = 16, col = alpha('blue', 0.5), ylim = c(-2, 2.5),
    # main = expression(paste(Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])", " vs ", Delta, "(log[L(Soft X-Rays)])")),
     xlab = expression(paste(Delta, "(log[L(Soft X-Rays)])")), 
     ylab = expression(paste(Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])")), type = 'n')
for (i in 1:len){
  if ((is.na(L_reducida_Hard[i] - L_reducida_Soft[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess_new)){
    if (((L_reducida_Hard[i] - L_reducida_Soft[i]) > threshold_colorX_excess_new)){
      points(L_reducida_Soft[i], L_reducida_Hard[i] - L_reducida_Soft[i], pch = 16, col = alpha('red4', 0.3))
    }
    else {
      points(L_reducida_Soft[i], L_reducida_Hard[i] - L_reducida_Soft[i], pch = 16, col = alpha('green', 0.3))
    }
  }
  else if ((is.na(L_reducida_Hard[i] - L_reducida_Soft[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess_new)){
    if (((L_reducida_Hard[i] - L_reducida_Soft[i]) > threshold_colorX_excess_new)){
      points(L_reducida_Soft[i], L_reducida_Hard[i] - L_reducida_Soft[i], pch = 16, col = alpha('gold', 0.3))
    }
    else {
      points(L_reducida_Soft[i], L_reducida_Hard[i] - L_reducida_Soft[i], pch = 16, col = alpha('blue', 0.3))
    }
  }
}
# abline(h = moda_colorX_red + 2*sigma_izquierda_colorX_red, lty = 1, lwd = 2)
# abline(h = moda_colorX_red - 2*sigma_izquierda_colorX_red, lty = 1, lwd = 2)
abline(h = threshold_colorX_excess_new, lty =1, lwd = 2)
abline(a = threshold_hard_excess_new, b = -1, lty = 4, lwd = 2)



##################################### LUMINOSITY EXCESS PLOTS ###########################################

plot(L_IR1, L_reducida_IR2, pch = 16, col = alpha('blue', 0.05),
    # main = expression(paste(Delta, "(log[L(IR W2)])"," vs ", "log[L(IR W1)]")),
     xlab = expression(paste("log[L(IR W1)]")), 
     ylab = expression(paste(Delta, "(log[L(IR W2)])")), type = 'n')
for (i in 1:len){
  if ((is.na(L_reducida_IR2[i] - L_reducida_IR1[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess_new)){
    if (((L_reducida_IR2[i] - L_reducida_IR1[i]) > threshold_colorIR_excess_new)){
      points(L_IR1[i], L_reducida_IR2[i], pch = 16, col = alpha('red4', 0.3))
    }
    else {
      points(L_IR1[i], L_reducida_IR2[i], pch = 16, col = alpha('green', 0.3))
    }
  }
  else if ((is.na(L_reducida_IR2[i] - L_reducida_IR1[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess_new)){
    if (((L_reducida_IR2[i] - L_reducida_IR1[i]) > threshold_colorIR_excess_new)){
      points(L_IR1[i], L_reducida_IR2[i], pch = 16, col = alpha('gold', 0.2))
    }
    else {
      points(L_IR1[i], L_reducida_IR2[i], pch = 16, col = alpha('blue', 0.05))
    }
  }
}
# abline(h = (moda_ir2_excess + 2*sigma_izquierda_ir2_excess), lty = 1, lwd = 2)
# abline(h = (moda_ir2_excess - 2*sigma_izquierda_ir2_excess), lty = 1, lwd = 2)
abline(h = threshold_ir2_excess_new, lty = 1, lwd = 2)


plot(L_IR, L_reducida_Radio, pch = 16, col = alpha('blue', 0.1),
    # main = expression(paste(Delta, "(log[L(Radio)])"," vs ", "log[L(IR W4)]")),
     xlab = expression(paste("log[L(W4)]")), 
     ylab = expression(paste(Delta, "(log[L(Radio)])")), type = 'n')
for (i in 1:len){
  if ((is.na(L_reducida_Radio[i] - L_reducida_IR[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess_new)){
    if ((color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess_new)){
      points(L_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('red4', 0.1))
    }
    else {
      points(L_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('green', 0.1))
    }
  }
  else if ((is.na(L_reducida_Radio[i] - L_reducida_IR[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess_new)){
    if (((L_reducida_Radio[i] - L_reducida_IR[i]) > threshold_color_Radio_IR_excess_new)){
      points(L_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('gold', 0.1))
    }
    else {
      points(L_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('blue', 0.1))
    }
  }
}
# abline(h = (moda_radio_excess - 2*sigma_izquierda_radio_excess), lty = 1, lwd = 2)
# abline(h = (moda_radio_excess + 2*sigma_izquierda_radio_excess), lty = 1, lwd = 2)
abline(h = threshold_radio_excess_new, lty = 1, lwd = 2)


plot(L_Soft, L_reducida_Hard, pch = 16, col = alpha('blue', 0.5),
    # main = expression(paste(Delta, "(log[L(Hard X-Rays)])", " vs ", "log[L(Soft X-Rays)]")),
     xlab = expression(paste("log[L(Soft X-Rays)]")), 
     ylab = expression(paste(Delta, "(log[L(Hard X-Rays)])")), type = 'n')
for (i in 1:len){
  if ((is.na(L_reducida_Hard[i] - L_reducida_Soft[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess_new)){
    if (((L_reducida_Hard[i] - L_reducida_Soft[i]) > threshold_colorX_excess_new)){
      points(L_Soft[i], L_reducida_Hard[i], pch = 16, col = alpha('red4', 0.3))
    }
    else {
      points(L_Soft[i], L_reducida_Hard[i], pch = 16, col = alpha('green', 0.3))
    }
  }
  else if ((is.na(L_reducida_Hard[i] - L_reducida_Soft[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess_new)){
    if (((L_reducida_Hard[i] - L_reducida_Soft[i]) > threshold_colorX_excess_new)){
      points(L_Soft[i], L_reducida_Hard[i], pch = 16, col = alpha('gold', 0.3))
    }
    else {
      points(L_Soft[i], L_reducida_Hard[i], pch = 16, col = alpha('blue', 0.3))
    }
  }
}
# abline(h = moda_hard_excess + 2*sigma_izquierda_hard_excess, lty = 1, lwd = 2)
# abline(h = moda_hard_excess - 2*sigma_izquierda_hard_excess, lty = 1, lwd = 2)
abline(h = threshold_hard_excess_new, lty = 1, lwd = 2)



################ CONTADOR DE OBJETOS PARA NUESTROS CRITERIOS Y LOS DE LA LITERATURA ###############

# IR_count_2s <- 0
# IR_count_stern <-0
# 
# for (i in 1:len){
#   if (((is.na(color_W2_W1[i])) == FALSE) & (color_W2_W1[i] > threshold_colorIR_lit) == TRUE){
#     IR_count_stern <- IR_count_stern + 1
#   }
#   if (((is.na(color_W2_W1[i])) == FALSE) & (color_W2_W1[i] > threshold_colorIR_ours) == TRUE){
#     IR_count_2s <- IR_count_2s + 1
#   }
# }
# 
# 
# Radio_count_2s <- 0
# Radio_count_Ibar <- 0
# 
# for (i in 1:len){
#   if (((is.na(color_Radio_IR[i])) == FALSE) & (color_Radio_IR[i] > threshold_color_Radio_IR_lit) == TRUE){
#     Radio_count_Ibar <- Radio_count_Ibar + 1
#   }
#   if (((is.na(color_Radio_IR[i])) == FALSE) & (color_Radio_IR[i] > threshold_color_Radio_IR_ours) == TRUE){
#     Radio_count_2s <- Radio_count_2s + 1
#   }
# }
# 
# 
# Radio_count_2s <- 0
# Radio_count_40 <- 0
# Radio_count_mix_2s_40 <- 0
# 
# for (i in 1:len){
#   if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(L_IR[i]) == FALSE)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) == TRUE){
#       Radio_count_2s <- Radio_count_2s + 1
#     }
#     if ((color_Radio_IR[i] > (40 - L_IR[i])) == TRUE){
#       Radio_count_40 <- Radio_count_40 + 1
#     }
#     if (((color_Radio_IR[i] > (40 - L_IR[i])) == TRUE) & ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) == TRUE)){
#       Radio_count_mix_2s_40 <- Radio_count_mix_2s_40 + 1
#     }
#   }
# }
# 
# Radio_count_Ibar <- 0
# Radio_count_40 <- 0
# Radio_count_mix_Ibar_40 <- 0
# 
# for (i in 1:len){
#   if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(L_IR[i]) == FALSE)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_lit) == TRUE){
#       Radio_count_Ibar <- Radio_count_Ibar + 1
#     }
#     if ((color_Radio_IR[i] > (40 - L_IR[i])) == TRUE){
#       Radio_count_40 <- Radio_count_40 + 1
#     }
#     if (((color_Radio_IR[i] > (40 - L_IR[i])) == TRUE) & ((color_Radio_IR[i] > threshold_color_Radio_IR_lit) == TRUE)){
#       Radio_count_mix_Ibar_40 <- Radio_count_mix_Ibar_40 + 1
#     }
#   }
# }
# 
# X_count_2s <- 0
# X_count_42 <- 0
# X_count_mix <- 0
# 
# for (i in 1:len){
#   if ((is.na(color_X[i]) == FALSE) & (is.na(L_Soft[i]) == FALSE)){
#     if ((color_X[i] > threshold_colorX_ours) == TRUE){
#       X_count_2s <- X_count_2s + 1
#     }
#     if ((color_X[i] > (42 - L_Soft[i])) == TRUE){
#       X_count_42 <- X_count_42 + 1
#     }
#     if (((color_X[i] > (42 - L_Soft[i])) == TRUE) & ((color_X[i] > threshold_colorX_ours) == TRUE)){
#       X_count_mix <- X_count_mix + 1
#     }
#   }
# }


####################################### MATEOS ET AL (2012) PLOT #############################################

# W1_W2 <- sdss_xmatch$W1mag - sdss_xmatch$W2mag
# 
# W2_W3 <- sdss_xmatch$W2mag - sdss_xmatch$W3mag
# 
# 
# plot(W2_W3, W1_W2, pch = 16, col = alpha('blue', 0.05),
#      main = "W2-W3 color vs W1-W2 color",
#      xlab = "W2-W3", 
#      ylab = "W1-W2", type = 'n', ylim = c(-0.2, 2), xlim = c(0, 5))
# for (i in 1:len){
#   if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess)){
#     if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(W2_W3[i], W1_W2[i], pch = 16, col = alpha('red4', 0.3))
#     }
#     else {
#       points(W2_W3[i], W1_W2[i], pch = 16, col = alpha('green', 0.3))
#     }
#   }
#   else if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess)){
#     if ((color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(W2_W3[i], W1_W2[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(W2_W3[i], W1_W2[i], pch = 16, col = alpha('blue', 0.01))
#     }
#   }
# }
# # abline(a = -0.222, b = 0.315, lty = 1, lwd = 2)
# # abline(a = 0.796, b = 0.315, lty = 1, lwd = 2)
# # abline(a = 7.624, b = -3.172, lty = 1, lwd = 2)
# 
# segments(x0 = 2.250, y0 = 0.486, x1 = 5.200, y1 = 1.416, lty = 1, lwd = 2)
# segments(x0 = 1.958, y0 = 1.413, x1 = 5.200, y1 = 2.434, lty = 1, lwd = 2)
# segments(x0 = 2.250, y0 = 0.486, x1 = 1.958, y1 = 1.413, lty = 1, lwd = 2)
# abline(h = 0.8, lty = 2, lwd= 2)
# abline(h = 0.5, lty = 3, lwd= 2)



############### HISTOGRAMAS CON NUESTROS THRESHOLDS Y LOS DE LA LITERATURA ################

# hist(color_W2_W1, col = 'darkred', breaks = 250, xlim = c(-0.7, 0.2), freq = FALSE,
#      main = ' Histogram of IR color', xlab = 'IR color')
# abline(v = moda_color_W2_W1 + 2*sigma_izquierda_color_W2_W1, lty = 1, lwd = 2)
# abline(v = moda_color_W2_W1 - 2*sigma_izquierda_color_W2_W1, lty = 1, lwd = 2)
# abline(v = 0.4*0.8 + W1_W2_norm, lty = 2, lwd = 2)
# # abline(v = 0.4*(moda_color_IR + 2*sigma_izquierda_colorIR) + W1_W2_norm, lty = 1, lwd = 2)
# # abline(v = 0.4*(moda_color_IR - 2*sigma_izquierda_colorIR) + W1_W2_norm, lty = 1, lwd = 2)
# 
# hist(color_Radio_IR, col = 'darkgreen', freq = FALSE, breaks = 70,
#      main = 'Histogram of Radio-IR color', xlab = 'Radio-IR color')
# abline(v = 0.23 + nuradio_nu24, lty = 2, lwd = 2)
# abline(v = (moda_color_Radio_IR - 2*sigma_izquierda_color_Radio_IR), lty = 1, lwd = 2)
# abline(v = (moda_color_Radio_IR + 2*sigma_izquierda_color_Radio_IR), lty = 1, lwd = 2)
# 
# hist(color_X, col = 'darkblue', freq = FALSE, breaks = 40,
#      main = ' Histogram of X-Rays color', xlab = 'X-Rays color')
# # abline(v = log10(5/2), lty = 2, lwd = 2)
# abline(v = moda_colorX + 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
# abline(v = moda_colorX - 2*sigma_izquierda_colorX, lty = 1, lwd = 2)
# 
# 
# hist(color_W2_W1_reducido, col = 'darkred', breaks = 250, xlim = c(-0.2, 0.3), freq = FALSE,
#      main = expression(paste('Histogram of ', Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])")),
#      xlab = expression(paste(Delta, "(log[L(IR W2)])", " - ", Delta, "(log[L(IR W1)])")))
# abline(v = (moda_color_W2_W1_red + 2*sigma_izquierda_color_W2_W1_red), lty = 1, lwd = 2)
# abline(v = (moda_color_W2_W1_red - 2*sigma_izquierda_color_W2_W1_red), lty = 1, lwd = 2)
# 
# hist(color_Radio_IR_reducido, col = 'darkgreen', freq = FALSE, breaks = 70,
#      main = expression(paste('Histogram of ', Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(IR W4)])")),
#      xlab = expression(paste(Delta, "(log[L(Radio)])", " - ", Delta, "(log[L(IR W4)])")))
# abline(v = moda_color_Radio_IR_red - 2*sigma_izquierda_color_Radio_IR_red, lty = 1, lwd = 2)
# abline(v = moda_color_Radio_IR_red + 2*sigma_izquierda_color_Radio_IR_red, lty = 1, lwd = 2)
# 
# hist(color_X_reducido, col = 'darkblue', freq = FALSE, breaks = 50,
#      main = expression(paste('Histogram of ', Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])")),
#      xlab = expression(paste(Delta, "(log[L(Hard X-Rays)])", " - ", Delta, "(log[L(Soft X-Rays)])")))
# abline(v = moda_colorX_red + 2*sigma_izquierda_colorX_red, lty = 1, lwd = 2)
# abline(v = moda_colorX_red - 2*sigma_izquierda_colorX_red, lty = 1, lwd = 2)
# 
# 
# hist(L_reducida_IR2, col = 'darkred', freq = FALSE, breaks = 250, xlim = c(-1,2),
#      main = expression(paste("Histogram of"," ", Delta, "(log[L(IR W2 band)])")), 
#      xlab = expression(paste(Delta, "(log[L(IR W2 band)])")))
# abline(v = moda_ir2_excess + 2*sigma_izquierda_ir2_excess, lty = 1, lwd = 2)
# abline(v = moda_ir2_excess - 2*sigma_izquierda_ir2_excess, lty = 1, lwd = 2)
# 
# 
# hist(L_reducida_Radio, col = 'darkgreen', freq = FALSE, breaks = 70,
#      main = expression(paste("Histogram of"," ", Delta, "(log[L(Radio)])")), 
#      xlab = expression(paste(Delta, "(log[L(Radio)])")))
# abline(v = moda_radio_excess + 2*sigma_izquierda_radio_excess, lty = 1, lwd = 2)
# abline(v = moda_radio_excess - 2*sigma_izquierda_radio_excess, lty = 1, lwd = 2)
# 
# 
# hist(L_reducida_Hard, col = 'darkblue', freq = FALSE, breaks = 30,
#      main = expression(paste("Histogram of"," ", Delta, "(log[L(Hard X-Rays)])")), 
#      xlab = expression(paste(Delta, "(log[L(Hard X-Rays)])")))
# abline(v = moda_hard_excess + 2*sigma_izquierda_hard_excess, lty = 1, lwd = 2)
# abline(v = moda_hard_excess - 2*sigma_izquierda_hard_excess, lty = 1, lwd = 2)



#####################################################   Todos vs Todos (IR)  ###############################################################

# plot(color_W2_W1, L_reducida_IR2, pch = 16, col = alpha('blue', 0.05), xlim = c(-0.6, 0.2),
#      main = expression(paste(Delta, "(log[L(IR W2)])", " vs ", "IR color")),
#      xlab = 'IR color', 
#      ylab = expression(paste(Delta, "(log[L(IR W2)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) & (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1[i], L_reducida_IR2[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_W2_W1[i], L_reducida_IR2[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) | (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1[i], L_reducida_IR2[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_W2_W1[i], L_reducida_IR2[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }
# 
# 
# 
# plot(color_W2_W1_reducido, L_reducida_IR2, pch = 16, col = alpha('blue', 0.05),
#      main = expression(paste(Delta, "(log[L(IR W2)])", " vs ", Delta, "(IR color)")),
#      xlab = expression(paste(Delta, "(IR color)")), 
#      ylab = expression(paste(Delta, "(log[L(IR W2)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) & (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1_reducido[i], L_reducida_IR2[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_W2_W1_reducido[i], L_reducida_IR2[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) | (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1_reducido[i], L_reducida_IR2[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_W2_W1_reducido[i], L_reducida_IR2[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }
# 
# 
# 
# plot(color_W2_W1, color_W2_W1_reducido, pch = 16, col = alpha('blue', 0.05), xlim = c(-1, 0.5), ylim = c(-0.5, 0.75),
#      main = expression(paste(Delta, "(IR color)", " vs ", "IR color")),
#      xlab = 'IR color', 
#      ylab = expression(paste(Delta, "(IR color)"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] > threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) & (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1[i], color_W2_W1_reducido[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_W2_W1[i], color_W2_W1_reducido[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_W2_W1[i]) == FALSE) & (is.na(color_W2_W1_reducido[i]) == FALSE) & (is.na(L_reducida_IR2[i]) == FALSE) & (L_reducida_IR2[i] < threshold_ir2_excess)){
#     if ((color_W2_W1[i] > threshold_colorIR_ours) | (color_W2_W1_reducido[i] > threshold_colorIR_excess)){
#       points(color_W2_W1[i], color_W2_W1_reducido[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_W2_W1[i], color_W2_W1_reducido[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }


#####################################################   Todos vs Todos (Radio)  ###############################################################

# plot(color_Radio_IR, L_reducida_Radio, pch = 16, col = alpha('blue', 0.05),
#      main = expression(paste(Delta, "(log[L(Radio)])", " vs ", "Radio-IR color")),
#      xlab = 'Radio-IR color', 
#      ylab = expression(paste(Delta, "(log[L(Radio)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) & (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_Radio_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) | (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_Radio_IR[i], L_reducida_Radio[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }
# 
# 
# 
# plot(color_Radio_IR_reducido, L_reducida_Radio, pch = 16, col = alpha('blue', 0.05),
#      main = expression(paste(Delta, "(log[L(Radio)])", " vs ", Delta, "(Radio-IR color)")),
#      xlab = expression(paste(Delta, "(Radio-IR color)")), 
#      ylab = expression(paste(Delta, "(log[L(Radio)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) & (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR_reducido[i], L_reducida_Radio[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_Radio_IR_reducido[i], L_reducida_Radio[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) | (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR_reducido[i], L_reducida_Radio[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_Radio_IR_reducido[i], L_reducida_Radio[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }
# 
# 
# plot(color_Radio_IR, color_Radio_IR_reducido, pch = 16, col = alpha('blue', 0.05),
#      main = expression(paste(Delta, "(Radio-IR color)", " vs ", "Radio-IR color")),
#      xlab = 'Radio-IR color', 
#      ylab = expression(paste(Delta, "(Radio-IR color)"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] > threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) & (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR[i], color_Radio_IR_reducido[i], pch = 16, col = alpha('red4', 0.05))
#     }
#     else {
#       points(color_Radio_IR[i], color_Radio_IR_reducido[i], pch = 16, col = alpha('green', 0.05))
#     }
#   }
#   else if ((is.na(color_Radio_IR[i]) == FALSE) & (is.na(color_Radio_IR_reducido[i]) == FALSE) & (is.na(L_reducida_Radio[i]) == FALSE) & (L_reducida_Radio[i] < threshold_radio_excess)){
#     if ((color_Radio_IR[i] > threshold_color_Radio_IR_ours) | (color_Radio_IR_reducido[i] > threshold_color_Radio_IR_excess)){
#       points(color_Radio_IR[i], color_Radio_IR_reducido[i], pch = 16, col = alpha('gold', 0.05))
#     }
#     else {
#       points(color_Radio_IR[i], color_Radio_IR_reducido[i], pch = 16, col = alpha('blue', 0.02))
#     }
#   }
# }


#####################################################   Todos vs Todos (X-Rays)  ###############################################################

# plot(color_X, L_reducida_Hard, pch = 16, col = alpha('blue', 0.05), xlim = c(-1, 3),
#      main = expression(paste(Delta, "(log[L(Hard X-Rays)])", " vs ", "X-Rays color")),
#      xlab = 'X-Rays color', 
#      ylab = expression(paste(Delta, "(log[L(Hard X-Rays)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) & (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X[i], L_reducida_Hard[i], pch = 16, col = alpha('red4', 0.3))
#     }
#     else {
#       points(color_X[i], L_reducida_Hard[i], pch = 16, col = alpha('green', 0.3))
#     }
#   }
#   else if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) | (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X[i], L_reducida_Hard[i], pch = 16, col = alpha('gold', 0.3))
#     }
#     else {
#       points(color_X[i], L_reducida_Hard[i], pch = 16, col = alpha('blue', 0.3))
#     }
#   }
# }
# 
# 
# 
# plot(color_X_reducido, L_reducida_Hard, pch = 16, col = alpha('blue', 0.05), xlim = c(-1.5, 3),
#      main = expression(paste(Delta, "(log[L(Hard X-Rays)])", " vs ", Delta, "(X-Rays color)")),
#      xlab = expression(paste(Delta, "(X_Rays color)")), 
#      ylab = expression(paste(Delta, "(log[L(Hard X-Rays)])"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) & (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X_reducido[i], L_reducida_Hard[i], pch = 16, col = alpha('red4', 0.3))
#     }
#     else {
#       points(color_X_reducido[i], L_reducida_Hard[i], pch = 16, col = alpha('green', 0.3))
#     }
#   }
#   else if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) | (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X_reducido[i], L_reducida_Hard[i], pch = 16, col = alpha('gold', 0.3))
#     }
#     else {
#       points(color_X_reducido[i], L_reducida_Hard[i], pch = 16, col = alpha('blue', 0.3))
#     }
#   }
# }
# 
# 
# 
# 
# plot(color_X, color_X_reducido, pch = 16, col = alpha('blue', 0.05), xlim = c(-1,3), ylim = c(-1.5, 3),
#      main = expression(paste(Delta, "(X-Rays color)", " vs ", "X-Rays color")),
#      xlab = "(X-Rays color)", 
#      ylab = expression(paste(Delta, "(X-Rays color)"), type = 'n'))
# for (i in 1:len){
#   if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] > threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) & (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X[i], color_X_reducido[i], pch = 16, col = alpha('red4', 0.3))
#     }
#     else {
#       points(color_X[i], color_X_reducido[i], pch = 16, col = alpha('green', 0.3))
#     }
#   }
#   else if ((is.na(color_X[i]) == FALSE) & (is.na(color_X_reducido[i]) == FALSE) & (is.na(L_reducida_Hard[i]) == FALSE) & (L_reducida_Hard[i] < threshold_hard_excess)){
#     if ((color_X[i] > threshold_colorX_ours) | (color_X_reducido[i] > threshold_colorX_excess)){
#       points(color_X[i], color_X_reducido[i], pch = 16, col = alpha('gold', 0.3))
#     }
#     else {
#       points(color_X[i], color_X_reducido[i], pch = 16, col = alpha('blue', 0.3))
#     }
#   }
# }


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



