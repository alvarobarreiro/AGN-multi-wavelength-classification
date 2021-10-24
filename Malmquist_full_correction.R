

########## CARGA DE LOS DATOS ##########

setwd('D:/AGN/data/final_plots') # directorio de trabajo
sdss_xmatch <- read.csv('final_xmatch.csv', header = TRUE) # Archivo con el cross-match entre SDSS y los demás catálogos
len <- length(sdss_xmatch$X) # longitud de los datos
redshift <- sdss_xmatch$Z 

SFR <- sdss_xmatch$sfr_tot_p50_1 # SFR (en logarítmico)
M <- sdss_xmatch$lgm_tot_p50 # Masa estelar de la galaxia (en logarítmico)
Delta_SFR = (sdss_xmatch$sfr_tot_p84_1 - sdss_xmatch$sfr_tot_p16_1)/2 # Sigma para la SFR
Delta_M = (sdss_xmatch$lgm_tot_p84 - sdss_xmatch$lgm_tot_p16)/2 # Sigma para la M*


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
  if ((is.na(SFR[i]) == FALSE & SFR[i] == -9999) | (is.na(M[i]) == FALSE & M[i] == -9999)){ 
    
    SFR[i] <- NA
    M[i] <- NA
    Delta_M[i] <- NA
    Delta_SFR[i] <- NA
    
    L_Radio[i] <- NA
    L_IR[i] <- NA
    L_IR1[i] <- NA
    L_IR2[i] <- NA
    L_Soft[i] <- NA
    L_Hard[i] <- NA
    color_W2_W1[i] <- NA
    color_Radio_IR[i] <- NA
    color_X[i] <- NA
    
    color_IR[i] <- NA
    q_24[i] <- NA
    
  }
  else if (((is.na(M[i]/Delta_M[i]) == FALSE) & ((M[i]/Delta_M[i]) < 2)) | ((is.na(SFR[i]/Delta_SFR[i]) == FALSE) & ((SFR[i]/Delta_SFR[i]) < 2))){ 
    
    SFR[i] <- NA
    M[i] <- NA
    Delta_M[i] <- NA
    Delta_SFR[i] <- NA
    
    L_Radio[i] <- NA
    L_IR[i] <- NA
    L_IR1[i] <- NA
    L_IR2[i] <- NA
    L_Soft[i] <- NA
    L_Hard[i] <- NA
    color_W2_W1[i] <- NA
    color_Radio_IR[i] <- NA
    color_X[i] <- NA
    
    color_IR[i] <- NA
    q_24[i] <- NA
    
  }
  if (redshift[i] < 0.01){
    
    SFR[i] <- NA
    M[i] <- NA
    Delta_M[i] <- NA
    Delta_SFR[i] <- NA
    
    L_Radio[i] <- NA
    L_IR[i] <- NA
    L_IR1[i] <- NA
    L_IR2[i] <- NA
    L_Soft[i] <- NA
    L_Hard[i] <- NA
    color_W2_W1[i] <- NA
    color_Radio_IR[i] <- NA
    color_X[i] <- NA
    
    color_IR[i] <- NA
    q_24[i] <- NA
    
  }
}

sigma_SoftX <- sdss_xmatch$Soft_flux_error # Noise
SN_SoftX <- sdss_xmatch$Soft_flux/sigma_SoftX # Signal-to-Noise ratio

sigma_HardX <- sdss_xmatch$Hard_flux_error 
SN_HardX <- sdss_xmatch$Hard_flux/sigma_HardX 

sigma_Radio <- sdss_xmatch$e_S1_4 
SN_Radio <- sdss_xmatch$S1_4/sigma_Radio 

sigma_IR4 <- sdss_xmatch$IR_flux_error..erg.s.cm.2.Hz. 
SN_IR4 <- sdss_xmatch$IR_flux..erg.s.cm.2.Hz./sigma_IR4 

sigma_IR1 <- sdss_xmatch$IR1_flux_error_.erg.s.cm2.Hz. 
SN_IR1 <- sdss_xmatch$IR1_flux_.erg.s.cm2.Hz./sigma_IR1 

sigma_IR2 <- sdss_xmatch$IR2_flux_error_.erg.s.cm2.Hz. 
SN_IR2 <- sdss_xmatch$IR2_flux_.erg.s.cm2.Hz./sigma_IR2 


# Filtro en flujo en rayos X blandos:
for (i in 1:len){ 
  if (is.na(SN_SoftX[i]) == TRUE){ 
    L_Soft[i] <- NA
  }
  else{
    if (SN_SoftX[i] < 2){
      L_Soft[i] <- NA
    }
  }
}

# Filtro en flujo en rayos X duros:
for (i in 1:len){ 
  if (is.na(SN_HardX[i]) == TRUE){ 
    L_Hard[i] <- NA
  }
  else{
    if (SN_HardX[i] < 2){
      L_Hard[i] <- NA
    }
  }
}

# Filtro en flujo en Radio:
for (i in 1:len){ 
  if (is.na(SN_Radio[i]) == TRUE){ 
    L_Radio[i] <- NA
  }
  else{
    if (SN_Radio[i] < 2){
      L_Radio[i] <- NA
    }
  }
}

# Filtro en flujo en IR (banda W4):
for (i in 1:len){ 
  if (is.na(SN_IR4[i]) == TRUE){ 
    L_IR[i] <- NA
  }
  else{
    if (SN_IR4[i] < 2){
      L_IR[i] <- NA
    }
  }
}

# Filtro en flujo en IR (banda W1):
for (i in 1:len){ 
  if (is.na(SN_IR1[i]) == TRUE){ 
    L_IR1[i] <- NA
  }
  else{
    if (SN_IR1[i] < 2){
      L_IR1[i] <- NA
    }
  }
}

# Filtro en flujo en IR (banda W2):
for (i in 1:len){ 
  if (is.na(SN_IR2[i]) == TRUE){ 
    L_IR2[i] <- NA
  }
  else{
    if (SN_IR2[i] < 2){
      L_IR2[i] <- NA
    }
  }
}


# Filtro en flujo en colorW2_W1:

for (i in 1:len){ 
  if (is.na(SN_IR1[i]) == TRUE | is.na(SN_IR2[i]) == TRUE){ 
    color_W2_W1[i] <- NA
    color_IR[i]
  }
  else if ((is.na(SN_IR1[i]) == FALSE & SN_IR1[i] < 2) | (is.na(SN_IR2[i]) == FALSE & SN_IR2[i] < 2)){
    color_W2_W1[i] <- NA
    color_IR[i] <- NA
  }
}

# Filtro en flujo en color_Radio_IR:

for (i in 1:len){ 
  if (is.na(SN_Radio[i]) == TRUE | is.na(SN_IR4[i]) == TRUE){ 
    color_Radio_IR[i] <- NA
    q_24[i] <- NA
  }
  else if ((is.na(SN_Radio[i]) == FALSE & SN_Radio[i] < 2) | (is.na(SN_IR4[i]) == FALSE & SN_IR4[i] < 2)){
    color_Radio_IR[i] <- NA
    q_24[i] <- NA
  }
}

# Filtro en flujo en color_X:

for (i in 1:len){ 
  if (is.na(SN_HardX[i]) == TRUE | is.na(SN_SoftX[i]) == TRUE){ 
    color_X[i] <- NA
  }
  else if ((is.na(SN_HardX[i]) == FALSE & SN_HardX[i] < 2) | (is.na(SN_SoftX[i]) == FALSE & SN_SoftX[i] < 2)){
    color_X[i] <- NA
  }
}


##########################   REMOVE EXTREME MASSES AND SFRs   ####################################


for (i in 1:len){
  if (is.na(M[i]) == FALSE & is.na(SFR[i]) == FALSE){
    
    if (M[i] < 8.5 | M[i] > 11.5 | SFR[i] < -2 | SFR[i] > 1){
      
      L_Radio[i] <- NA
      L_IR[i] <- NA
      L_IR1[i] <- NA
      L_IR2[i] <- NA
      L_Soft[i] <- NA
      L_Hard[i] <- NA
      color_W2_W1[i] <- NA
      color_Radio_IR[i] <- NA
      color_X[i] <- NA
      
      color_IR[i] <- NA
      q_24[i] <- NA
      
    }
  }
}


############# GRIDS Y FUNCIONES PARA EL CÁLCULO DE LOS KERNEL ESTIMATORS ##############

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



get_Vmax_inverse = function(L, sensitivity){
  z_max = 10^((L - sensitivity)/2) 
  z_max = apply(z_max, 2, min)
  z_max = ifelse(z_max > 0.01, z_max, 0.01)
  z_max = ifelse(z_max < 0.07, z_max, 0.07)
  return(1/z_max^3)
}

# Sensitivity (i.e. flux limit above which we have data for each band): it's been visually estimated from
# plot(redshift, L - 2*log10(redshift)) for each frequency

sensitivity_Radio = 40.85
sensitivity_IR = 44.8
sensitivity_IR1 = 44.6
sensitivity_IR2 = 44.2
sensitivity_Soft = 42.8
sensitivity_Hard = 43.5
sensitivity_SDSS = 12 # plot(redshift, M - 2*log10(redshift)) in this case



matrix_Radio = matrix(c(L_Radio, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_Radio_vector = c(sensitivity_Radio, sensitivity_SDSS)
V_max_Radio_inverse = get_Vmax_inverse(matrix_Radio, sens_Radio_vector)

matrix_IR = matrix(c(L_IR, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_IR_vector = c(sensitivity_IR, sensitivity_SDSS)
V_max_IR_inverse = get_Vmax_inverse(matrix_IR, sens_IR_vector)

matrix_IR1 = matrix(c(L_IR1, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_IR1_vector = c(sensitivity_IR1, sensitivity_SDSS)
V_max_IR1_inverse = get_Vmax_inverse(matrix_IR1, sens_IR1_vector)

matrix_IR2 = matrix(c(L_IR2, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_IR2_vector = c(sensitivity_IR2, sensitivity_SDSS)
V_max_IR2_inverse = get_Vmax_inverse(matrix_IR2, sens_IR2_vector)

matrix_Soft = matrix(c(L_Soft, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_Soft_vector = c(sensitivity_Soft, sensitivity_SDSS)
V_max_Soft_inverse = get_Vmax_inverse(matrix_Soft, sens_Soft_vector)

matrix_Hard = matrix(c(L_Hard, M), nrow = 2, ncol = length(M), byrow = TRUE)
sens_Hard_vector = c(sensitivity_Hard, sensitivity_SDSS)
V_max_Hard_inverse = get_Vmax_inverse(matrix_Hard, sens_Hard_vector)



matrix_color_IR = matrix(c(L_IR1, L_IR2, M), nrow = 3, ncol = length(M), byrow = TRUE)
sens_color_IR_vector = c(sensitivity_IR1, sensitivity_IR2, sensitivity_SDSS)
V_max_color_IR_inverse = get_Vmax_inverse(matrix_color_IR, sens_color_IR_vector)

matrix_color_Radio = matrix(c(L_Radio, L_IR, M), nrow = 3, ncol = length(M), byrow = TRUE)
sens_color_Radio_vector = c(sensitivity_Radio, sensitivity_IR, sensitivity_SDSS)
V_max_color_Radio_inverse = get_Vmax_inverse(matrix_color_Radio, sens_color_Radio_vector)

matrix_color_X = matrix(c(L_Hard, L_Soft, M), nrow = 3, ncol = length(M), byrow = TRUE)
sens_color_X_vector = c(sensitivity_Hard, sensitivity_Soft, sensitivity_SDSS)
V_max_color_X_inverse = get_Vmax_inverse(matrix_color_X, sens_color_X_vector)


# Getting densities, averages and dispersions for each luminosity and color:


get_outputs = function(V_max_inverse, L, sensitivity, peso_total, densidad_pesada, densidad_pesada_2){
  
  suavizado_minimo = 0.1
  peso_total <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  densidad_pesada <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  densidad_pesada_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  
  for (p in 1:len){
    L_p = L[p]
    
    if ((is.na(L_p) == FALSE) & (L_p - 2*log10(redshift[p])) > sensitivity){
      densidad_p = densidad(SFR_i = SFR_grid, M_j = M_grid, SFR_p = SFR[p], M_p = M[p],
                            Delta_SFR_p = suavizado_minimo + Delta_SFR[p], Delta_M_p = suavizado_minimo + Delta_M[p])
      w_p = densidad_p*V_max_inverse[p]
      peso_total = peso_total + w_p
      densidad_pesada = densidad_pesada + L_p*w_p
      densidad_pesada_2 = densidad_pesada_2 + L_p^2*w_p
    }
  }
  L_mean = densidad_pesada/peso_total
  L_mean_squared = densidad_pesada_2/peso_total
  Dispersion = sqrt(L_mean_squared - L_mean^2)
  
  outputs = list("Density" = peso_total, "L_mean" = L_mean, "Dispersion" = Dispersion)
  return(outputs)
  
}


get_color_outputs = function(V_max_inverse, L1, L2, color, sensitivity1, sensitivity2,
                                peso_total, densidad_pesada, densidad_pesada_2){
  
  suavizado_minimo = 0.1
  peso_total <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  densidad_pesada <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  densidad_pesada_2 <- matrix(0, nrow = 250, ncol = 250, byrow = TRUE)
  
  for (p in 1:len){
    L1_p = L1[p]
    L2_p = L2[p]
    color_p = color[p]
    
    if (((is.na(L1_p) == FALSE) & (L1_p - 2*log10(redshift[p])) > sensitivity1) & 
        ((is.na(L2_p) == FALSE) & (L2_p - 2*log10(redshift[p])) > sensitivity2)){
      
      densidad_p = densidad(SFR_i = SFR_grid, M_j = M_grid, SFR_p = SFR[p], M_p = M[p],
                            Delta_SFR_p = suavizado_minimo + Delta_SFR[p], Delta_M_p = suavizado_minimo + Delta_M[p])
      w_p = densidad_p*V_max_inverse[p]
      peso_total = peso_total + w_p
      densidad_pesada = densidad_pesada + color_p*w_p
      densidad_pesada_2 = densidad_pesada_2 + color_p^2*w_p
      
    }
  }
  
  color_mean = densidad_pesada/peso_total
  color_mean_squared = densidad_pesada_2/peso_total
  Dispersion = sqrt(color_mean_squared - color_mean^2)
  
  outputs = list("Density" = peso_total, "color_mean" = color_mean, "Dispersion" = Dispersion)
  return(outputs)
  
}


# Getting results from the previous functions:

out_Radio = get_outputs(V_max_Radio_inverse, L_Radio, sensitivity_Radio, peso_total, densidad_pesada, densidad_pesada_2)
out_IR = get_outputs(V_max_IR_inverse, L_IR, sensitivity_IR, peso_total, densidad_pesada, densidad_pesada_2)
out_IR1 = get_outputs(V_max_IR1_inverse, L_IR1, sensitivity_IR1, peso_total, densidad_pesada, densidad_pesada_2)
out_IR2 = get_outputs(V_max_IR2_inverse, L_IR2, sensitivity_IR2, peso_total, densidad_pesada, densidad_pesada_2)
out_Soft = get_outputs(V_max_Soft_inverse, L_Soft, sensitivity_Soft, peso_total, densidad_pesada, densidad_pesada_2)
out_Hard = get_outputs(V_max_Hard_inverse, L_Hard, sensitivity_Hard, peso_total, densidad_pesada, densidad_pesada_2)


out_color_W2_W1 = get_color_outputs(V_max_color_IR_inverse, L_IR1, L_IR2, color_W2_W1, 
                                  sensitivity_IR1, sensitivity_IR2, peso_total, densidad_pesada, densidad_pesada_2)

out_color_Radio_IR = get_color_outputs(V_max_color_Radio_inverse, L_Radio, L_IR, color_Radio_IR, 
                                    sensitivity_Radio, sensitivity_IR, peso_total, densidad_pesada, densidad_pesada_2)

out_color_X = get_color_outputs(V_max_color_X_inverse, L_Soft, L_Hard, color_X, 
                                    sensitivity_Soft, sensitivity_Hard, peso_total, densidad_pesada, densidad_pesada_2)


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
  fraction_inside = 1 - fraction # dará la fracción de masa DENTRO del contorno
  
  image.plot(M_X, SFR_y, Density, col = matlab.like(1000),
             xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Kernel Estimator: Density')
  contour(M_X, SFR_y, fraction_inside, levels = .9, add = TRUE, lty = 2)
  contour(M_X, SFR_y, fraction_inside, levels = .5, add = TRUE, lty = 3)
  contour(M_X, SFR_y, fraction_inside, levels = .1, add = TRUE, lty = 1)
  
  image.plot(M_X, SFR_y, fraction, col = matlab.like(1000),
             xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = 'Kernel Estimator: Fraction')
  contour(M_X, SFR_y, fraction_inside, levels = .9, add = TRUE, lty = 2)
  contour(M_X, SFR_y, fraction_inside, levels = .5, add = TRUE, lty = 3)
  contour(M_X, SFR_y, fraction_inside, levels = .1, add = TRUE, lty = 1)
  
  return(fraction_inside)
}

contours_IR1 = get_isocontour_levels(out_IR1$Density, d_M, d_SFR)
contours_IR2 = get_isocontour_levels(out_IR2$Density, d_M, d_SFR)
contours_IR4 = get_isocontour_levels(out_IR$Density, d_M, d_SFR)
contours_Radio = get_isocontour_levels(out_Radio$Density, d_M, d_SFR)
contours_Soft = get_isocontour_levels(out_Soft$Density, d_M, d_SFR)
contours_Hard = get_isocontour_levels(out_Hard$Density, d_M, d_SFR)

contours_color_W2_W1 = get_isocontour_levels(out_color_W2_W1$Density, d_M, d_SFR)
contours_color_Radio_IR = get_isocontour_levels(out_color_Radio_IR$Density, d_M, d_SFR)
contours_color_X = get_isocontour_levels(out_color_X$Density, d_M, d_SFR)



############ HISTOGRAMAS 2D DE LOS NUEVOS <L> Y <COLOR> Y DE SUS RESPECTIVAS DISPERSIONES ################

get_plots = function(L, contours, title){
  
  image.plot(M_X, SFR_y, L, col = matlab.like(15), ylim = c(-2, 1), xlim = c(8.5, 11.5),
             xlab = 'log[M* (solar masses)]', ylab = 'log[SFR (solar masses/yr)]', main = title)
  contour(M_X, SFR_y, contours, levels = .9, add = TRUE, lty = 3)
  contour(M_X, SFR_y, contours, levels = .5, add = TRUE, lty = 2)
  contour(M_X, SFR_y, contours, levels = .1, add = TRUE, lty = 1)
  
}

plot_Radio = get_plots(out_Radio$L_mean, contours_Radio, '<L> Radio')
plot_Radio_disp = get_plots(out_Radio$Dispersion, contours_Radio, 'Dispersion Radio')

plot_IR = get_plots(out_IR$L_mean, contours_IR4, '<L> IR (W4 band)')
plot_IR_disp = get_plots(out_IR$Dispersion, contours_IR4, 'Dispersion IR (W4 band)')

plot_IR1 = get_plots(out_IR1$L_mean, contours_IR1, '<L> IR (W1 band)')
plot_IR1_disp = get_plots(out_IR1$Dispersion, contours_IR1, 'Dispersion IR (W1 band)')

plot_IR2 = get_plots(out_IR2$L_mean, contours_IR2, '<L> IR (W2 band)')
plot_IR2_disp = get_plots(out_IR2$Dispersion, contours_IR2, 'Dispersion IR (W2 band)')

plot_Soft = get_plots(out_Soft$L_mean, contours_Soft, '<L> Soft X-Rays')
plot_Soft_disp = get_plots(out_Soft$Dispersion, contours_Soft, 'Dispersion Soft X-Rays')

plot_Hard = get_plots(out_Hard$L_mean, contours_Hard, '<L> Hard X-Rays')
plot_Hard_disp = get_plots(out_Hard$Dispersion, contours_Hard, 'Dispersion Hard X-Rays')



plot_color_W2_W1 = get_plots(out_color_W2_W1$color_mean, contours_color_W2_W1, '<IR color>')
plot_color_W2_W1_disp = get_plots(out_color_W2_W1$Dispersion, contours_color_W2_W1, 'Dispersion IR color')

plot_color_Radio_IR = get_plots(out_color_Radio_IR$color_mean, contours_color_Radio_IR, '<Radio - IR (W4 band) color>')
plot_color_Radio_IR_disp = get_plots(out_color_Radio_IR$Dispersion, contours_color_Radio_IR, 'Dispersion Radio - IR(W4) color')

plot_color_X = get_plots(out_color_X$color_mean, contours_color_X, '<X-Rays color>')
plot_color_X_disp = get_plots(out_color_X$Dispersion, contours_color_X, 'Dispersion X-Rays color')



############ CÁLCULO DE LAS LUMINOSIDADES Y COLORES MEDIOS Y DE SUS EXCESOS ############

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


L_media_soft <- L_mean_vector(L = L_Soft, L_media = L_media_soft, L_mean = out_Soft$L_mean)
L_media_hard <- L_mean_vector(L = L_Hard, L_media = L_media_hard, L_mean = out_Hard$L_mean)
L_media_radio <- L_mean_vector(L = L_Radio, L_media = L_media_radio, L_mean = out_Radio$L_mean)
L_media_ir <- L_mean_vector(L = L_IR, L_media = L_media_ir, L_mean = out_IR$L_mean)
L_media_ir1 <- L_mean_vector(L = L_IR1, L_media = L_media_ir1, L_mean = out_IR1$L_mean)
L_media_ir2 <- L_mean_vector(L = L_IR2, L_media = L_media_ir2, L_mean = out_IR2$L_mean)
color_W2_W1_media <- L_mean_vector(L = color_W2_W1, L_media = color_W2_W1_media, L_mean = out_color_W2_W1$color_mean)
color_Radio_IR_media <- L_mean_vector(L = color_Radio_IR, L_media = color_Radio_IR_media, L_mean = out_color_Radio_IR$color_mean)
colorX_media <- L_mean_vector(L = color_X, L_media = colorX_media, L_mean = out_color_X$color_mean)

L_excess_Soft <- L_Soft - L_media_soft
L_excess_Hard <- L_Hard - L_media_hard
L_excess_Radio <- L_Radio - L_media_radio
L_excess_IR <- L_IR - L_media_ir
L_excess_IR1 <- L_IR1 - L_media_ir1
L_excess_IR2 <- L_IR2 - L_media_ir2
color_W2_W1_excess <- color_W2_W1 - color_W2_W1_media
color_Radio_IR_excess <- color_Radio_IR - color_Radio_IR_media
color_X_excess <- color_X - colorX_media



############## CÁLCULO DE LOS PARÁMETROS DE CONVERSIÓN ENTRE NUESTROS COLORES Y LOS DE LA LITERATURA ###############

W1_W2_norm <- log10((171*3.4*10^-6)/(309*4.6*10^-6)) # cociente entre flujos cero de WISE
nuradio_nu24 <- log10((1.4*10^9)/((3*10^8)/(22.1*10^-6))) # cociente de frecuencias (1.4 GHz/22.1 micras)

################################## ACTIVENESS PLOTS #####################################

L_excess_IR2_copy = L_excess_IR2
V_max_IR2_inverse_copy = V_max_IR2_inverse
non_valid_IR2 = which(is.na(L_excess_IR2)==TRUE | is.na(V_max_IR2_inverse)==TRUE)
L_excess_IR2_copy[non_valid_IR2] = NA
V_max_IR2_inverse_copy[non_valid_IR2] = NA

L_excess_Radio_copy = L_excess_Radio
V_max_Radio_inverse_copy = V_max_Radio_inverse
non_valid_Radio = which(is.na(L_excess_Radio)==TRUE | is.na(V_max_Radio_inverse)==TRUE)
L_excess_Radio_copy[non_valid_Radio] = NA
V_max_Radio_inverse_copy[non_valid_Radio] = NA

L_excess_Hard_copy = L_excess_Hard
V_max_Hard_inverse_copy = V_max_Hard_inverse
non_valid_Hard = which(is.na(L_excess_Hard)==TRUE | is.na(V_max_Hard_inverse)==TRUE)
L_excess_Hard_copy[non_valid_Hard] = NA
V_max_Hard_inverse_copy[non_valid_Hard] = NA

color_W2_W1_excess_copy = color_W2_W1_excess
V_max_color_IR_inverse_copy = V_max_color_IR_inverse
non_valid_color_IR = which(is.na(color_W2_W1_excess)==TRUE | is.na(V_max_color_IR_inverse)==TRUE)
color_W2_W1_excess_copy[non_valid_color_IR] = NA
V_max_color_IR_inverse_copy[non_valid_color_IR] = NA

color_Radio_IR_excess_copy = color_Radio_IR_excess
V_max_color_Radio_inverse_copy = V_max_color_Radio_inverse
non_valid_color_Radio_IR = which(is.na(color_Radio_IR_excess)==TRUE | is.na(V_max_color_Radio_inverse)==TRUE)
color_Radio_IR_excess_copy[non_valid_color_Radio_IR] = NA
V_max_color_Radio_inverse_copy[non_valid_color_Radio_IR] = NA

color_X_excess_copy = color_X_excess
V_max_color_X_inverse_copy = V_max_color_X_inverse
non_valid_color_X = which(is.na(color_X_excess)==TRUE | is.na(V_max_color_X_inverse)==TRUE)
color_X_excess_copy[non_valid_color_X] = NA
V_max_color_X_inverse_copy[non_valid_color_X] = NA


# sigmas = function(lum, V_max_inverse){
#   
#   x_ordenado = na.omit(sort.int(lum, na.last = TRUE, index.return = TRUE)$x)
#   ranking = na.omit(sort.int(V_max_inverse, na.last = TRUE, index.return = TRUE)$x)
#   ranking = cumsum(ranking)
#   ranking_interp = approxfun(x_ordenado, ranking)
#   
#   ranking_interp2 = approxfun(ranking, x_ordenado)
#   total_dens = max(ranking)
#   mid_dens = total_dens/2
#   pivot_point = ranking_interp2(mid_dens)
#   
#   rho = total_dens
#   z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
#   rho_greater = rho - ranking_interp(z)
#   rho_normal = rep(0, length(rho_greater))
#   #pivot_point = median(x_ordenado)
#   
#   for (i in 1:length(z)){
#     if (z[i] <= pivot_point){
#       rho_symmetric = 2*ranking_interp(pivot_point) - ranking_interp(z[i])
#     }
#     else {
#       rho_symmetric = ranking_interp(2*pivot_point - z[i])
#     }
#     rho_normal[i] = min(rho_symmetric, rho_greater[i])
#   }
#   
#   rho_active = rho_greater - rho_normal
#   f_active = rho_active/rho_greater
#   z_active = z[min(which(f_active >= 0.5))]
#   n_sigmas = (z_active - quantile(na.omit(lum), probs = 0.50))/
#     (quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16))
#   
#   return(as.numeric(n_sigmas))
# }
# 
# 
# sigma_minima = min(sigmas(L_excess_IR2_copy, V_max_IR2_inverse_copy), sigmas(L_excess_Radio_copy, V_max_Radio_inverse_copy),
#                    sigmas(L_excess_Hard_copy, V_max_Hard_inverse_copy), sigmas(color_W2_W1_excess_copy, V_max_color_IR_inverse_copy), 
#                    sigmas(color_Radio_IR_excess_copy, V_max_color_Radio_inverse_copy), sigmas(color_X_excess_copy, V_max_color_X_inverse_copy))
#  
# sigma_maxima = max(sigmas(L_excess_IR2_copy, V_max_IR2_inverse_copy), sigmas(L_excess_Radio_copy, V_max_Radio_inverse_copy),
#                    sigmas(L_excess_Hard_copy, V_max_Hard_inverse_copy), sigmas(color_W2_W1_excess_copy, V_max_color_IR_inverse_copy), 
#                    sigmas(color_Radio_IR_excess_copy, V_max_color_Radio_inverse_copy), sigmas(color_X_excess_copy, V_max_color_X_inverse_copy))
# 
# 
# activeness = function(lum, xlabel, colour, V_max_inverse){
# 
#   indice_ordenado = na.omit(sort.int(lum, na.last = TRUE, index.return = TRUE)$ix)
#   x_ordenado = lum[indice_ordenado]
#   #ranking = na.omit(sort.int(V_max_inverse, na.last = TRUE, index.return = TRUE)$x)
#   ranking = V_max_inverse[indice_ordenado]
#   ranking = cumsum(ranking)
#   ranking = ranking/ranking[length(ranking)]
#   ranking_interp = approxfun(x_ordenado, ranking)
#   
#   ranking_interp2 = approxfun(ranking, x_ordenado)
#   total_dens = max(ranking)
#   mid_dens = total_dens/2
#   pivot_point = ranking_interp2(mid_dens)
#   
#   rho = total_dens
#   z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
#   rho_greater = rho - ranking_interp(z)
#   rho_normal = rep(0, length(rho_greater))
#   #pivot_point = median(x_ordenado)
#   
#   for (i in 1:length(z)){
#     if (z[i] <= pivot_point){
#       rho_symmetric = 2*ranking_interp(pivot_point) - ranking_interp(z[i])
#     }
#     else {
#       rho_symmetric = ranking_interp(2*pivot_point - z[i])
#     }
#     rho_normal[i] = min(rho_symmetric, rho_greater[i])
#   }
#   
#   rho_active = rho_greater - rho_normal
#   f_active = rho_active/rho_greater
#   z_active = z[min(which(f_active >= 0.5))]
#   
#   minimo = sigma_minima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
#     quantile(na.omit(lum), probs = 0.50)
#   maximo = sigma_maxima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
#     quantile(na.omit(lum), probs = 0.50)
#   
#   library(ggplot2)
#   
#   df = data.frame(z, log10(rho_normal))
#   df2 = data.frame(z, log10(rho_greater))
#   
#   print(ggplot(df, aes(x = df$z)) +
#           ylim(0, max(na.omit(log10(rho_greater)))) +
#           geom_line(aes(y = df$log10.rho_normal.), color = colour, linetype = 'dashed') +
#           geom_line(aes(y = df2$log10.rho_greater.), color = 'black') +
#           geom_rect(aes(xmin = minimo, xmax = maximo, ymin = -Inf, ymax = Inf), color = 'grey', alpha = 0.005) +
#           geom_vline(xintercept = z_active) +
#           theme_bw() + 
#           xlab(xlabel) + 
#           ylab('log (rho > x)'))
# 
#   return(z_active)
# }



activeness = function(lum, xlabel, colour, V_max_inverse){
  
  indice_ordenado = sort.int(lum, na.last = TRUE, index.return = TRUE)$ix[1:length(na.omit(lum))]
  
  x_ordenado = lum[indice_ordenado]
  #ranking = na.omit(sort.int(V_max_inverse, na.last = TRUE, index.return = TRUE)$x)
  ranking = V_max_inverse[indice_ordenado]
  ranking = cumsum(ranking)
  ranking = ranking/ranking[length(ranking)]
  ranking_interp = approxfun(x_ordenado, ranking, yleft = 0, yright = 1)
  
  ranking_interp2 = approxfun(ranking, x_ordenado)
  total_dens = max(ranking)
  mid_dens = total_dens/2
  pivot_point = ranking_interp2(mid_dens)
  #pivot_point = 0
  
  
  rho = total_dens
  z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
  rho_greater = rho - ranking_interp(z)
  rho_normal = rep(0, length(rho_greater))
  
  
  for (i in 1:length(z)){
    if (z[i] <= pivot_point){
      rho_symmetric = 2*ranking_interp(pivot_point) - ranking_interp(z[i])
    }
    else {
      rho_symmetric = ranking_interp(2*pivot_point - z[i])
    }
    rho_normal[i] = min(rho_symmetric, rho_greater[i])
  }
  
  rho_active = rho_greater - rho_normal
  f_active = rho_active/rho_greater
  z_active = z[min(which(f_active >= 0.5))]
  
  plot(z, rho_greater, pch='.', log = 'y', ylim = c(0.001, 1), type='n')
  lines(z, rho_greater, ylim = c(0.001, 1))
  lines(z, rho_normal, col='blue')
  lines(z, rho_active, col='red')
  abline(v=z_active)
  return(z_active)
}




activeness(L_excess_IR2_copy, 'Luminosity excess (IR W2 band)', 'red', V_max_IR2_inverse_copy)
activeness(L_excess_Radio_copy, 'Luminosity excess (Radio)', 'green', V_max_Radio_inverse_copy)
activeness(L_excess_Hard, 'Luminosity excess (Hard X-Rays)', 'blue', V_max_Hard_inverse)

activeness(color_W2_W1_excess_copy, 'IR color excess', 'red', V_max_color_IR_inverse_copy)
activeness(color_Radio_IR_excess_copy, 'Radio-IR color excess', 'green', V_max_color_Radio_inverse_copy)
activeness(color_X_excess_copy, 'X-Rays color excess', 'blue', V_max_color_X_inverse_copy)


############## CÁLCULO DE LOS THRESHOLDS PARA CADA LUMINOSIDAD Y COLOR ################

threshold_colorIR_lit = 0.4*0.8 + W1_W2_norm # Stern threshold ---> W1-W2 > 0.8
threshold_colorIR_lit2 = 0.4*0.5 + W1_W2_norm # Assef threshold ---> W1-W2 > 0.5
threshold_color_Radio_IR_lit = 0.23 + nuradio_nu24 # Ibar threshold ---> q_24 < -0.23
threshold_Radio_lit = 40 # L > 10^40 erg/s ---> Radio-loud AGN
threshold_Hard_lit = 42 # L > 10^42 erg/s ---> AGN
threshold_colorX_lit = log10(7)

threshold_ir2_excess = as.numeric(activeness(L_excess_IR2, 'Luminosity excess (IR W2 band)', 'red', V_max_IR2_inverse_copy))
threshold_colorIR_excess = as.numeric(activeness(color_W2_W1_excess, 'IR color excess', 'red', V_max_color_IR_inverse_copy))
threshold_radio_excess = as.numeric(activeness(L_excess_Radio, 'Luminosity excess (Radio)', 'green', V_max_Radio_inverse_copy))
threshold_color_Radio_IR_excess = as.numeric(activeness(color_Radio_IR_excess, 'Radio-IR color excess', 'green', V_max_color_Radio_inverse_copy))
threshold_hard_excess = as.numeric(activeness(L_excess_Hard, 'Luminosity excess (Hard X-Rays)', 'blue', V_max_Hard_inverse_copy))
threshold_colorX_excess = as.numeric(activeness(color_X_excess, 'X-Rays color excess', 'blue', V_max_color_X_inverse_copy))


################################### FINAL PLOTS ########################################

# Función que dados los excesos de color y luminosidad de una galaxia le asigna
# un color de acuerdo con nuestra clasificación:
# Rojo: active galaxy. Verde: luminosity extreme. Amarillo: color extreme. Azul: normal galaxy
# get_galaxy_type hace esto para un objeto dados sus excesos de luminosidad y color

color_active = alpha('red4', 0.2)


get_galaxy_type = function(L_i, color_i, threshold_L, threshold_color, alpha_value){
  galaxy_type = alpha("blue", alpha_value)
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

galaxy_class_IR = galaxy_classification(L_excess_IR2, color_W2_W1_excess,
                                        threshold_ir2_excess, threshold_colorIR_excess, 0.2)

plot(L_excess_IR2, color_W2_W1_excess, col = galaxy_class_IR, pch = 16, ylim = c(-0.5, 1),
     xlab = "IR (W2 band) luminosity excess", ylab = "IR color excess")
abline(h = threshold_colorIR_excess, lty = 2, lwd = 2)
abline(v = threshold_ir2_excess, lty = 4, lwd = 2)

data_IR = which(is.na(L_excess_IR2 & color_W2_W1_excess)==FALSE)
#plot(M[data_IR], SFR[data_IR], col = galaxy_class_IR[data_IR], pch = 16, xlim = c(8.5, 11.5), ylim = c(-2, 1),
#     xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")


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



galaxy_class_Radio = galaxy_classification(L_excess_Radio, color_Radio_IR_excess,
                                           threshold_radio_excess, threshold_color_Radio_IR_excess, 0.3)

plot(L_excess_Radio, color_Radio_IR_excess, col = galaxy_class_Radio, pch = 16, xlim = c(-1.5, 3), ylim = c(-2, 3),
     xlab = "Radio luminosity excess", ylab = "Radio-IR color excess")
abline(h = threshold_color_Radio_IR_excess, lty = 2, lwd = 2)
abline(v = threshold_radio_excess, lty = 4, lwd = 2)

data_Radio = which(is.na(L_excess_Radio & color_Radio_IR_excess)==FALSE)
#plot(M[data_Radio], SFR[data_Radio], col = galaxy_class_Radio[data_Radio], pch = 16, xlim = c(8.5, 11.5), ylim = c(-2, 1),
#     xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")

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

##################### X-RAYS DON'T HAVE THRESHOLD IN COLOR SO WE DEFINE THEIR OWN FUNCTIONS ####################


get_galaxy_type_X = function(L_i, threshold_L, alpha_value){
  
  galaxy_type = alpha("blue", alpha_value)
  if ((is.na(L_i) == FALSE) & (L_i > threshold_L)){
    galaxy_type = alpha("green", alpha_value)
  }

  return (galaxy_type)
}
# galaxy_classification lo generaliza y recorre todos los objetos para cada banda

galaxy_classification_X = function(L, threshold_L, alpha_value){
  N = length(L)
  galaxy_class = as.vector(rep(alpha("black", 1), N))
  for (i in 1:N){
    galaxy_class[i] = get_galaxy_type_X(L[i], threshold_L, alpha_value)
  }
  return(galaxy_class)
}


galaxy_class_X = galaxy_classification_X(L_excess_Hard, threshold_hard_excess, 1)

plot(L_excess_Hard, color_X_excess, col = galaxy_class_X, pch = 16, ylim = c(-1, 2),
     xlab = "Hard X-Rays luminosity excess", ylab = "X-Rays color excess")
#abline(h = threshold_colorX_excess, lty = 2, lwd = 2)
abline(v = threshold_hard_excess, lty = 4, lwd = 2)

#data_X = which(is.na(L_excess_Hard & color_X_excess)==FALSE)
#plot(M[data_X], SFR[data_X], col = galaxy_class_X[data_X], pch = 16, xlim = c(8.5, 11.5), ylim = c(-2, 1),
#     xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")


plot(L_Hard, color_X, col = galaxy_class_X, pch = 16,
     xlab = "log[L (Hard X-Rays)]", ylab = "X-Rays color", ylim = c(-1.5, 2.5))
#abline(h = threshold_colorX_lit, lty = 2, lwd = 2)
abline(v = 42, lty = 4, lwd = 2)
axis(2, col = "black", las = 0) 

#######################################################################################################################


# Now we'll plot our different galaxy classes separately in the M*-SFR diagram:

#IR

data_IR_active = which((galaxy_class_IR == alpha('red4', 0.2)))

plot(M[data_IR_active], SFR[data_IR_active], col = galaxy_class_IR[data_IR_active], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_W2_W1, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .1, add = TRUE, lty = 1)


data_IR_lum = which((galaxy_class_IR == alpha('green', 0.2)))
plot(M[data_IR_lum], SFR[data_IR_lum], col = galaxy_class_IR[data_IR_lum], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_W2_W1, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .1, add = TRUE, lty = 1)


data_IR_col = which((galaxy_class_IR == alpha('gold', 0.2)))
plot(M[data_IR_col], SFR[data_IR_col], col = galaxy_class_IR[data_IR_col], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_W2_W1, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .1, add = TRUE, lty = 1)


data_IR_normal = which((galaxy_class_IR == alpha('blue', 0.2)))
plot(M[data_IR_normal], SFR[data_IR_normal], col = galaxy_class_IR[data_IR_normal], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_W2_W1, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_W2_W1, levels = .1, add = TRUE, lty = 1)

#RADIO
data_Radio_nondet = which((is.na(L_excess_Radio)==FALSE) & (is.na(color_Radio_IR_excess)==TRUE))

data_Radio_active = which((is.na(L_excess_Radio)==FALSE) & (is.na(color_Radio_IR_excess)==FALSE) &
                         (galaxy_class_Radio == alpha('red4', 0.3)))
plot(M[data_Radio_active], SFR[data_Radio_active], col = galaxy_class_Radio[data_Radio_active], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .1, add = TRUE, lty = 1)


data_Radio_lum = which((is.na(L_excess_Radio)==FALSE) & (is.na(color_Radio_IR_excess)==FALSE) &
                         (galaxy_class_Radio == alpha('green', 0.3)))
plot(M[data_Radio_lum], SFR[data_Radio_lum], col = galaxy_class_Radio[data_Radio_lum], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .1, add = TRUE, lty = 1)


data_Radio_col = which((is.na(L_excess_Radio)==FALSE) & (is.na(color_Radio_IR_excess)==FALSE) &
                            (galaxy_class_Radio == alpha('gold', 0.3)))
plot(M[data_Radio_col], SFR[data_Radio_col], col = galaxy_class_Radio[data_Radio_col], pch = 16, 
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .1, add = TRUE, lty = 1)


data_Radio_normal = which((is.na(L_excess_Radio)==FALSE) & (is.na(color_Radio_IR_excess)==FALSE) &
                         (galaxy_class_Radio == alpha('blue', 0.03)))
plot(M[data_Radio_normal], SFR[data_Radio_normal], col = galaxy_class_Radio[data_Radio_normal], pch = 16,
     xlim = c(8.5, 11.5), ylim = c(-2, 1), xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_Radio_IR, levels = .1, add = TRUE, lty = 1)


#X-RAYS
data_X_nondet = which((is.na(L_excess_Hard)==FALSE) & (is.na(color_X_excess)==TRUE))

data_X_lum = which((is.na(L_excess_Hard)==FALSE) & (galaxy_class_X == alpha('green', 1)))
plot(M[data_X_lum], SFR[data_X_lum], col = galaxy_class_X[data_X_lum], pch = 16, xlim = c(8.5, 11.5), ylim = c(-2, 1),
     xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_X, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_X, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_X, levels = .1, add = TRUE, lty = 1)


data_X_normal = which((is.na(L_excess_Hard)==FALSE) & (galaxy_class_X == alpha('blue', 1)))
plot(M[data_X_normal], SFR[data_X_normal], col = galaxy_class_X[data_X_normal], pch = 16, xlim = c(8.5, 11.5), ylim = c(-2, 1),
     xlab = "log[M* (solar masses)]", ylab = "log[SFR (solar masses/yr)")
contour(M_X, SFR_y, contours_color_X, levels = .9, add = TRUE, lty = 3)
contour(M_X, SFR_y, contours_color_X, levels = .5, add = TRUE, lty = 2)
contour(M_X, SFR_y, contours_color_X, levels = .1, add = TRUE, lty = 1)

##########################################################################################################


sets = function(v, set){
  v_new = v
  v_new[set] = NA
  return(v_new)
}


# IR color

color_W2_W1_excess_grey = color_W2_W1_excess
color_W2_W1_excess_red = sets(color_W2_W1_excess_grey, data_IR_nondet)
color_W2_W1_excess_yellow = sets(color_W2_W1_excess_red, data_IR_active)
color_W2_W1_excess_green = sets(color_W2_W1_excess_yellow, data_IR_col)
color_W2_W1_excess_blue = sets(color_W2_W1_excess_green, data_IR_lum)

V_max_color_IR_inverse_grey = V_max_color_IR_inverse
V_max_color_IR_inverse_red = sets(V_max_color_IR_inverse_grey, data_IR_nondet)
V_max_color_IR_inverse_yellow = sets(V_max_color_IR_inverse_red, data_IR_active)
V_max_color_IR_inverse_green = sets(V_max_color_IR_inverse_yellow, data_IR_col)
V_max_color_IR_inverse_blue = sets(V_max_color_IR_inverse_green, data_IR_lum)

color_W2_W1_excess_grey = na.omit(color_W2_W1_excess_grey)
color_W2_W1_excess_red = na.omit(color_W2_W1_excess_red)
color_W2_W1_excess_yellow = na.omit(color_W2_W1_excess_yellow)
color_W2_W1_excess_green = na.omit(color_W2_W1_excess_green)
color_W2_W1_excess_blue = na.omit(color_W2_W1_excess_blue)

V_max_color_IR_inverse_grey = na.omit(V_max_color_IR_inverse_grey)
V_max_color_IR_inverse_red = na.omit(V_max_color_IR_inverse_red)
V_max_color_IR_inverse_yellow = na.omit(V_max_color_IR_inverse_yellow)
V_max_color_IR_inverse_green = na.omit(V_max_color_IR_inverse_green)
V_max_color_IR_inverse_blue = na.omit(V_max_color_IR_inverse_blue)

range_color_IR = seq(min(na.omit(color_W2_W1_excess)), max(na.omit(color_W2_W1_excess)), 
                     (max(na.omit(color_W2_W1_excess)) - min(na.omit(color_W2_W1_excess)))/50)


weighted.hist(color_W2_W1_excess_grey, V_max_color_IR_inverse_grey, col = alpha('grey', 1), breaks=seq(-0.25,0.25,0.5/50), xaxis = TRUE)
weighted.hist(color_W2_W1_excess_red, V_max_color_IR_inverse_red, col = alpha('red', 1), breaks=seq(-0.25,0.25,0.5/50), add=T, xaxis = FALSE)
weighted.hist(color_W2_W1_excess_yellow, V_max_color_IR_inverse_yellow, breaks=seq(-0.25,0.25,0.5/50),col = alpha('yellow',1), add=T, xaxis = FALSE)
weighted.hist(color_W2_W1_excess_green, V_max_color_IR_inverse_green, breaks=seq(-0.25,0.25,0.5/50), col = alpha('green', 1), add=T, xaxis = FALSE)
weighted.hist(color_W2_W1_excess_blue, V_max_color_IR_inverse_blue, breaks=seq(-0.25,0.25,0.5/50), col = alpha('blue', 1),add=T, xaxis = FALSE)


# L_IR2 luminosity

L_excess_IR2_grey = L_excess_IR2
L_excess_IR2_red = sets(L_excess_IR2_grey, data_IR_nondet)
L_excess_IR2_yellow = sets(L_excess_IR2_red, data_IR_active)
L_excess_IR2_green = sets(L_excess_IR2_yellow, data_IR_col)
L_excess_IR2_blue = sets(L_excess_IR2_green, data_IR_lum)

V_max_IR2_inverse_grey = V_max_IR2_inverse
V_max_IR2_inverse_red = sets(V_max_IR2_inverse_grey, data_IR_nondet)
V_max_IR2_inverse_yellow = sets(V_max_IR2_inverse_red, data_IR_active)
V_max_IR2_inverse_green = sets(V_max_IR2_inverse_yellow, data_IR_col)
V_max_IR2_inverse_blue = sets(V_max_IR2_inverse_green, data_IR_lum)

L_excess_IR2_grey = na.omit(L_excess_IR2_grey)
L_excess_IR2_red = na.omit(L_excess_IR2_red)
L_excess_IR2_yellow = na.omit(L_excess_IR2_yellow)
L_excess_IR2_green = na.omit(L_excess_IR2_green)
L_excess_IR2_blue = na.omit(L_excess_IR2_blue)

V_max_IR2_inverse_grey = na.omit(V_max_IR2_inverse_grey)
V_max_IR2_inverse_red = na.omit(V_max_IR2_inverse_red)
V_max_IR2_inverse_yellow = na.omit(V_max_IR2_inverse_yellow)
V_max_IR2_inverse_green = na.omit(V_max_IR2_inverse_green)
V_max_IR2_inverse_blue = na.omit(V_max_IR2_inverse_blue)

range_IR2 = seq(min(na.omit(L_excess_IR2)), max(na.omit(L_excess_IR2)), 
                     (max(na.omit(L_excess_IR2)) - min(na.omit(L_excess_IR2)))/50)

weighted.hist(L_excess_IR2_grey, V_max_IR2_inverse_grey, col = alpha('grey', 1), breaks = seq(-1,1, 2/50), xaxis = TRUE)
weighted.hist(L_excess_IR2_red, V_max_IR2_inverse_red, col = alpha('red', 1), breaks = seq(-1,1, 2/50), add=T, xaxis = FALSE)
weighted.hist(L_excess_IR2_yellow, V_max_IR2_inverse_yellow, col = alpha('yellow', 1), breaks = seq(-1,1, 2/50), add=T, xaxis = FALSE)
weighted.hist(L_excess_IR2_green, V_max_IR2_inverse_green, col = alpha('green', 1), breaks = seq(-1,1, 2/50), add=T, xaxis = FALSE)
weighted.hist(L_excess_IR2_blue, V_max_IR2_inverse_blue, col = alpha('blue', 1), breaks = seq(-1,1, 2/50), add=T, xaxis = FALSE)


# color Radio-IR

color_Radio_IR_excess_grey = color_Radio_IR_excess
color_Radio_IR_excess_red = sets(color_Radio_IR_excess_grey, data_Radio_nondet)
color_Radio_IR_excess_yellow = sets(color_Radio_IR_excess_red, data_Radio_active)
color_Radio_IR_excess_green = sets(color_Radio_IR_excess_yellow, data_Radio_col)
color_Radio_IR_excess_blue = sets(color_Radio_IR_excess_green, data_Radio_lum)

V_max_color_Radio_IR_inverse_grey = V_max_color_Radio_inverse
V_max_color_Radio_IR_inverse_red = sets(V_max_color_Radio_IR_inverse_grey, data_Radio_nondet)
V_max_color_Radio_IR_inverse_yellow = sets(V_max_color_Radio_IR_inverse_red, data_Radio_active)
V_max_color_Radio_IR_inverse_green = sets(V_max_color_Radio_IR_inverse_yellow, data_Radio_col)
V_max_color_Radio_IR_inverse_blue = sets(V_max_color_Radio_IR_inverse_green, data_Radio_lum)

color_Radio_IR_excess_grey = na.omit(color_Radio_IR_excess_grey)
color_Radio_IR_excess_red = na.omit(color_Radio_IR_excess_red)
color_Radio_IR_excess_yellow = na.omit(color_Radio_IR_excess_yellow)
color_Radio_IR_excess_green = na.omit(color_Radio_IR_excess_green)
color_Radio_IR_excess_blue = na.omit(color_Radio_IR_excess_blue)

V_max_color_Radio_IR_inverse_grey = na.omit(V_max_color_Radio_IR_inverse_grey)
V_max_color_Radio_IR_inverse_red = na.omit(V_max_color_Radio_IR_inverse_red)
V_max_color_Radio_IR_inverse_yellow = na.omit(V_max_color_Radio_IR_inverse_yellow)
V_max_color_Radio_IR_inverse_green = na.omit(V_max_color_Radio_IR_inverse_green)
V_max_color_Radio_IR_inverse_blue = na.omit(V_max_color_Radio_IR_inverse_blue)

range_color_Radio = seq(min(na.omit(color_Radio_IR_excess)), max(na.omit(color_Radio_IR_excess)), 
                     (max(na.omit(color_Radio_IR_excess)) - min(na.omit(color_Radio_IR_excess)))/20)

weighted.hist(color_Radio_IR_excess_grey, V_max_color_Radio_IR_inverse_grey, col = alpha('grey', 1), breaks = seq(-2,2.5,4.5/20), xaxis = TRUE)
weighted.hist(color_Radio_IR_excess_red, V_max_color_Radio_IR_inverse_red, col = alpha('red', 1), breaks = seq(-2,2.5,4.5/20),  xaxis = FALSE)
weighted.hist(color_Radio_IR_excess_yellow, V_max_color_Radio_IR_inverse_yellow, breaks = seq(-2,2.5,4.5/20),col = alpha('yellow',1), add=T, xaxis = FALSE)
weighted.hist(color_Radio_IR_excess_green, V_max_color_Radio_IR_inverse_green, breaks = seq(-2,2.5,4.5/20), col = alpha('green', 1), add=T, xaxis = FALSE)
weighted.hist(color_Radio_IR_excess_blue, V_max_color_Radio_IR_inverse_blue, breaks = seq(-2,2.5,4.5/20), col = alpha('blue', 1), add=T, xaxis = FALSE)


# Radio luminosity:

L_excess_Radio_grey = L_excess_Radio
L_excess_Radio_red = sets(L_excess_Radio_grey, data_Radio_nondet)
L_excess_Radio_yellow = sets(L_excess_Radio_red, data_Radio_active)
L_excess_Radio_green = sets(L_excess_Radio_yellow, data_Radio_col)
L_excess_Radio_blue = sets(L_excess_Radio_green, data_Radio_lum)

V_max_Radio_inverse_grey = V_max_Radio_inverse
V_max_Radio_inverse_red = sets(V_max_Radio_inverse_grey, data_Radio_nondet)
V_max_Radio_inverse_yellow = sets(V_max_Radio_inverse_red, data_Radio_active)
V_max_Radio_inverse_green = sets(V_max_Radio_inverse_yellow, data_Radio_col)
V_max_Radio_inverse_blue = sets(V_max_Radio_inverse_green, data_Radio_lum)

L_excess_Radio_grey = na.omit(L_excess_Radio_grey)
L_excess_Radio_red = na.omit(L_excess_Radio_red)
L_excess_Radio_yellow = na.omit(L_excess_Radio_yellow)
L_excess_Radio_green = na.omit(L_excess_Radio_green)
L_excess_Radio_blue = na.omit(L_excess_Radio_blue)

V_max_Radio_inverse_grey = na.omit(V_max_Radio_inverse_grey)
V_max_Radio_inverse_red = na.omit(V_max_Radio_inverse_red)
V_max_Radio_inverse_yellow = na.omit(V_max_Radio_inverse_yellow)
V_max_Radio_inverse_green = na.omit(V_max_Radio_inverse_green)
V_max_Radio_inverse_blue = na.omit(V_max_Radio_inverse_blue)

range_Radio = seq(min(na.omit(L_excess_Radio)), max(na.omit(L_excess_Radio)), 
                (max(na.omit(L_excess_Radio)) - min(na.omit(L_excess_Radio)))/50)

weighted.hist(L_excess_Radio_grey, V_max_Radio_inverse_grey, col = alpha('grey', 1), breaks = seq(-1,2,3/20), xaxis = TRUE)
weighted.hist(L_excess_Radio_red, V_max_Radio_inverse_red, col = alpha('red', 1), breaks = seq(-1,2,3/20), add=T, xaxis = FALSE)
weighted.hist(L_excess_Radio_yellow, V_max_Radio_inverse_yellow, col = alpha('yellow', 1), breaks = seq(-1,2,3/20), add=T, xaxis = FALSE)
weighted.hist(L_excess_Radio_green, V_max_Radio_inverse_green, col = alpha('green', 1), breaks = seq(-1,2,3/20), add=T, xaxis = FALSE)
weighted.hist(L_excess_Radio_blue, V_max_Radio_inverse_blue, col = alpha('blue', 1), breaks = seq(-1,2,3/20), add=T, xaxis = FALSE)


# X-Rays color

color_X_excess_grey = color_X_excess
color_X_excess_green = sets(color_X_excess_grey, data_X_nondet)
color_X_excess_blue = sets(color_X_excess_green, data_X_lum)

V_max_color_X_inverse_grey = V_max_color_X_inverse
V_max_color_X_inverse_green = sets(V_max_color_X_inverse_grey, data_X_nondet)
V_max_color_X_inverse_blue = sets(V_max_color_X_inverse_green, data_X_lum)

color_X_excess_grey = na.omit(color_X_excess_grey)
color_X_excess_green = na.omit(color_X_excess_green)
color_X_excess_blue = na.omit(color_X_excess_blue)

V_max_color_X_inverse_grey = na.omit(V_max_color_X_inverse_grey)
V_max_color_X_inverse_green = na.omit(V_max_color_X_inverse_green)
V_max_color_X_inverse_blue = na.omit(V_max_color_X_inverse_blue)

range_color_X = seq(min(na.omit(color_X_excess)), max(na.omit(color_X_excess)), 
                     (max(na.omit(color_X_excess)) - min(na.omit(color_X_excess)))/50)

weighted.hist(color_X_excess_grey, V_max_color_X_inverse_grey, col = alpha('grey', 1), breaks=seq(-1.5,2,3.5/25), xaxis = TRUE)
weighted.hist(color_X_excess_green, V_max_color_X_inverse_green, col = alpha('green', 1), breaks=seq(-1.5,2,3.5/25), add=T, xaxis = FALSE)
weighted.hist(color_X_excess_blue, V_max_color_X_inverse_blue, breaks=seq(-1.5,2,3.5/25),col = alpha('blue', 1), add=T, xaxis = FALSE)


# Hard luminosity
L_excess_Hard_grey = L_excess_Hard
L_excess_Hard_green = sets(L_excess_Hard_grey, data_X_nondet)
L_excess_Hard_blue = sets(L_excess_Hard_green, data_X_lum)

V_max_Hard_inverse_grey = V_max_Hard_inverse
V_max_Hard_inverse_green = sets(V_max_Hard_inverse_grey, data_X_nondet)
V_max_Hard_inverse_blue = sets(V_max_Hard_inverse_green, data_X_lum)

L_excess_Hard_grey = na.omit(L_excess_Hard_grey)
L_excess_Hard_green = na.omit(L_excess_Hard_green)
L_excess_Hard_blue = na.omit(L_excess_Hard_blue)

V_max_Hard_inverse_grey = na.omit(V_max_Hard_inverse_grey)
V_max_Hard_inverse_green = na.omit(V_max_Hard_inverse_green)
V_max_Hard_inverse_blue = na.omit(V_max_Hard_inverse_blue)

range_Hard = seq(min(na.omit(L_excess_Hard)), max(na.omit(L_excess_Hard)), 
                (max(na.omit(L_excess_Hard)) - min(na.omit(L_excess_Hard)))/50)

weighted.hist(L_excess_Hard_grey, V_max_Hard_inverse_grey, col = alpha('grey', 1), breaks = seq(-2,2,4/10), xaxis = TRUE)
weighted.hist(L_excess_Hard_green, V_max_Hard_inverse_green, col = alpha('green', 1), breaks = seq(-2,2,4/10), add=T, xaxis = FALSE)
weighted.hist(L_excess_Hard_blue, V_max_Hard_inverse_blue, col = alpha('blue', 1), breaks = seq(-2,2, 4/10), add=T, xaxis = FALSE)



