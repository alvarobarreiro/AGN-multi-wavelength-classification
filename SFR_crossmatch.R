
#############################    STATISTICAL STUFF    #################################

setwd('D:/AGN') #directorio de trabajo
luminosities <- read.csv('Luminosities.csv', header=TRUE) #léeme el archivo con las luminosidades 
#de "califa_multiwavelength"
SFR <- read.csv('SFR.csv', header=FALSE) #léeme el archivo con las SFR de Catalán-Torrecilla
names(SFR) <- c('ID', 'Name', 'Dist', 'FHaobs', 'LHaobs', 'SFRHacorr', 'e_SFRHacorr',
                'L22um', 'e_L22um', 'RA', 'DEC', 'recno') #pon a SFR una cabecera con los nombres de las columnas

a <- as.character(luminosities$Name) #Vector con los nombres de las galaxias de "luminosities"
#(son exactamente las mismas que en "califa_multiwavelength")
b <- as.character(gsub(" ", "", SFR$Name, fixed = TRUE)) #Vector con los nombres de las galaxias en "SFR"
#Esta última línea de código simplemente elimina los espacios en los nombres para "SFR". De lo contrario no
#podemos comparar en el bucle de abajo entre los nombres de los objetos en "luminosities" (que no tienen espacios)
#y en "SFR" (que sí los tienen) Ejemplo: IC 5376 en SFR sería IC5376 en CALIFA. Se trata de quitar el espacio.

len1 <- length(a)
len2 <- length(b)
SFR1 <- vector() #inicializamos en vector vacío una variable SFR1 que nos dará la tasa de formación
#estelar para cada objeto de "califa_multiwavelength" (si es que está disponible en "SFR")

#Bucle:
for (i in 1:len1){ #para cada galaxia de "califa_multiwavelength" (1062)
  for (j in 1:len2){ #recórreme todas las de "SFR" (271)
    if (a[i]==b[j]){ #si el objeto de CALIFA está en SFR (mismo nombre)
      SFR1[i] <- SFR$SFRHacorr[j] #guárdame en el vector SFR1 el valor de su tasa de formación estelar
    }
  }
}#Por alguna puñetera razón, SFR1 sólo tiene 1039 entradas en lugar de las 1062 que debería tener.
#Inspeccionando los datos a mano se ve que las galaxias de la 1040 a la 1062 en "luminosities" no están
#en "SFR", y por tanto podemos asignarle un 'NA' (not available):
for (k in 1:23){ 
  SFR1<-c(SFR1, NA)
} 
#Ahora SFR1 sí tiene 1062 entradas con información de las tasas de formación estelar
#269 objetos de los 271 de Catalán-Torrecilla estaban en nuestro califa_multiwavelength
#Por tanto SFR1 contiene 269 valores de la tasa de formación estelar. El resto de 
#las 1069 entradas contiene un 'NA' (not available).



##############################     INFRARRED 22 MICRONS     ###################################

lumir<-luminosities$IR.luminosity_221.erg.#vector con las luminosidades en el IR a 22.1 micras
lumir<-lumir*(3*10^8/(22.1*10^-6))#pasamos de erg/s/Hz a erg/s
lumir<-lumir/(10^43) #Así es como lo hace Catalán-Torrecilla...
sfr<-SFR1 #renombramos el vector SFR1 
loglumir<-log10(lumir) #tomamos logaritmos
logsfr<-log10(sfr) #tomamos logaritmos
plot(loglumir, logsfr, xlim = c(-3, 2), pch = 16, col = "orange", main = "IR (22.1 microns)",
     xlab = "log [L(22 microns)/10^43 (erg/s)]", ylab = "log [SFR (solar mass/yr)]") 
#representamos sfr frente a lumir (logaritmos)

ajuste1<- lm(logsfr ~ loglumir) #hacemos el ajuste lineal por mínimos cuadrados
ord<-as.numeric(ajuste1$coefficients[1]) #sacamos la ordenada en el origen
pend<-as.numeric(ajuste1$coefficients[2]) #y la pendiente
abline(a=0.377, b=0.564, lw=2, col="black") #Aquí iría una línea con los valores de los autores
abline(a=ord, b=pend, lw=2, col="blue") #representamos la recta del ajuste lineal
legend("topleft", legend = c("Catalán-Torrecilla", "Our fit"), fill = c("black", "blue"))



#Cálculo de la S:

S_IR<-as.numeric() #inicializamos en vector vacío nuestro parámetro clave "S" para el infrarrojo
for (i in 1:len1){ #para cada galaxia de CALIFA
  S_IR[i]<-(ord+pend*loglumir[i]-logsfr[i]) #calcúlame el valor de "S" para el infrarrojo
}
S_IR_raw <- S_IR #nos será útil luego
S_IR<-na.omit(S_IR) #na.omit es para quitar NaN, NA y demás purria (si no, no nos deja calcular nada)



#Ahora tomamos la variable S normal, luego centramos la variable S en la moda, y finalmente
#dividimos la variable centrada entre una desviación típica (lo que viene siendo tipificar
#una variable en estadística, vamos, sólo que en este caso estamos utilizando la moda en lugar
#de la media como parámetro de centralización: S(normalizada) = (S - moda)/sigma).
#Histogramamos cada paso:

histograma <-hist(S_IR, breaks = seq(-0.75, 1.35, 0.1), freq=FALSE, col = "orange", main = "Histogram of S (IR)",
                  xlab = "S (IR)")#histograma de la variable S
moda_IR <- histograma$mids[which.max(histograma$counts)] #cálculo de la moda usando el histograma



# Normalizaciones (estimaciones de sigma, calculamos 3 distintas):

sigma_moda_IR <- sqrt(sum((S_IR - moda_IR)^2)/(length(S_IR) - 1))
#definición tradicional de la desv. típica pero usando la moda en lugar de la media, como se dijo arriba

sigma_percentil_IR <- as.numeric((quantile(S_IR, probs = 0.84) - quantile(S_IR, probs = 0.16))/2) 
#dev. típica estimada mediante el rango intercuartílico (prescripción de los percentiles)

S_IR_izq <- as.numeric()
for (i in 1:length(S_IR)){
  if (S_IR[i] <= moda_IR)
    S_IR_izq <- c(S_IR_izq, S_IR[i])
}
sigma_izquierda_IR <- sqrt(sum((S_IR_izq - moda_IR)^2)/(length(S_IR_izq) - 1))
#dev. típica estimada usando sólo la parte de los datos a la izquierda de la moda





#Histogramamos ahora la VARIABLE CENTRADA y le plantamos la correspondiente gaussiana por encima

S_IR_centrada <- S_IR - moda_IR #centramos S y su histograma en la moda
histograma_cent <-hist(S_IR_centrada, breaks = seq(-0.75, 1.35, 0.1), freq=FALSE, col = "orange",
                       main = "Histogram of S (IR, mode-centered)", xlab = "S (IR, mode-centered)")
xdens <- seq(-1.5, 1.5, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = moda_IR, sd = sigma_izquierda_IR), col = "black", lwd = 2)#gaussiana centrada N(0,s)
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 0.290"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))





#Histogramamos ahora la VARIABLE TIPIFICADA (o normalizada) y le plantamos la correspondiente gaussiana por encima

S_IR_norm <- S_IR_centrada/sigma_izquierda_IR
histograma_norm <- hist(S_IR_norm, 
                        breaks = seq(-0.75/sigma_izquierda_IR, 1.35/sigma_izquierda_IR, 0.10/sigma_izquierda_IR), 
                        main = "Histogram of S (IR, normalized)", freq=FALSE, 
                        col = "orange",xlab = "S (IR, normalized)")

#Lo de los breaks aquí parece una locura, pero es para que el bineado coincida exactamente
#con el del histograma anterior de la variable centrada "histograma_cent". Si no hacemos esto,
#el bienado cambia y la moda deja de estar en el 0 para la normalización

xdens <- seq(-5, 5, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = 0, sd = 1), col = "black", lwd = 2)#gaussiana tipificada (N(0,1))
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 1"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))




# ----------------------------------------------------------------------------------
#Construimos ahora un histograma simétrico a la parte SF de la izquierda de la moda
#Es decir, el siguiente histograma será sólo la parte SF (moda y su izquierda)
#y su reflejo en la parte a la derecha de la moda
# ----------------------------------------------------------------------------------

histograma_sim<-histograma_norm
longitud<-length(histograma_sim$mids)
counter <- 0
for (i in 1:longitud){
  if (histograma_sim$mids[i]<0){
    counter<-counter+1
  }
}

for (j in (counter+2):(2*counter+1)) {
  histograma_sim$counts[j]<-histograma_sim$counts[2*counter - j + 2] 
  histograma_sim$density[j]<-histograma_sim$density[2*counter - j + 2]
}

for (j in (2*counter+2):(longitud)) {
  histograma_sim$counts[j]<-0
  histograma_sim$density[j]<-0
}

plot(histograma_sim, freq = FALSE, col = "orange", main = "Histogram of S (IR, SF + symmetric)",
     xlab = "S (IR, symmetric)")




# ---------------------------------------------------------------------------------------------------------
#Ahora podemos jugar a poner juntos ambos histogramas: el histograma total, con SF + AGN, y éste simétrico
#último que sería una especie de "sólo SF". Así tendremos una idea del exceso correspondiente a AGN.
# ---------------------------------------------------------------------------------------------------------

library(scales)
plot(histograma_norm, freq = FALSE, col='skyblue', border=F,
     main = "Histogram of S (IR)", xlab = "S (IR)")
plot(histograma_sim, freq = FALSE, add=T, col=scales::alpha('red', .5), border=F)
legend("topright", legend = c("SF + symmetric", "SF + AGN", "Overlap"), 
       fill = c(scales::alpha('red', .35), "skyblue", scales::alpha('red', .5)))

kernel<-density(S_IR_norm)
plot(kernel, main = "Kernel density of S (IR, normalized)")
polygon(kernel, col = "red", border = "black")





# ---------------------------------------------------------------------------------------------------------
# El siguiente paso va a ser histogramar exactamente lo mismo pero restando a la parte positiva del histograma
# las cuentas de la parte simétrica negativa. Así tendremos una idea del exceso correspondiente a AGN.
# ---------------------------------------------------------------------------------------------------------

histograma_diff<-histograma_norm
longitud<-length(histograma_diff$mids) #longitud es el número de intervalos del histograma
counter<-0 #counter van a ser el número de breaks que tenemos a la izquierda del cero
for (i in 1:longitud){
  if (histograma_diff$mids[i]<0){
    counter<-counter+1
  }
}

for (j in (counter+2):(2*counter+1)) { #Ahora a las cuentas de la parte positiva del histograma...
  histograma_diff$counts[j]<-histograma_diff$counts[j]-histograma_diff$counts[2*counter - j + 2] 
  histograma_diff$density[j]<-histograma_diff$density[j]-histograma_diff$density[2*counter - j + 2]
} # les restamos las de la parte simetrica negativa

plot(histograma_diff, main = "IR: [SF + AGN] - symmetric contribution",
     xlab = "S (IR, normalized)", freq = FALSE, col = "orange") # (Vaya potra)
xdens <- seq(-5, 5, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = 0, sd = 1), col = "black", lwd = 2)#gaussiana
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 1"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))







##############################     RADIO     ###################################

#Se procede igual para el radio

lumradio<-luminosities$Radio.luminosity.erg.
#lumradio<-lumradio*(1.4*10^9) #pasamos de erg/s/Hz a erg/s (nuestros catálogos de radio operan a 1.4 GHz)
sfr<-SFR1
loglumradio<-log10(lumradio)
logsfr<-log10(sfr)
plot(loglumradio, logsfr, xlim = c(27, 30.5), pch = 16, col = "green", main = "Radio", 
     xlab = "log [L(Radio) (erg/s)]", ylab = "log [SFR (solar mass/yr)]")

ajuste2<- lm(logsfr ~ loglumradio)
ord2<-as.numeric(ajuste2$coefficients[1])
pend2<-as.numeric(ajuste2$coefficients[2])
abline(a=ord2, b=pend2, lw=2, col="blue")
abline(a=log10(3.8*10^-29), b=1, lw=2, col="black")
legend("topleft", legend = c("Otí-Floranes", "Our fit"), fill = c("black", "blue"))




#Cálculo de la S:

S_radio<-as.numeric()#cálculo de la S
for (i in 1:len1){
  S_radio[i]<-(ord2+pend2*loglumradio[i]-logsfr[i])
}
S_radio_raw<-S_radio
S_radio<-na.omit(S_radio) #na.omit es para quitar NaN, NA y demás purria
mean_radio<-mean(S_radio) #media de S_radio
sig2<-sqrt(var(S_radio)) #desviación típica de S_radio





#Ahora tomamos los datos en crudo, luego centramos la variable S en la moda, y finalmente
#dividimos la variable centrada entre una desviación típica (lo que viene siendo tipificar
#una variable en estadística, vamos, sólo que en este caso estamos utilizando la moda en lugar
#de la media como medida de centralización: S(normalizada) = (S - moda)/sigma).
#Histogramamos cada paso:

histograma2<-hist(S_radio, breaks = seq(-0.95, 1.55, 0.1), freq=FALSE, col = "green", 
                  main = "Histogram of S (Radio)", xlab = "S (Radio)")#histograma de los datos en crudo
moda_radio <- histograma2$mids[which.max(histograma2$counts)] #cálculo de la moda usando el histograma





# -------- Normalizaciones (estimaciones de sigma, calculamos 3 distintas): ----------

sigma_moda_radio <- sqrt(sum((S_radio - moda_radio)^2)/(length(S_radio) - 1))
#definición tradicional de la desv. típica pero usando la moda en lugar de la media

sigma_percentil_radio <- as.numeric((quantile(S_radio, probs = 0.84) - quantile(S_radio, probs = 0.16))/2) 
#dev. típica estimada mediante el rango intercuartílico (prescripción de los percentiles)

S_radio_izq <- as.numeric()
for (i in 1:length(S_radio)){
  if (S_radio[i] <= moda_radio)
    S_radio_izq <- c(S_radio_izq, S_radio[i])
}
sigma_izquierda_radio <- sqrt(sum((S_radio_izq - moda_radio)^2)/(length(S_radio_izq) - 1))
#estimación de la desv. típica usando sólo la parte a la izquierda de la moda

# --------------------------------------------------------------------------------------------------





#Histogramamos ahora la VARIABLE CENTRADA y le plantamos la correspondiente gaussiana por encima

S_radio_centrada <- S_radio - moda_radio
histograma2_cent<-hist(S_radio_centrada, breaks = seq(-0.85, 1.65, 0.1), freq=FALSE, col = "green",
                       main = "Histogram of S (Radio, mode-centered)", xlab = "S (Radio, mode-centered)")
#los "breaks" están trucados para que haya un "mid" (es decir, el centro de una de las barras)
#en el 0 del histograma. De esta manera la moda coincide con el 0.
xdens <- seq(-2, 2, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = 0, sd = sigma_izquierda_radio), col = "black", lwd = 2) #gaussiana 
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 0.238"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))





#Histogramamos ahora la VARIABLE TIPIFICADA y le plantamos la correspondiente gaussiana por encima

S_radio_norm <- S_radio_centrada/sigma_izquierda_radio
histograma2_norm <- hist(S_radio_norm, 
                         breaks = seq(-0.85/sigma_izquierda_radio, 1.65/sigma_izquierda_radio, 
                                      0.1/sigma_izquierda_radio), 
                         main = "Histogram of S (Radio, normalized)", 
                         freq=FALSE, col = "green", xlab = "S (Radio, normalized)")
xdens <- seq(-7, 7, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = 0, sd = 1), col = "black", lwd = 2)#gaussiana
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 1"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))

#lines(density(S_radio_norm))
#polygon(density(S_radio_norm), col = "skyblue", border = "black")



# ----------------------------------------------------------------------------------
#Construimos ahora un histograma simétrico a la parte SF de la izquierda de la moda
#Es decir, el siguiente histograma será sólo la parte SF (moda y su izquierda)
#y su reflejo en la parte a la derecha de la moda
# ----------------------------------------------------------------------------------

histograma_sim2<-histograma2_norm
longitud2<-length(histograma_sim2$mids)
counter2 <- 0
for (i in 1:longitud2){
  if (histograma_sim2$mids[i]<0){
    counter2<-counter2+1
  }
}

for (j in (counter2+2):(2*counter2+1)) {
  histograma_sim2$counts[j]<-histograma_sim2$counts[2*counter2 - j + 2] 
  histograma_sim2$density[j]<-histograma_sim2$density[2*counter2 - j + 2]
}

for (j in (2*counter2+2):(longitud2)) {
  histograma_sim2$counts[j]<-0
  histograma_sim2$density[j]<-0
}

plot(histograma_sim2, freq = FALSE, col = "green", main = "Histogram of S (Radio, SF + symmetric)",
     xlab = "S (Radio, symmetric)")





# ---------------------------------------------------------------------------------------------------------
#Ahora podemos jugar a poner juntos ambos histogramas: el histograma total, con SF + AGN, y éste simétrico
#último que sería una especie de "sólo SF". Así tendremos una idea del exceso correspondiente a AGN.
# ---------------------------------------------------------------------------------------------------------

library(scales)
plot(histograma2_norm, freq = FALSE,  col='skyblue', border=F,
     main = "Histogram of S (Radio)", xlab = "S (Radio)")
plot(histograma_sim2, freq = FALSE, add=T, col=scales::alpha('red', .5), border=F)
legend("topright", legend = c("SF + symmetric", "SF + AGN", "Overlap"), 
       fill = c(scales::alpha('red', .35), "skyblue", scales::alpha('red', .5)))

kernel2<-density(S_radio_norm)
plot(kernel2, main = "Kernel density of S (Radio, normalized)")
polygon(kernel2, col = "skyblue", border = "black")





# ---------------------------------------------------------------------------------------------------------
# El siguiente paso va a ser histogramar exactamente lo mismo pero restando a la parte positiva del histograma
# las cuentas de la parte simétrica negativa. Así tendremos una idea del exceso correspondiente a AGN.
# ---------------------------------------------------------------------------------------------------------


histograma2_diff<-histograma2_norm
longitud2<-length(histograma2_diff$mids) #longitud es el número de intervalos del histograma
counter2<-0 #counter van a ser el número de breaks que tenemos a la izquierda del cero
for (i in 1:longitud2){
  if (histograma2_diff$mids[i]<0){
    counter2<-counter2+1
  }
}

for (j in (counter2+2):(2*counter2 + 1)) { #Ahora a las cuentas de la parte positiva del histograma...
  histograma2_diff$counts[j]<-histograma2_diff$counts[j]-histograma2_diff$counts[2*counter2 - j + 2] 
  histograma2_diff$density[j]<-histograma2_diff$density[j]-histograma2_diff$density[2*counter2 - j + 2]
} # les restamos las de la parte simetrica negativa. Estas últimas líneas de código
# parecen muy liosas, pero con el histograma delante se puede comprobar que esto funciona

plot(histograma2_diff, freq = FALSE, col = "green",
     main = "Radio: [SF + AGN] - symmetric contribution", xlab = "S (Radio, normalized)") 
# Aquí hemos representado lo que hemos hecho (vaya potra lo del plot de un histograma)...

xdens <- seq(-7, 7, 0.02) #intervalo simétrico para la gaussiana
lines(xdens, dnorm(xdens, mean = 0, sd = 1), col = "black", lwd = 2)#gaussiana
legend("topright", legend = c("Gaussian fit", "mean = 0", "sd = 1"), lty = c(1, 1, 1), 
       col = c("black", "white", "white"))







##############################     X-RAYS     ###################################

lumx<-luminosities$X.Ray.luminosity.erg.s.
sfr<-SFR1
loglumx<-log10(lumx)
logsfr<-log10(sfr)
plot(loglumx, logsfr, main = "X-Rays", xlim = c(38, 43), ylim = c(-0.5, 1.5), pch = 16, col = "red",
     xlab = "log [L(X-Rays) (erg/s)]", ylab = "log [SFR (solar mass/yr)]")

ajuste3<- lm(logsfr ~ loglumx)
ord3<-as.numeric(ajuste3$coefficients[1])
pend3<-as.numeric(ajuste3$coefficients[2])
abline(a=ord3, b=pend3, lw=2, col="blue")
abline(a=log10(3.51*2*10^-41), b=1, lw=2, col="black")
legend("topleft", legend = c("Otí-Floranes, soft range-only (0.2 - 2 keV)", "Our fit"), fill = c("black", "blue"))

#Cálculo de la S:

S_x<-as.numeric()#cálculo de la S
for (i in 1:len1){
  S_x[i]<-(ord3+pend3*loglumx[i]-logsfr[i])
}
S_x_raw<-S_x
S_x<-na.omit(S_x) #na.omit es para quitar NaN, NA y demás purria
mean_x<-mean(S_x) #media de S_radio
sig3<-sqrt(var(S_x)) #desviación típica de S_radio


# Normalizamos:

hist(S_x, freq = FALSE, breaks = 6, col = "darkblue")

kernel3<-density(S_x)
plot(kernel3, main = "Kernel density of S (X)")
polygon(kernel3, col = "darkblue", border = "black")







############################    Pendiente que minimiza dispersión de los datos     ##############################

histograma_norm <- hist(S_IR_norm, 
                        breaks = seq(-0.75/sigma_izquierda_IR, 1.35/sigma_izquierda_IR, 0.10/sigma_izquierda_IR), 
                        main = "Histogram and KDE of S (IR, normalized)", freq=FALSE, 
                        col = rgb(0.25,0.25,0.25,0.5), xlab = "S (IR, normalized)")
kernel<-density(S_IR_norm)
polygon(kernel, col = rgb(0,0,1,0.4), border = "black")
legend("topright", legend = c("Histogram", "KDE"), 
       fill = c(rgb(0.25,0.25,0.25,0.5), rgb(0,0,1,0.4)))


histograma2_norm <- hist(S_radio_norm, 
                         breaks = seq(-0.85/sigma_izquierda_radio, 1.65/sigma_izquierda_radio, 
                                      0.1/sigma_izquierda_radio), 
                         main = "Histogram and KDE of S (Radio, normalized)", 
                         freq=FALSE, col = rgb(0.25,0.25,0.25,0.5), xlab = "S (Radio, normalized)")
kernel2<-density(S_radio_norm)
polygon(kernel2, col = rgb(0,0,1,0.4), border = "black")
legend("topright", legend = c("Histogram", "KDE"), 
       fill = c(rgb(0.25,0.25,0.25,0.5), rgb(0,0,1,0.4)))

###########################################  KDE-only  ######################################################


#Lo hacemos para el IR:

kernel_S_IR <- density(S_IR)
max_IR <- which.max(kernel_S_IR$y)
maximo_IR <- kernel_S_IR$x[max_IR]
kernel_S_IR_centrado <- kernel_S_IR
kernel_S_IR_centrado$x <- kernel_S_IR_centrado$x - maximo_IR
kernel_S_IR_norm <- kernel_S_IR_centrado
kernel_S_IR_norm$x <- kernel_S_IR_norm$x/sigma_izquierda_IR
kernel_S_IR_norm$y <- kernel_S_IR_norm$y*sigma_izquierda_IR
plot(kernel_S_IR_norm, main = "Kernel density of S (IR, normalized)")
polygon(kernel_S_IR_norm, col = "red", border = "black")


len_ker_Ir <- length(kernel_S_IR_norm$y)
kernel_S_IR_sym <- kernel_S_IR_norm
for (i in (max_IR + 1):(2*max_IR - 1)){
  kernel_S_IR_sym$y[i] <- kernel_S_IR_sym$y[(2*max_IR) - 1 - i + 1]
}
for (i in (2*max_IR):(len_ker_Ir)){
  kernel_S_IR_sym$y[i] <- 0
}

plot(kernel_S_IR_sym, main = "Kernel density of S (IR, symmetric)")
polygon(kernel_S_IR_sym, col = "red", border = "black")

kernel_S_IR_diff <- kernel_S_IR_norm

for (i in 1:len_ker_Ir){
  kernel_S_IR_diff$y[i] <- kernel_S_IR_norm$y[i] - kernel_S_IR_sym$y[i]
}

plot(kernel_S_IR_diff, ylim = c(0, 0.37), main = "Kernel density of S (IR, difference)")
polygon(kernel_S_IR_diff, col = "red", border = "black")


#Ahora para el RADIO:

kernel_S_radio <- density(S_radio)
max_radio <- which.max(kernel_S_radio$y)
maximo_radio <- kernel_S_radio$x[max_radio]
kernel_S_radio_centrado <- kernel_S_radio
kernel_S_radio_centrado$x <- kernel_S_radio_centrado$x - maximo_radio
kernel_S_radio_norm <- kernel_S_radio_centrado
kernel_S_radio_norm$x <- kernel_S_radio_norm$x/sigma_izquierda_radio
kernel_S_radio_norm$y <- kernel_S_radio_norm$y*sigma_izquierda_radio
plot(kernel_S_radio_norm, main = "Kernel density of S (Radio, normalized)")
polygon(kernel_S_radio_norm, col = "skyblue", border = "black")

len_ker_radio <- length(kernel_S_radio_norm$y)
kernel_S_radio_sym <- kernel_S_radio_norm
for (i in (max_radio + 1):(2*max_radio - 1)){
  kernel_S_radio_sym$y[i] <- kernel_S_radio_sym$y[(2*max_radio) - 1 - i + 1]
}
for (i in (2*max_radio):(len_ker_radio)){
  kernel_S_radio_sym$y[i] <- 0
}

plot(kernel_S_radio_sym, main = "Kernel density of S (Radio, symmetric)")
polygon(kernel_S_radio_sym, col = "skyblue", border = "black")

kernel_S_radio_diff <- kernel_S_radio_norm

for (i in 1:len_ker_radio){
  kernel_S_radio_diff$y[i] <- kernel_S_radio_norm$y[i] - kernel_S_radio_sym$y[i]
}

plot(kernel_S_radio_diff, ylim = c(0, 0.37), main = "Kernel density of S (Radio, difference)")
polygon(kernel_S_radio_diff, col = "skyblue", border = "black")

#Ahora para el X:

kernel_S_x <- density(S_x)
max_x <- which.max(kernel_S_x$y)
maximo_x <- kernel_S_x$x[max_x]
kernel_S_x_centrado <- kernel_S_x
kernel_S_x_centrado$x <- kernel_S_x_centrado$x - maximo_x
kernel_S_x_norm <- kernel_S_x_centrado

S_x_izq <- as.numeric()
for (i in 1:length(S_x)){
  if (S_x[i] <= maximo_x)
    S_x_izq <- c(S_x_izq, S_x[i])
}
sigma_izquierda_x <- sqrt(sum((S_x_izq - maximo_x)^2)/(length(S_x_izq) - 1))
#dev. típica estimada usando sólo la parte de los datos a la izquierda de la moda
S_x_centrada <- S_x - maximo_x
S_x_norm <- S_x_centrada/sigma_izquierda_x

kernel_S_x_norm$x <- kernel_S_x_norm$x/sigma_izquierda_x
kernel_S_x_norm$y <- kernel_S_x_norm$y*sigma_izquierda_x
plot(kernel_S_x_norm, main = "Kernel density of S (X, normalized)")
polygon(kernel_S_x_norm, col = "darkblue", border = "black")


len_ker_x <- length(kernel_S_x_norm$y)
kernel_S_x_sym <- kernel_S_x_norm
for (i in (max_x + 1):(2*max_x - 1)){
  kernel_S_x_sym$y[i] <- kernel_S_x_sym$y[(2*max_x) - 1 - i + 1]
}
for (i in (2*max_x):(len_ker_x)){
  kernel_S_x_sym$y[i] <- 0
}

plot(kernel_S_x_sym, main = "Kernel density of S (X, symmetric)")
polygon(kernel_S_x_sym, col = "darkblue", border = "black")

kernel_S_x_diff <- kernel_S_x_norm

for (i in 1:len_ker_x){
  kernel_S_x_diff$y[i] <- kernel_S_x_norm$y[i] - kernel_S_x_sym$y[i]
}

plot(kernel_S_x_diff, main = "Kernel density of S (X, difference)")
polygon(kernel_S_x_diff, col = "darkblue", border = "black")


#KERNEL DENSITIES SUPERIMPOSED:

#Para el IR:

plot(kernel_S_IR_norm, lty = 1, lwd = 2, cex = 2, col = "red", main = "Kernel densities of S (IR, normalized)")
lines(kernel_S_IR_sym, lty = 2, lwd = 2, cex = 2, col = "red")
lines(kernel_S_IR_diff, lty = 3, lwd = 2, cex = 2, col = "red")

leg1 <- expression(SF + AGN, Symmetric, Difference)
legend('topright', col = "red", leg1, lty = c(1,2,3))

#Ahora el Radio:

plot(kernel_S_radio_norm, lty = 1, lwd = 2, cex = 2, col = "blue", 
     main = "Kernel densities of S (Radio, normalized)")
lines(kernel_S_radio_sym, lty = 2, lwd = 2, cex = 2, col = "blue")
lines(kernel_S_radio_diff, lty = 3, lwd = 2, cex = 2, col = "blue")

leg2 <- expression(SF + AGN, Symmetric, Difference)
legend('topright', col = "blue", leg2, lty = c(1,2,3))

#Ahora el X:

plot(kernel_S_x_norm, lty = 1, lwd = 2, cex = 2, col = "black", 
     main = "Kernel densities of S (X, normalized)")
lines(kernel_S_x_sym, lty = 2, lwd = 2, cex = 2, col = "black")
lines(kernel_S_x_diff, lty = 3, lwd = 2, cex = 2, col = "black")

leg2 <- expression(SF + AGN, Symmetric, Difference)
legend('topright', col = "black", leg2, lty = c(1,2,3))

#Y ahora... ¡todos juntos! :D

plot(kernel_S_IR_norm, lty = 1, lwd = 2, cex = 2, col = "orange", 
     main = "Kernel densities of S (IR + Radio, normalized)", xlab = "S (IR + Radio, normalized)")
lines(kernel_S_IR_sym, lty = 2, lwd = 2, cex = 2, col = "orange")
lines(kernel_S_IR_diff, lty = 3, lwd = 2, cex = 2, col = "orange")
lines(kernel_S_radio_norm, lty = 1, lwd = 2, cex = 2, col = "skyblue")
lines(kernel_S_radio_sym, lty = 2, lwd = 2, cex = 2, col = "skyblue")
lines(kernel_S_radio_diff, lty = 3, lwd = 2, cex = 2, col = "skyblue")

leg3 <- expression(SF + AGN, Symmetric, Difference, IR, Radio)
legend('topright', col = "black", leg3, lty = c(1,2,3,NA,NA), fill = c(NA, NA, NA, "orange", "skyblue"),
       border = NA)


# GRACIAS STACKOVERFLOW:

#Distribuciones superpuestas para el IR:

## calculate the range of the graph
xlim <- range(kernel_S_IR_norm$x, kernel_S_IR_sym$x, kernel_S_IR_diff$x)
ylim <- range(0, kernel_S_IR_norm$y, kernel_S_IR_sym$y, kernel_S_IR_diff$y)
#pick the colours
kernel1Col <- rgb(1,0,0,0.25)
kernel2Col <- rgb(0,0,1,0.25)
kernel3Col <- rgb(0,1,0,0.5)
## plot the carrots and set up most of the plot parameters
plot(kernel_S_IR_norm, xlim = xlim, ylim = ylim, xlab = 'S (IR, normalized)',
     main = 'Kernel densities of S (IR, normalized)', 
     panel.first = grid())
#put our density plots in
polygon(kernel_S_IR_norm, density = -1, col = kernel1Col)
polygon(kernel_S_IR_sym, density = -1, col = kernel2Col)
polygon(kernel_S_IR_diff, density = -1, col = kernel3Col)
## add a legend in the corner
legend('topright',c('SF + AGN','Symmetric', 'Difference'),
       fill = c(kernel1Col, kernel2Col, kernel3Col), bty = 'n',
       border = NA)

# Y ahora para el Radio...

## calculate the range of the graph
xlim <- range(kernel_S_radio_norm$x, kernel_S_radio_sym$x, kernel_S_radio_diff$x)
ylim <- range(0, kernel_S_radio_norm$y, kernel_S_radio_sym$y, kernel_S_radio_diff$y)
#pick the colours
kernel4Col <- rgb(1,0,0,0.25)
kernel5Col <- rgb(0,0,1,0.25)
kernel6Col <- rgb(0,1,0,0.5)
## plot the carrots and set up most of the plot parameters
plot(kernel_S_radio_norm, xlim = xlim, ylim = ylim, xlab = 'S (Radio, normalized)',
     main = 'Kernel densities of S (Radio, normalized)',
     panel.first = grid())
#put our density plots in
polygon(kernel_S_radio_norm, density = -1, col = kernel4Col)
polygon(kernel_S_radio_sym, density = -1, col = kernel5Col)
polygon(kernel_S_radio_diff, density = -1, col = kernel6Col)
legend('topright',c('SF + AGN','Symmetric', 'Difference'),
       fill = c(kernel4Col, kernel5Col, kernel6Col), bty = 'n',
       border = NA)

# Y ahora para el X...

## calculate the range of the graph
xlim <- range(kernel_S_x_norm$x, kernel_S_x_sym$x, kernel_S_x_diff$x)
ylim <- range(0, kernel_S_x_norm$y, kernel_S_x_sym$y, kernel_S_x_diff$y)
#pick the colours
kernel7Col <- rgb(1,0,0,0.25)
kernel8Col <- rgb(0,0,1,0.25)
kernel9Col <- rgb(0,1,0,0.5)
## plot the carrots and set up most of the plot parameters
plot(kernel_S_x_norm, xlim = xlim, ylim = ylim, xlab = 'S (X, normalized)',
     main = 'Kernel densities of S (X, normalized)',
     panel.first = grid())
#put our density plots in
polygon(kernel_S_x_norm, density = -1, col = kernel7Col)
polygon(kernel_S_x_sym, density = -1, col = kernel8Col)
polygon(kernel_S_x_diff, density = -1, col = kernel9Col)
legend('topright',c('SF + AGN','Symmetric', 'Difference'),
       fill = c(kernel7Col, kernel8Col, kernel9Col), bty = 'n',
       border = NA)

############################## COLOUR-RAMPPALETTE ###################################

S_IR_raw_centered <- S_IR_raw - maximo_IR
S_IR_raw_norm<-S_IR_raw_centered/sigma_izquierda_IR
S_IR_raw_norm<-as.numeric(S_IR_raw_norm)


S_radio_raw_centered <- S_radio_raw - maximo_radio
S_radio_raw_norm<-S_radio_raw_centered/sigma_izquierda_radio
S_radio_raw_norm<-as.numeric(S_radio_raw_norm)

S_x_raw_centered <- S_x_raw - maximo_x
S_x_raw_norm<-S_x_raw_centered/sigma_izquierda_x
S_x_raw_norm<-as.numeric(S_x_raw_norm)

paso <- kernel_S_IR_norm$x[2] - kernel_S_IR_norm$x[1]
paso <-paso/2
paso2 <- kernel_S_radio_norm$x[2] - kernel_S_radio_norm$x[1]
paso2<-paso2/2
paso3 <- kernel_S_x_norm$x[2] - kernel_S_x_norm$x[1]
paso3<-paso3/2

#------------------
#Código para el IR:
#------------------

D_S_IR<-as.numeric()
D_S_IR_sym<-as.numeric()

for (i in 1:len1){
  if (is.na(S_IR_raw_norm[i])==TRUE){
    S_IR_raw_norm[i]<-10000
  }
}

for (i in 1:len1){
  for (j in 1:len_ker_Ir){
    if (abs(S_IR_raw_norm[i] - kernel_S_IR_norm$x[j]) < paso){
      D_S_IR[i]<-kernel_S_IR_norm$y[j]
      D_S_IR_sym[i]<-kernel_S_IR_sym$y[j]
    }
  }
}

p_SF_IR<-D_S_IR_sym/D_S_IR
p_AGN_IR<-(1-p_SF_IR)

for (k in 1:23){ 
  D_S_IR<-c(D_S_IR, NA)
  D_S_IR_sym<-c(D_S_IR_sym, NA)
  p_SF_IR<-c(p_SF_IR, NA)
  p_AGN_IR<-c(p_AGN_IR, NA)
} 


plot(loglumir, logsfr, type = "n", xlim = c(-3,2), pch = 16, col = "orange", main = "IR (22.1 microns)",
     xlab = "log [L(22 microns)/10^43 (erg/s)]", ylab = "log [SFR (solar mass/yr)]") 
library(colorRamps)
matcol<-matlab.like(11)
for (i in 1:len1){
  if (is.na(p_AGN_IR[i])==FALSE){
    points(loglumir[i], logsfr[i], pch = 16, col =  matcol[round(p_AGN_IR[i]*10, digits = 0)+1])
  }
}
abline(a=0.377, b=0.564, lw=2, col="black") #Aquí iría una línea con los valores de los autores
abline(a=ord, b=pend, lw=2, col="magenta") #representamos la recta del ajuste lineal
legend("topleft", legend = c("Catalán-Torrecilla et al.", "Our fit"), fill = c("black", "magenta"))

#---------------
#AHORA EL RADIO:
#---------------

D_S_radio<-as.numeric()
D_S_radio_sym<-as.numeric()

for (i in 1:len1){
  if (is.na(S_radio_raw_norm[i])==TRUE){
    S_radio_raw_norm[i]<-10
  }
}

for (i in 1:len1){
  for (j in 1:len_ker_radio){
    if (abs(S_radio_raw_norm[i] - kernel_S_radio_norm$x[j]) < paso2){
      D_S_radio[i]<-kernel_S_radio_norm$y[j]
      D_S_radio_sym[i]<-kernel_S_radio_sym$y[j]
    }
  }
}

p_SF_radio<-D_S_radio_sym/D_S_radio
p_AGN_radio<-(1-p_SF_radio)

for (k in 1:23){ 
  D_S_radio<-c(D_S_radio, NA)
  D_S_radio_sym<-c(D_S_radio_sym, NA)
  p_SF_radio<-c(p_SF_radio, NA)
  p_AGN_radio<-c(p_AGN_radio, NA)
} 


plot(loglumradio, logsfr, type = "n", xlim = c(27, 30.5), pch = 16, col = "orange", main = "Radio",
     xlab = "log [L(Radio) (erg/s/Hz)]", ylab = "log [SFR (solar mass/yr)]") 
library(colorRamps)
matcol<-matlab.like(11)
for (i in 1:len1){
  if (is.na(p_AGN_radio[i])==FALSE){
    points(loglumradio[i], logsfr[i], pch = 16, col =  matcol[round(p_AGN_radio[i]*10, digits = 0)+1])
  }
}
abline(a=ord2, b=pend2, lw=2, col="magenta")
abline(a=log10(3.8*10^-29), b=1, lw=2, col="black")
legend("topleft", legend = c("Otí-Floranes et al.", "Our fit"), fill = c("black", "magenta"))

#-----------
#AHORA EL X:
#-----------

D_S_x<-as.numeric()
D_S_x_sym<-as.numeric()

for (i in 1:len1){
  if (is.na(S_x_raw_norm[i])==TRUE){
    S_x_raw_norm[i]<-10
  }
}

for (i in 1:len1){
  for (j in 1:len_ker_x){
    if (abs(S_x_raw_norm[i] - kernel_S_x_norm$x[j]) < paso3){
      D_S_x[i]<-kernel_S_x_norm$y[j]
      D_S_x_sym[i]<-kernel_S_x_sym$y[j]
    }
  }
}

p_SF_x<-D_S_x_sym/D_S_x
p_AGN_x<-(1-p_SF_x)

for (k in 1:135){ 
  D_S_x<-c(D_S_x, NA)
  D_S_x_sym<-c(D_S_x_sym, NA)
  p_SF_x<-c(p_SF_x, NA)
  p_AGN_x<-c(p_AGN_x, NA)
} 

plot(loglumx, logsfr, main = "X-Rays", type = 'n', xlim = c(38, 43), ylim = c(-0.5, 1.5), pch = 16, col = "red",
     xlab = "log [L(X-Rays) (erg/s)]", ylab = "log [SFR (solar mass/yr)]")
library(colorRamps)
matcol<-matlab.like(11)
for (i in 1:len1){
  if (is.na(p_AGN_x[i])==FALSE){
    points(loglumx[i], logsfr[i], pch = 16, col =  matcol[round(p_AGN_x[i]*10, digits = 0)+1])
  }
}
abline(a=ord3, b=pend3, lw=2, col="magenta")
#abline(a=log10(3.51*2*10^-41), b=1, lw=2, col="black")
#legend("topleft", legend = c("Otí-Floranes, soft range-only (0.2 - 2 keV)", "Our fit"), fill = c("black", "magenta"))
abline(v=42, lw=2, col="black")
legend("topleft", legend = c("Traditional criterion", "Our fit"), fill = c("black", "magenta"))


############################ PROBABILITY PLOTS ###################################
#par(mfrow=c(3,1))
#IR:
xnew<-seq(min(S_IR_norm), max(S_IR_norm), (max(S_IR_norm)-min(S_IR_norm))/511)
psf<-kernel_S_IR_sym$y/kernel_S_IR_norm$y
pagn<-(1-psf)
plot(xnew, psf, pch = 16, col = "skyblue",
     main = "SF and AGN probabilities (IR)", xlab = "S (IR, normalized)", ylab = "Probability")
points(xnew, pagn, pch = 16, col = "orange")
legend('right', 0.2, legend = c("SF probability", "AGN probability"), 
       fill = c("skyblue", "orange"))

#RADIO:
xnew2<-seq(min(S_radio_norm), max(S_radio_norm), (max(S_radio_norm)-min(S_radio_norm))/511)
psf2<-kernel_S_radio_sym$y/kernel_S_radio_norm$y
pagn2<-(1-psf2)
plot(xnew2, psf2, pch = 16, col = "skyblue",
     main = "SF and AGN probabilities (Radio)", xlab = "S (Radio, normalized)", ylab = "Probability")
points(xnew2, pagn2, pch = 16, col = "orange")
legend('right', 0.2, legend = c("SF probability", "AGN probability"), 
       fill = c("skyblue", "orange"))

#X-RAYS:
xnew3<-seq(min(S_x_norm), max(S_x_norm), (max(S_x_norm)-min(S_x_norm))/511)
psf3<-kernel_S_x_sym$y/kernel_S_x_norm$y
pagn3<-(1-psf3)
plot(xnew3, psf3, pch = 16, col = "skyblue",
     main = "SF and AGN probabilities (X-rays)", xlab = "S (X-rays, normalized)", ylab = "Probability")
points(xnew3, pagn3, pch = 16, col = "orange")
legend('right', 0.2, legend = c("SF probability", "AGN probability"), 
       fill = c("skyblue", "orange"))

####################################################################################

h3<-hist(S_x_norm, freq = FALSE, breaks = seq(-2.5,4.5, 1), col = rgb(0.25,0.25,0.25,0.5),
         main = "Histogram and KDE of S (X-rays, normalized)", xlab = "S (X-rays, normalized)")
polygon(kernel_S_x_norm, col = rgb(0,0,1,0.4), border = "black")
legend("topleft", legend = c("Histogram", "KDE"), 
       fill = c(rgb(0.25,0.25,0.25,0.5), rgb(0,0,1,0.4)))