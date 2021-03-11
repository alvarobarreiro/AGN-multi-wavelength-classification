
###################################### KDEs COLOR ###########################################


ro = density(na.omit(color_W2_W1))
ro$bw = 0.1 * sigma_izquierda_color_W2_W1
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 10), main = "KDE IR color")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'red')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(color_Radio_IR))
ro$bw = 0.1 * sigma_izquierda_color_Radio_IR
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE Radio-IR color")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'green')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(color_X))
ro$bw = 0.1 * sigma_izquierda_colorX
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE X-Rays color")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'blue')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))


####################################### KDEs COLOR EXCESS ################################################


ro = density(na.omit(color_W2_W1_reducido))
ro$bw = 0.1 * sigma_izquierda_color_W2_W1_red
xmax = ro$x[which.max(ro$y)]
#buscar mediana y dibujar simétrica
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 12), main = "KDE IR color excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'red')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(color_Radio_IR_reducido))
ro$bw = 0.1 * sigma_izquierda_color_Radio_IR_red
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE Radio-IR color excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'green')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(color_X_reducido))
ro$bw = 0.1 * sigma_izquierda_colorX_red
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE X-Rays color excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'blue')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))
abline(v = 0)
abline(v = xmax)

###################################### KDEs LUMINOSITY EXCESS #########################################

ro = density(na.omit(L_reducida_IR2))
ro$bw = 0.1 * sigma_izquierda_ir2_excess
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 2), main = "KDE IR luminosity excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'red')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(L_reducida_Radio))
ro$bw = 0.1 * sigma_izquierda_radio_excess
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE Radio luminosity excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'green')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))

ro = density(na.omit(L_reducida_Hard))
ro$bw = 0.1 * sigma_izquierda_hard_excess
xmax = ro$x[which.max(ro$y)]
x_inverted = 2*xmax - ro$x
symmetric = which(x_inverted > xmax)
plot(ro, log = 'y', ylim = c(0.001, 1), main = "KDE Hard X-Rays luminosity excess")
lines(x_inverted[symmetric], ro$y[symmetric], col = 'blue')
# ratio = ro$y/approx(x_inverted[symmetric], ro$y[symmetric], ro$x)$y
# plot(ro$x, ratio, ylim = c(0, 3))
