activeness2 = function(lum, V_max_inverse, xlabel, colour){
  
  sorted_luminosity = sort.int(lum, na.last = TRUE, index.return = TRUE)
  x_ordenado = na.omit(sorted_luminosity$x)
  N = length(x_ordenado)
  sorted_index = sorted_luminosity$ix[1:N]
  sorted_V_max_inverse = V_max_inverse[sorted_index]
  density_lower = cumsum(sorted_V_max_inverse)
  total_density = density_lower[N]
  density_lower_interp = approxfun(x_ordenado, density_lower, yleft = 0, yright = total_density)
  z = seq(min(x_ordenado), max(x_ordenado), (max(x_ordenado) - min(x_ordenado))/300)
  density_greater = total_density - density_lower_interp(z)
  normal = rep(0, length(z))
  sym = rep(0, length(z))
  pivot_interp = approxfun(density_lower, x_ordenado)
  pivot_point = pivot_interp(total_density/2)
  sigma = pivot_point - pivot_interp(total_density*.16)
  #pivot_point=0
  # pivot_point = moda_color_W2_W1_red
  for (i in 1:length(z)){
    if (z[i] <= pivot_point){
      symmetric = 2*density_lower_interp(pivot_point) - density_lower_interp(z[i])
    }
    else {
      symmetric = density_lower_interp(2*pivot_point - z[i])
    }
    normal[i] = min(symmetric, density_greater[i])
    sym[i] = symmetric
  }
  plot(z, 2*density_lower_interp(pivot_point) - density_lower_interp(z), col='red', log='y', ylim = c(1, 1e9))
  points(z, density_lower_interp(2*pivot_point - z), col='yellow')
  points(z, density_greater, col='black')
  lines(z, normal, col='blue')
  lines(z, sym, col = 'green', lty=2)
  lines(z, total_density*0.5*(1 - erf((z - pivot_point)/sigma)))
  print(sigma)
  print(pivot_point)
  
  active = density_greater - normal
  f_active = active/density_greater
  z_active = z[min(which(f_active >= 0.5))]
  
  plot(z, f_active)
  
  minimo = sigma_minima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
    quantile(na.omit(lum), probs = 0.50)
  maximo = sigma_maxima*(quantile(na.omit(lum), probs = 0.50) - quantile(na.omit(lum), probs = 0.16)) +
    quantile(na.omit(lum), probs = 0.50)
  
  library(ggplot2)
  
  df = data.frame(z, log10(normal))
  df2 = data.frame(z, log10(density_greater))
  df3 =data.frame(z, log10(symmetric))
  # print(ggplot(df, aes(x = df$z)) +
  #         ylim(0, max(na.omit(log10(density_greater)))) +
  #         geom_line(aes(y = df$log10.normal.), color = colour, linetype = 'dashed') +
  #         geom_line(aes(y = df2$log10.density_greater.), color = 'black') +
  #         geom_line(aes(y = df3$log10.symmetric.), color = 'blue') +
  #         geom_rect(aes(xmin = minimo, xmax = maximo, ymin = -Inf, ymax = Inf), color = 'grey', alpha = 0.005) +
  #         geom_vline(xintercept = z_active) +
  #         theme_bw() +
  #         xlab(xlabel) +
  #         ylab('log (N > x)'))
  
  return(z_active)
}
#activeness2(color_W2_W1_excess, V_max_color_W2_W1_inverse, 'Luminosity excess (IR W2 band)', 'red')
#activeness2(L_excess_IR2, V_max_IR2_inverse, 'Luminosity excess (IR W2 band)', 'red')
#activeness2(L_excess_Radio, V_max_Radio_inverse, 'Luminosity excess (IR W2 band)', 'red')
activeness2(color_Radio_IR_excess, V_max_color_Radio_IR_inverse, 'Luminosity excess (IR W2 band)', 'red')