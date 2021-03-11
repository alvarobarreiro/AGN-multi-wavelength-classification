x_ordenado = na.omit(sort.int(L_reducida_Hard, na.last = TRUE, index.return = TRUE)$x)
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

plot(z, N_normal, pch = 16, col = 'blue', type = 'n', log = "y", ylab = 'N (> x)', xlab = 'Luminosity excess (Hard X-Rays)')
lines(z, N_normal, lty = 2, col = 'blue')
lines(z, N_greater, pch = 16, col = 'black')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
# lines(z, N_active/N_greater, col = 'red') # 


# z_active = z[which.max(N_active)]
print(z_active)
abline(v = z_active)
abline(v = 2*(quantile(na.omit(L_reducida_Hard), probs = 0.5) - quantile(na.omit(L_reducida_Hard), probs = 0.5))
       + quantile(na.omit(L_reducida_Hard), probs = 0.5),
       lty = 2)
# abline(h = 0.5)

