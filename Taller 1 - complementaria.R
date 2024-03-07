.initial_time <- Sys.time() # Cronometrar corrida
set.seed(123) # Fijar semilla

# Taller 1 ----------------------------------------------------------------
pacman::p_load(readxl, dplyr, stats, forecast)

# 1.a ---------------------------------------------------------------------

bonds <-
  read_excel(
    path = file.path("Datos taller 1.xlsx"),
    sheet = "Datos semanales"
  ) %>% as.data.frame()

rownames(bonds) <- bonds[, 1]
# Eliminar la primera columna del dataframe
bonds <- bonds[, -1]

colnames(bonds) <- c(.25, .5, 1, 2, 3, 5, 7, 10)

last_curve <- tail(bonds, 1)


get_flows <- function(last_curve, bond) {
  semesters <- max(1, 2 * (colnames(last_curve)[bond] %>% as.numeric()))
  rate <- last_curve[bond]
  cupon_semesters <- rep(rate * 100 * 182.5 / 360, semesters) %>% unlist()
  names(cupon_semesters) <- seq(from = 1, to = length(cupon_semesters)) / 2
  cupon_semesters[length(cupon_semesters)] <- tail(cupon_semesters, 1) + 100
  if (colnames(last_curve)[bond] %>% as.numeric() == 0.25) {
    names(cupon_semesters) <- 0.25
    cupon_semesters <- rate * 100 * 91.25 / 360 + 100
  }
  cupon_semesters
}

flows <- sapply(c(1:8), function(x) get_flows(last_curve, x))
names(flows) <- paste("Bono de", colnames(bonds), "años")
flows
#' Hasta ahora no hemos hecho nada, solo representar en un objeto bonito (flows)
#' los flujos (y las periodicidades) de los 8 bonos del tesoro

# 1.b ---------------------------------------------------------------------

# Crear matriz
n_alphas <- 7
alphas_names <- paste("Alpha ", c(1:n_alphas))
alphas_iniciales <- runif(n = 6, min = 0.01, max = 0.1)
Desde <- c(0, .25, 1, 2, 3, 5, 7)
Hasta <- c(tail(Desde, length(Desde) - 1), 10)

matriz <- cbind(alphas_iniciales, Desde, Hasta)
rownames(matriz) <- alphas_names
matriz

# Factores de descuento
get_df <- function(matriz, t) {
  # Comprobar si x está en el rango permitido
  if (t < 0 || t > max(matriz[, "Hasta"])) {
    stop("El valor de x debe estar dentro del rango de la matriz.")
  }

  # Inicializar la suma de la integral
  integral <- 0

  # Iterar sobre los intervalos
  for (i in 1:nrow(matriz)) {
    # Obtener los valores de alpha y los intervalos
    alpha <- matriz[i, 1]
    desde <- matriz[i, 2]
    hasta <- matriz[i, 3]

    # Si x está dentro del intervalo actual, sumar la parte proporcional y salir del bucle
    if (t <= hasta) {
      integral <- integral + alpha * (t - desde)
      break
    }

    # Si el intervalo completo está antes de x, sumar toda el área del rectángulo
    integral <- integral + alpha * (hasta - desde)
  }

  df <- exp(-integral)
  df
}

# Un plot que muestra que a mayor T, mayor
datos_de_prueba <- seq(from = 0, to = 10, length.out = 100)
plot(datos_de_prueba, sapply(datos_de_prueba, function(x) get_df(matriz, x)))

# Valores presentes
get_pv_error <- function(flujo, matriz) {
  periodos <- names(flujo) %>% as.numeric()
  factors <- sapply(periodos, function(x) get_df(matriz, x))
  present_value <- sum(factors * flujo)
  error <- (100 - present_value)^2
  error
}

# Crear función objetivo
funcion_objetivo <- function(alphas, flows) {
  matriz <- cbind(alphas, Desde, Hasta)
  objetivo <- lapply(flows, function(x) get_pv_error(x, matriz)) %>%
    unlist() %>%
    sum()
  objetivo
}

environment(funcion_objetivo) <- environment()

# Optimizar

alphas_elegidos <- optim(
  par = rep(0.01, n_alphas), # Los valores iniciales de los parámetros
  fn = funcion_objetivo, # La función objetivo a minimizar
  flows = flows,
  method = "L-BFGS-B", # Método de optimización
  lower = rep(0, 6), # Límites inferiores para los alphas (todos 0, es decir, positivos)
  upper = rep(2, 6) # Límites superiores para los alphas (Infinito)
)

alphas_optimos_7_feb <- alphas_elegidos$par

# Punto 3 ----------------------------------------------------------

encontrar_alphas <- function(row) {
  curva_spot <- bonds[row, ]
  flows <- sapply(c(1:8), function(x) {
    get_flows(curva_spot, x)
  })
  alphas_elegidos <- optim(
    par = rep(0.01, n_alphas),
    # Los valores iniciales de los parámetros
    fn = funcion_objetivo,
    # La función objetivo a minimizar
    flows = flows,
    method = "L-BFGS-B",
    # Método de optimización
    lower = rep(0, 6),
    # Límites inferiores para los alphas (todos 0, es decir, positivos)
    upper = rep(2, 6) # Límites superiores para los alphas (Infinito)
  )
  alphas_elegidos$par
}

serie_alphas <- sapply(seq_len(nrow(bonds)), encontrar_alphas) %>%
  t() %>%
  diff()

# Punto 4 ---------------------------------------------------------------------

pca_result <- prcomp(serie_alphas, scale. = FALSE, center = FALSE)

# Punto 5 -----------------------------------------------------------------

var_explicada <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Encontrar el número de componentes que explican más del 90% de la varianza
var_acumulada <- cumsum(var_explicada)
n_componentes <- which.max(var_acumulada > 0.9)

# Seleccionar los primeros n componentes principales
pca_scores <- pca_result$x[, 1:n_componentes]

# Punto 6 -----------------------------------------------------------------

# Ajustar un modelo ARIMA a cada componente principal
arima_models <- lapply(1:n_componentes, function(i) {
  auto.arima(pca_scores[, i])
})

# Simulo
n_semanas <- 52
n_simulaciones <- 10000

# Simular sendas futuras para cada componente principal
simulaciones <- lapply(arima_models, function(model) {
  matrix(replicate(n_simulaciones, simulate(model, nsim = n_semanas)), ncol = n_simulaciones)
})

# Espacio original
# Preparar un array para almacenar los resultados
resultados <- array(dim = c(n_semanas, ncol(serie_alphas), n_simulaciones))
rotation_matrix <- pca_result$rotation[, 1:n_componentes]

# Transformar las simulaciones de vuelta al espacio original
for (i in 1:n_simulaciones) {
  # Combinar las simulaciones de todos los componentes en una matriz
  sim_matrix <- sapply(simulaciones, function(x) x[, i])

  # Transformar de vuelta al espacio original usando la matriz de carga
  resultados_original_space <- rotation_matrix %*% t(sim_matrix)

  # Agregar los resultados al array
  resultados[, , i] <- t(resultados_original_space)
}

resultados <- apply(resultados, c(2, 3), cumsum)

tasas_simuladas <- sweep(
  t(resultados[n_semanas, , ]),
  2,
  alphas_elegidos$par,
  "+"
)

# Punto 2 -----------------------------------------------------------------

# Función que valora
nocional <- 35

get_error_2 <- function(C, alphas) {
  flujos <- rep(nocional * C, 10)
  flujos[10] <- flujos[10] + nocional
  names(flujos) <- c(seq(1:10))

  matriz <- cbind(alphas, Desde, Hasta) %>% as.data.frame()
  periodos <- names(flujos) %>% as.numeric()
  factors <- sapply(periodos, function(x) get_df(matriz, x))
  present_value <- sum(factors * flujos)

  error <- (present_value - nocional)^2
  return(error)
}
environment(get_error_2) <- environment()

# Punto 2
C_optimo <- optim(
  par = 0.05, # Los valores iniciales de los parámetros
  fn = get_error_2, # La función objetivo a minimizar
  alphas = alphas_elegidos$par,
  method = "L-BFGS-B", # Método de optimización
  lower = 0, # Límites inferiores para los alphas (todos 0, es decir, positivos)
  upper = 100 # Límites superiores para los alphas (Infinito)
)$par

get_error_2(C_optimo, alphas_elegidos$par) # Da 0, entonces quedó bien

# Punto 7 -----------------------------------------------------------------
# Crear los flujos
flujos <- rep(nocional * C_optimo, 10)
flujos[10] <- flujos[10] + nocional
names(flujos) <- c(seq(1:10)) - n_semanas / 52 # OJO

get_pv_1y <- function(C = C_optimo, alphas) {
  matriz <- cbind(alphas, Desde, Hasta) %>% as.data.frame()
  periodos <- names(flujos) %>% as.numeric()
  factors <- sapply(periodos, function(x) get_df(matriz, x))
  present_value <- sum(factors * flujos)

  present_value
}

# Hoy
pv_1Y <- sapply(seq_len(n_simulaciones), function(i) {
  get_pv_1y(alphas = tasas_simuladas[i, ])
}) # En millones
pv_1Y <- ((pv_1Y - nocional) / nocional) * 100 # En porcentaje

# Resultado final
hist(pv_1Y,
  col = "skyblue", # Color of the bars
  border = "white", # Color of the border of the bars
  xlim = c(min(pv_1Y), max(pv_1Y)), # X-axis limits
  main = "Histograma VP", # Main title
  xlab = "Valores", # X-axis label
  ylab = "Frecuencia", # Orientation of axis labels: always parallel to the axis,
  breaks = 100
)

Sys.time() - .initial_time # Reporta cuánto toma corriendo
