source("main_functions.R")
library(readxl)

DATA_FILE <- "Time-to-transition.xlsx"

T_DATES   <- as.Date(c("2031-01-01", "2041-01-01"))   # T1, T2
TAU_DATES <- as.Date(c("2024-01-01", "2036-09-04"))   # tau1, tau2

FIX_B_PLUS     <- c(NA, 1)     
FIX_THETA_PLUS <- c(NA, 1)

START_RDCM <- c(a = 0.03*365, b_minus = 0.01, b_plus = 1, 
                theta_minus = 0.02 * sqrt(365), theta_plus = 1)


RUN_SRDCM        <- TRUE
RUN_BOOT_RDCM    <- FALSE
RUN_BOOT_SRDCM   <- FALSE
N_BOOT           <- 1000

dat <- read_excel(DATA_FILE, col_types = c("date", "numeric"))
y    <- dat$`Time-to-transition`
date <- as.Date(dat$Date)


Ti    <- as.numeric(date - date[1])/365
TimeF <- as.numeric(T_DATES - date[1])/365
tau   <- as.numeric(TAU_DATES - date[1])/365
X     <- y

fit1 <- fit_rdcm(X, Ti, TimeF[1], tau[1], start = START_RDCM,
                 fix_b_plus = FIX_B_PLUS, fix_theta_plus = FIX_THETA_PLUS,
                 lower = c(0.01 * 365, 0.01,0.5, 0.001*sqrt(365), 0.5), 
                 upper = c(0.1 * 365, 0.1, 100, 0.1 * sqrt(365), 20))



print(fit1)

z1 <- rdcm_residuals(fit1)
hist(z1, freq = FALSE, breaks = 100, main = "", xlab = "", ylab = "",
     col = "gray90", border = "gray60", xlim = c(-6, 8))
curve(dnorm(x), col = "red", lwd = 2, add = TRUE)

print(ks.test(z1, "pnorm"))


if (RUN_BOOT_RDCM) {
  boot1 <- bootstrap_rdcm(fit1, B = N_BOOT, seed = 1, save_file = "boot_RDCM.rds")
  cat(N_BOOT, "\n")
  print(param_table(fit1, boot1), digits = 4)
}



if (RUN_SRDCM) {
  
  p1 <- fit1$par
  
  start2 <- c(a = c(p1["a1"], p1["a1"]),
              b_minus = rep(p1["b_minus1"], 2), b_plus = rep(p1["b_plus1"], 2),
              theta_minus = c(p1["theta_minus1"], p1["theta_minus1"]),
              theta_plus = rep(p1["theta_plus1"], 2))
  
  fit2 <- fit_srdcm(X, Ti, TimeF, tau, start = start2,
                    pi0 = c(0.5, 0.5), Q = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
                    fix_b_plus = FIX_B_PLUS, fix_theta_plus = FIX_THETA_PLUS,
                    lower = c(0.01*365, 0.001*365, 0.01, 0.01, 0.5, 0.5, 0.001*sqrt(365), 
                              0.00001*sqrt(365), 0.5, 0.5),
                    upper = c(0.2*365, 0.1*365,   0.15, 0.1,  50,  50,  0.1*sqrt(365),   
                              0.01*sqrt(365),    50,   50 ))
  
  print(fit2)
  
  z2 <- srdcm_residuals(fit2)
  hist(z2, freq = FALSE, breaks = 100, main = "", xlab = "", ylab = "",
       col = "gray90", border = "gray60", xlim = c(-6, 8))
  curve(dnorm(x), col = "red", lwd = 2, add = TRUE)
  
  
  print(ks.test(z2, "pnorm"))
  
  g1 <- fit2$gamma[, 1]
  col6 <- ifelse(g1 >= 0.75, "red", ifelse(g1 >= 0.5, "green",
                                           ifelse(g1 > 0.25, "blue", "black")))
  plot(date[-1], X[-1], type = "h", col = col6, xlab = "", ylab = expression(X[t]))
  
  if (RUN_BOOT_SRDCM) {
    boot2 <- bootstrap_srdcm(fit2, B = N_BOOT, seed = 1, save_file = "boot_SRDCM.rds")
    cat(N_BOOT, "\n")
    print(param_table(fit2, boot2), digits = 4)
    print(round(boot2$sd[grep("^p", names(boot2$sd))], 4))
  }
}

