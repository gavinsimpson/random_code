## For GAMM models from mgcv:::gamm

## Model Checking function
tsDiagGamm <- function(x, timevar, observed, f = 0.3, type = "normalized") {
    resi <- resid(x$lme, type = type)
    fits <- fitted(x$lme)
    on.exit(layout(1))
    layout(matrix(1:6, ncol = 3, byrow = TRUE))
    plot(resi ~ fits, ylab = "Normalized Residuals",
         xlab = "Fitted Values", main = "Fitted vs. Residuals")
    lines(lowess(x = fits, y = resi, f = f), col = "blue",
          lwd = 2)
    plot(resi ~ timevar, ylab = "Normalized Residuals",
         xlab = "Time", main = "Time series of residuals")
    lines(lowess(x = timevar, y = resi, f = f), col = "blue", lwd = 2)
    plot(observed ~ fits, ylab = "Observed",
         xlab = "Fitted Values", main = "Fitted vs. Observed",
         type = "n")
    abline(a = 0, b = 1, col = "red")
    points(observed ~ fits)
    lines(lowess(x = fits, y = observed, f = f), col = "blue",
          lwd = 2)
    hist(resi, freq = FALSE, xlab = "Normalized Residuals")
    qqnorm(resi)
    qqline(resi)
    acf(resi, main = "ACF of Residuals")
}
