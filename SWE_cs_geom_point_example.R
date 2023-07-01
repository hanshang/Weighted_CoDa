###########################################
# Compositional Data Analysis (CoDa) model
# An example of model fitting
###########################################

R_square_fit_ex <- function(dat, weighting, geom_weight, fh,
                         ncomp_selection = c("eigen_ratio", "fixed"),
                         forecasting_method)
{
  n_year = nrow(dat)
  n_age = ncol(dat)

  # compute geometric mean with geometrically decaying weights

  if(weighting == "simple")
  {
    wei = rep(1/n_year, n_year)
  }
  else if(weighting == "geom")
  {
    wei = vector("numeric", n_year)
    for(ik in 1:n_year)
    {
      wei[ik] = geom_weight * (1 - geom_weight)^(n_year - ik)
    }
  }

  # compute geometric mean

  geom_mean = vector("numeric", n_age)
  geom_mean = rep(0, n_age)
  for(ik in 1:n_year)
  {
    geom_mean = geom_mean + wei[ik] * log(dat[ik,])
  }
  alpha_x = exp(geom_mean)

  # centered log ratio

  h_x_t = matrix(NA, n_year, n_age)
  for(ik in 1:n_year)
  {
    h_x_t[ik,] = log(dat[ik,]/alpha_x)
  }

  wq = diag(wei)
  h_x_t_weight = wq %*% h_x_t
  SVD_decomp = svd(h_x_t_weight)
  if(ncomp_selection == "eigen_ratio")
  {
    # select ncomp by LRS(JASA, 2020) eigen-value ratio tests

    eigen_values = SVD_decomp$d^2
    lambda_val = eigen_values[which(eigen_values > 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))

    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for(ik in 1:k_max)
    {
      eigen_val_ratio[ik] = lambda_val[ik+1]/lambda_val[ik] * ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 0) + ifelse(lambda_val[ik]/lambda_val[1] < tau, 1, 0)
    }
    ncomp = which.min(eigen_val_ratio)
  }
  else if(ncomp_selection == "fixed")
  {
    ncomp = 6
  }
  if(class(ncomp) != "numeric")
  {
    ncomp = 1
  }

  basis = as.matrix(SVD_decomp$v[,1:ncomp])
  score = h_x_t %*% basis
  recon = basis %*% t(score)

  # reconstruction (model in-sample fitting)

  f_x_t_star_recon = d_x_t_star_recon = matrix(NA, n_age, n_year)
  for(ik in 1:n_year)
  {
    f_x_t_star_recon[,ik] = exp(recon[,ik])/sum(exp(recon[,ik]))
    d_x_t_star_recon[,ik] = (f_x_t_star_recon[,ik] * alpha_x)/sum(f_x_t_star_recon[,ik] * alpha_x)
  }
  R2 = 1 - sum((t(d_x_t_star_recon) - dat)^2)/sum((dat - colMeans(dat))^2)

  # forecasts of principal component scores

  score_fore = matrix(NA, ncomp, fh)
  for(ik in 1:ncomp)
  {
    if(forecasting_method == "RWF_no_drift")
    {
      score_fore[ik,] = rwf(as.numeric(score[,ik]), h = fh, drift = FALSE)$mean
    }
    else if(forecasting_method == "RWF_drift")
    {
      score_fore[ik,] = rwf(as.numeric(score[,ik]), h = fh, drift = TRUE)$mean
    }
    else if(forecasting_method == "ETS")
    {
      score_fore[ik,] = forecast(ets(as.numeric(score[,ik])), h = fh)$mean
    }
    else if(forecasting_method == "ARIMA")
    {
      score_fore[ik,] = forecast(auto.arima(as.numeric(score[,ik])), h = fh)$mean
    }
  }

  # obtain forecasts in real-valued space

  fore_val = basis %*% score_fore

  # back-transformation

  f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, fh)
  for(ik in 1:fh)
  {
    f_x_t_star_fore[,ik] = exp(fore_val[,ik])/sum(exp(fore_val[,ik]))
    d_x_t_star_fore[,ik] = (f_x_t_star_fore[,ik] * alpha_x)/sum((f_x_t_star_fore[,ik] * alpha_x))
  }
  return(list(R2 = R2, recon = d_x_t_star_recon, fore_count = d_x_t_star_fore,
              alpha_x = alpha_x, h_x_t = h_x_t, basis_fore = basis, score_fore = score_fore,
              score = score, ncomp = ncomp))
}

#############
# An example
#############

## female

# weighted CoDa

SWE_female_pop_ex_weighted <- R_square_fit_ex(dat = SWE_female_pop[1:269,], weighting = "geom",
                                              geom_weight = 0.054,
                                              fh = 1, ncomp_selection = "eigen_ratio",
                                              forecasting_method = "ARIMA")

# unweighted (standard) CoDa

SWE_female_pop_ex_unweighted <- R_square_fit_ex(dat = SWE_female_pop[1:269,], weighting = "simple",
                                                fh = 1, ncomp_selection = "eigen_ratio",
                                                forecasting_method = "ARIMA")

# Graphical displays

age = 0:110
savepdf("Fig_4a", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_female_pop_ex_weighted$alpha_x * 10^5, ylim = c(0, 5500),
     type = "l", ylab = expression(alpha[x]), xlab = "",
     col = 4, lty = 4, main = "Swedish female data")
lines(age, SWE_female_pop_ex_unweighted$alpha_x * 10^5, col = 2, lty = 2)
dev.off()

savepdf("Fig_4c", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_female_pop_ex_unweighted$basis_fore, col = 2, lty = 2,
     type = "l", ylim = c(-0.2, 0.2), xlab = "", ylab = "First functional principal component")
lines(age, SWE_female_pop_ex_weighted$basis_fore, col = 4, lty = 4)
dev.off()

savepdf("Fig_4e", width = 12, height = 10, toplines = 0.8)
plot(fts(age, t(SWE_female_pop[1:269,])*10^5), col = "gray",
     xlab = "Age", ylab = "Life-table death count")
lines(age, SWE_female_pop[270,] * 10^5, col = 1, lty = 1)
lines(age, SWE_female_pop_ex_unweighted$fore * 10^5, col = 2, lty = 2)
lines(age, SWE_female_pop_ex_weighted$fore * 10^5, col = 4, lty = 4)
legend("topright", c("Observed data", "Holdout data", "Standard CoDa forecast", "Weighted CoDa forecast"),
       col = c("gray", 1, 2, 4), lty = c(1, 2, 4), cex = 0.8)
dev.off()

savepdf("Fig_4g", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_female_pop[270,] * 10^5, col = 1, lty = 1, type = "l", xlab = "Age", ylab = "Life-table death count")
lines(age, SWE_female_pop_ex_unweighted$fore * 10^5, col = 2, lty = 2)
lines(age, SWE_female_pop_ex_weighted$fore * 10^5, col = 4, lty = 4)
dev.off()

## male

# weighted CoDa

SWE_male_pop_ex_weighted <- R_square_fit_ex(dat = SWE_male_pop[1:269,], weighting = "geom",
                                            geom_weight = 0.079,
                                            fh = 1, ncomp_selection = "eigen_ratio",
                                            forecasting_method = "ARIMA")

# unweighted (standard) CoDa

SWE_male_pop_ex_unweighted <- R_square_fit_ex(dat = SWE_male_pop[1:269,], weighting = "simple",
                                              fh = 1, ncomp_selection = "eigen_ratio",
                                              forecasting_method = "ARIMA")

# Graphical displays

savepdf("Fig_4b", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_male_pop_ex_weighted$alpha_x * 10^5, type = "l", ylim = c(0, 6500),
     ylab = "", col  = 4, xlab = "",
     main = "Swedish male data", lty = 4)
lines(age, SWE_male_pop_ex_unweighted$alpha_x * 10^5, col = 2, lty = 2)
dev.off()

savepdf("Fig_4d", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_male_pop_ex_unweighted$basis_fore, col = 2, lty = 2,
     type = "l", ylim = c(-0.25, 0.2), xlab = "", ylab = "")
lines(age, SWE_male_pop_ex_weighted$basis_fore, col = 4, lty = 4)
dev.off()

savepdf("Fig_4f", width = 12, height = 10, toplines = 0.8)
plot(fts(age, t(SWE_male_pop[1:269,])*10^5), col = "gray",
     xlab = "Age", ylab = "")
lines(age, SWE_male_pop[270,] * 10^5, col = 1, lty = 1)
lines(age, SWE_male_pop_ex_unweighted$fore * 10^5, col = 2, lty = 2)
lines(age, SWE_male_pop_ex_weighted$fore * 10^5, col = 4, lty = 4)
dev.off()

savepdf("Fig_4h", width = 12, height = 10, toplines = 0.8)
plot(age, SWE_male_pop[270,] * 10^5, col = 1, lty = 1, type = "l", xlab = "Age", ylab = "Life-table death count")
lines(age, SWE_male_pop_ex_unweighted$fore * 10^5, col = 2, lty = 2)
lines(age, SWE_male_pop_ex_weighted$fore * 10^5, col = 4, lty = 4)
dev.off()
