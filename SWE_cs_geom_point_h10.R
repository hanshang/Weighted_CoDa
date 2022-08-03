#############
# CoDa model
#############

R_square_fit <- function(dat, weighting, geom_weight, fh,
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

    # log ratio

    h_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
      h_x_t[ik,] = log(dat[ik,]/alpha_x)
    }

    wq = diag(wei)
    h_x_t_weight = wq %*% h_x_t

    # static FPCA
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

#################
# select weights
#################

# training data: 1 to (n - 20)
# validation data: (n - 19) to (n - 10)
# testing data: (n - 9) to n

# lambda_para: parameter used in the geometrically decaying weights
# dat: (n by p) data matrix
# horizon: forecast horizon, chosen from 1 to 10
# criterion: error criterion = c("KL_div", "JS_div_simple", "JS_div_geo")

CoDa_geom_weight_select <- function(lambda_para, dat, horizon, order_selection, fore_method,
                                    criterion = c("KL_div", "JS_div_simple", "JS_div_geo"))
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), (11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        dum <- R_square_fit(dat = dat[1:(n - 21 + iw),], weighting = "geom",
                            geom_weight = lambda_para, fh = horizon,
                            ncomp_selection = order_selection, forecasting_method = fore_method)
        den_fore[,iw] = dum$fore_count[,horizon]
        rm(dum)
    }

    true_dat = replace(dat, which(dat == 10^-5), NA)

    err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        data_c = cbind(True = true_dat[(n - 21 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
        data_c = data_c[!is.na(data_c[,1]),]
        if(criterion == "KL_div")
        {
            err[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
        }
        else if(criterion == "JS_div_simple")
        {
            M = rowMeans(data_c)
            P_M = cbind(data_c[,1], M)
            E_M = cbind(as.numeric(data_c[,2]), M)
            colnames(E_M) = colnames(P_M) = c("True", "M")
            err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
            rm(M); rm(P_M); rm(E_M)
        }
        else if(criterion == "JS_div_geo")
        {
            M = apply(data_c, 1, geometric.mean)
            P_M = cbind(data_c[,1], M)
            E_M = cbind(as.numeric(data_c[,2]), M)
            colnames(E_M) = colnames(P_M) = c("True", "M")
            err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
            rm(M); rm(P_M); rm(E_M)
        }
        else
        {
            warning("Please specify a criterion from the list.")
        }
        print(ik); rm(ik)
    }
    rm(n); rm(true_dat); rm(den_fore)
    return(mean(err, na.rm = TRUE))
}

##########
## female
##########

## fore_method = "ARIMA"

# order_selection = "eigen_ratio"

CoDa_geo_weight_select_ARIMA_KL_div_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_geo_weight_select_ARIMA_KL_div_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_simple_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_geo_weight_select_ARIMA_JS_div_simple_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_geo_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_geo")
    CoDa_geo_weight_select_ARIMA_JS_div_geo_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

colnames(CoDa_geo_weight_select_ARIMA_KL_div_female) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_simple_female) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_geo_female) = c("Opt_para", "Opt_obj")

# order_selection = "fixed"

CoDa_geo_weight_select_ARIMA_KL_div_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_geo_weight_select_ARIMA_KL_div_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_simple_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_geo_weight_select_ARIMA_JS_div_simple_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_geo_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_geo")
    CoDa_geo_weight_select_ARIMA_JS_div_geo_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
colnames(CoDa_geo_weight_select_ARIMA_KL_div_female_ncomp_fixed) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_simple_female_ncomp_fixed) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_geo_female_ncomp_fixed) = c("opt_para", "opt_obj")

## fore_method = "RWD"

# order_selection = "eigen_ratio"

CoDa_geo_weight_select_RWD_KL_div_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_geo_weight_select_RWD_KL_div_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_simple_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_geo_weight_select_RWD_JS_div_simple_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_geo_female = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_geo")
    CoDa_geo_weight_select_RWD_JS_div_geo_female[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

colnames(CoDa_geo_weight_select_RWD_KL_div_female) =
colnames(CoDa_geo_weight_select_RWD_JS_div_simple_female) =
colnames(CoDa_geo_weight_select_RWD_JS_div_geo_female) = c("Opt_para", "Opt_obj")

# order_selection = "fixed"

CoDa_geo_weight_select_RWD_KL_div_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_geo_weight_select_RWD_KL_div_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_simple_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_geo_weight_select_RWD_JS_div_simple_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_geo_female_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_female_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_geo")
    CoDa_geo_weight_select_RWD_JS_div_geo_female_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
colnames(CoDa_geo_weight_select_RWD_KL_div_female_ncomp_fixed) =
colnames(CoDa_geo_weight_select_RWD_JS_div_simple_female_ncomp_fixed) =
colnames(CoDa_geo_weight_select_RWD_JS_div_geo_female_ncomp_fixed) = c("opt_para", "opt_obj")


########
## male
########

## fore_method = "ARIMA"

# order_selection = "eigen_ratio"

CoDa_geo_weight_select_ARIMA_KL_div_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_geo_weight_select_ARIMA_KL_div_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_simple_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_geo_weight_select_ARIMA_JS_div_simple_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_geo_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_geo")
    CoDa_geo_weight_select_ARIMA_JS_div_geo_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

colnames(CoDa_geo_weight_select_ARIMA_KL_div_male) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_simple_male) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_geo_male) = c("Opt_para", "Opt_obj")

# order_selection = "fixed"

CoDa_geo_weight_select_ARIMA_KL_div_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_geo_weight_select_ARIMA_KL_div_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_simple_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_geo_weight_select_ARIMA_JS_div_simple_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_ARIMA_JS_div_geo_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_geo")
    CoDa_geo_weight_select_ARIMA_JS_div_geo_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
colnames(CoDa_geo_weight_select_ARIMA_KL_div_male_ncomp_fixed) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_simple_male_ncomp_fixed) =
colnames(CoDa_geo_weight_select_ARIMA_JS_div_geo_male_ncomp_fixed) = c("opt_para", "opt_obj")

## fore_method = "RWD"

# order_selection = "eigenvalue ratio"

CoDa_geo_weight_select_RWD_KL_div_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_geo_weight_select_RWD_KL_div_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_simple_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_geo_weight_select_RWD_JS_div_simple_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_geo_male = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_geo")
    CoDa_geo_weight_select_RWD_JS_div_geo_male[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

colnames(CoDa_geo_weight_select_RWD_KL_div_male) =
colnames(CoDa_geo_weight_select_RWD_JS_div_simple_male) =
colnames(CoDa_geo_weight_select_RWD_JS_div_geo_male) = c("Opt_para", "Opt_obj")

# order_selection = "fixed"

CoDa_geo_weight_select_RWD_KL_div_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_geo_weight_select_RWD_KL_div_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_simple_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_geo_weight_select_RWD_JS_div_simple_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_geo_weight_select_RWD_JS_div_geo_male_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_geom_weight_select, method = "L-BFGS-B", lower = 0, upper = 1,
                 dat = SWE_male_pop, horizon = iwk, order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_geo")
    CoDa_geo_weight_select_RWD_JS_div_geo_male_ncomp_fixed[iwk,] <- c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
colnames(CoDa_geo_weight_select_RWD_KL_div_male_ncomp_fixed) =
colnames(CoDa_geo_weight_select_RWD_JS_div_simple_male_ncomp_fixed) =
colnames(CoDa_geo_weight_select_RWD_JS_div_geo_male_ncomp_fixed) = c("opt_para", "opt_obj")

###########################################
# tabulate the estimated weight parameters
###########################################

CoDa_geo_weight_select_para = cbind(CoDa_geo_weight_select_ARIMA_KL_div_female[,1],
                                    CoDa_geo_weight_select_ARIMA_JS_div_simple_female[,1],
                                    CoDa_geo_weight_select_ARIMA_JS_div_geo_female[,1],
                                    CoDa_geo_weight_select_RWD_KL_div_female[,1],
                                    CoDa_geo_weight_select_RWD_JS_div_simple_female[,1],
                                    CoDa_geo_weight_select_RWD_JS_div_geo_female[,1],

                                    CoDa_geo_weight_select_ARIMA_KL_div_male[,1],
                                    CoDa_geo_weight_select_ARIMA_JS_div_simple_male[,1],
                                    CoDa_geo_weight_select_ARIMA_JS_div_geo_male[,1],
                                    CoDa_geo_weight_select_RWD_KL_div_male[,1],
                                    CoDa_geo_weight_select_RWD_JS_div_simple_male[,1],
                                    CoDa_geo_weight_select_RWD_JS_div_geo_male[,1])


CoDa_geo_weight_select_para_ncomp_fixed = cbind(CoDa_geo_weight_select_ARIMA_KL_div_female_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_ARIMA_JS_div_simple_female_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_ARIMA_JS_div_geo_female_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_KL_div_female_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_JS_div_simple_female_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_JS_div_geo_female_ncomp_fixed[,1],

                                                CoDa_geo_weight_select_ARIMA_KL_div_male_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_ARIMA_JS_div_simple_male_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_ARIMA_JS_div_geo_male_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_KL_div_male_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_JS_div_simple_male_ncomp_fixed[,1],
                                                CoDa_geo_weight_select_RWD_JS_div_geo_male_ncomp_fixed[,1])


colnames(CoDa_geo_weight_select_para) = colnames(CoDa_geo_weight_select_para_ncomp_fixed) = rep(c("KLD", "JSD(simple)", "JSD(geo)"), 4)
rownames(CoDa_geo_weight_select_para) = rownames(CoDa_geo_weight_select_para_ncomp_fixed) = 1:10

#####################
# CoDa test function
#####################

CoDa_geom_test <- function(lambda_para_optim, dat, horizon, weit_choice, order_selection, fore_method,
                           criterion = c("KL_div", "JS_div_simple", "JS_div_geo"))
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), (11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        dum <- R_square_fit(dat = dat[1:(n - 11 + iw),], weighting = weit_choice, geom_weight = lambda_para_optim,
                            fh = horizon, ncomp_selection = order_selection, forecasting_method = fore_method)
        den_fore[,iw] = dum$fore_count[,horizon]
    }

    true_dat = replace(dat, which(dat == 10^-5), NA)

    err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        data_c = cbind(True = true_dat[(n - 11 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
        data_c = data_c[!is.na(data_c[,1]),]
        if(criterion == "KL_div")
        {
            err[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
        }
        else if(criterion == "JS_div_simple")
        {
            M = rowMeans(data_c)
            P_M = cbind(data_c[,1], M)
            E_M = cbind(as.numeric(data_c[,2]), M)
            colnames(E_M) = colnames(P_M) = c("True", "M")
            err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
            rm(M); rm(P_M); rm(E_M)
        }
        else if(criterion == "JS_div_geo")
        {
            M = apply(data_c, 1, geometric.mean)
            P_M = cbind(data_c[,1], M)
            E_M = cbind(as.numeric(data_c[,2]), M)
            colnames(E_M) = colnames(P_M) = c("True", "M")
            err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
            rm(M); rm(P_M); rm(E_M)
        }
        else
        {
            warning("Please specify a criterion from the list.")
        }
        rm(data_c)
    }
    rm(den_fore)
    return(mean(err, na.rm = TRUE))
}

##################
# female (simple)
##################

## ncomp_selection = "eigen_ratio"

# fore_method = "ARIMA"

CoDa_simple_test_ARIMA_KL_div_female =
CoDa_simple_test_ARIMA_JS_div_simple_female = CoDa_simple_test_ARIMA_JS_div_geo_female = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_ARIMA_KL_div_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_simple_test_ARIMA_JS_div_simple_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_simple_test_ARIMA_JS_div_geo_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_ARIMA_female = cbind(CoDa_simple_test_ARIMA_KL_div_female,
                                      CoDa_simple_test_ARIMA_JS_div_simple_female, CoDa_simple_test_ARIMA_JS_div_geo_female)
colnames(CoDa_simple_test_ARIMA_female) = c("KL_div", "JS_div_simple", "JS_div_geo")

# fore_method = "RWD"

CoDa_simple_test_RWD_KL_div_female =
CoDa_simple_test_RWD_JS_div_simple_female = CoDa_simple_test_RWD_JS_div_geo_female = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_RWD_KL_div_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_simple_test_RWD_JS_div_simple_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_simple_test_RWD_JS_div_geo_female[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_RWD_female = cbind(CoDa_simple_test_RWD_KL_div_female,
                                    CoDa_simple_test_RWD_JS_div_simple_female, CoDa_simple_test_RWD_JS_div_geo_female)
colnames(CoDa_simple_test_RWD_female) = c("KL_div", "JS_div_simple", "JS_div_geo")

## ncomp_selection = "fixed"

# fore_method = "ARIMA"

CoDa_simple_test_ARIMA_KL_div_female_ncomp_fixed =
CoDa_simple_test_ARIMA_JS_div_simple_female_ncomp_fixed = CoDa_simple_test_ARIMA_JS_div_geo_female_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_ARIMA_KL_div_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_simple_test_ARIMA_JS_div_simple_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_simple_test_ARIMA_JS_div_geo_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_ARIMA_female_ncomp_fixed_results = cbind(CoDa_simple_test_ARIMA_KL_div_female_ncomp_fixed,
                                                          CoDa_simple_test_ARIMA_JS_div_simple_female_ncomp_fixed,
                                                          CoDa_simple_test_ARIMA_JS_div_geo_female_ncomp_fixed)
colnames(CoDa_simple_test_ARIMA_female_ncomp_fixed_results) = c("KL div", "JS div (simple)", "JS div (geo)")

# fore_method = "RWD"

CoDa_simple_test_RWD_KL_div_female_ncomp_fixed =
CoDa_simple_test_RWD_JS_div_simple_female_ncomp_fixed = CoDa_simple_test_RWD_JS_div_geo_female_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_RWD_KL_div_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_simple_test_RWD_JS_div_simple_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_simple_test_RWD_JS_div_geo_female_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_female_pop, horizon = iwk, weit_choice = "simple", order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_RWD_female_ncomp_fixed_results = cbind(CoDa_simple_test_RWD_KL_div_female_ncomp_fixed,
                                                        CoDa_simple_test_RWD_JS_div_simple_female_ncomp_fixed,
                                                        CoDa_simple_test_RWD_JS_div_geo_female_ncomp_fixed)
colnames(CoDa_simple_test_RWD_female_ncomp_fixed_results) = c("KL div", "JS div (simple)", "JS div (geo)")


###############
# female (geo)
###############

## fore_method = "ARIMA"

# order_selection = "eigenvalue ratio"

CoDa_geo_test_ARIMA_KL_div_female =
CoDa_geo_test_ARIMA_JS_div_simple_female = CoDa_geo_test_ARIMA_JS_div_geo_female = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_ARIMA_KL_div_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_KL_div_female[iwk,1],
                                                            dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                            fore_method = "ARIMA", criterion = "KL_div")

    CoDa_geo_test_ARIMA_JS_div_simple_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_simple_female[iwk,1],
                                                                   dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                                   fore_method = "ARIMA", criterion = "JS_div_simple")

    CoDa_geo_test_ARIMA_JS_div_geo_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_geo_female[iwk,1],
                                                                dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                                fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_ARIMA_female = cbind(CoDa_geo_test_ARIMA_KL_div_female,
                                   CoDa_geo_test_ARIMA_JS_div_simple_female, CoDa_geo_test_ARIMA_JS_div_geo_female)
colnames(CoDa_geo_test_ARIMA_female) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_geo_test_ARIMA_KL_div_female_ncomp_fixed =
CoDa_geo_test_ARIMA_JS_div_simple_female_ncomp_fixed = CoDa_geo_test_ARIMA_JS_div_geo_female_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_ARIMA_KL_div_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_KL_div_female_ncomp_fixed[iwk,1],
                                                                        dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                        fore_method = "ARIMA", criterion = "KL_div")

    CoDa_geo_test_ARIMA_JS_div_simple_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_simple_female_ncomp_fixed[iwk,1],
                                                                               dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                               fore_method = "ARIMA", criterion = "JS_div_simple")

    CoDa_geo_test_ARIMA_JS_div_geo_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_geo_female_ncomp_fixed[iwk,1],
                                                                            dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                            fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_ARIMA_female_ncomp_fixed_results = cbind(CoDa_geo_test_ARIMA_KL_div_female_ncomp_fixed,
                                                       CoDa_geo_test_ARIMA_JS_div_simple_female_ncomp_fixed,
                                                       CoDa_geo_test_ARIMA_JS_div_geo_female_ncomp_fixed)
colnames(CoDa_geo_test_ARIMA_female_ncomp_fixed_results) = c("KL div", "JS div (simple)", "JS div (geo)")

## fore_method = "RWD"

# order_selection = "eigenvalue ratio"

CoDa_geo_test_RWD_KL_div_female =
CoDa_geo_test_RWD_JS_div_simple_female = CoDa_geo_test_RWD_JS_div_geo_female = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_RWD_KL_div_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_KL_div_female[iwk,1],
                                                          dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                          fore_method = "RWF_drift", criterion = "KL_div")

    CoDa_geo_test_RWD_JS_div_simple_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_simple_female[iwk,1],
                                                                 dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                                 fore_method = "RWF_drift", criterion = "JS_div_simple")

    CoDa_geo_test_RWD_JS_div_geo_female[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_geo_female[iwk,1],
                                                              dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                              fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_RWD_female = cbind(CoDa_geo_test_RWD_KL_div_female,
                                 CoDa_geo_test_RWD_JS_div_simple_female, CoDa_geo_test_RWD_JS_div_geo_female)
colnames(CoDa_geo_test_RWD_female) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_geo_test_RWD_KL_div_female_ncomp_fixed =
CoDa_geo_test_RWD_JS_div_simple_female_ncomp_fixed = CoDa_geo_test_RWD_JS_div_geo_female_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_RWD_KL_div_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_KL_div_female_ncomp_fixed[iwk,1],
                                                                      dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                      fore_method = "RWF_drift", criterion = "KL_div")

    CoDa_geo_test_RWD_JS_div_simple_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_simple_female_ncomp_fixed[iwk,1],
                                                                             dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                             fore_method = "RWF_drift", criterion = "JS_div_simple")

    CoDa_geo_test_RWD_JS_div_geo_female_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_geo_female_ncomp_fixed[iwk,1],
                                                                          dat = SWE_female_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                          fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}
CoDa_geo_test_RWD_female_ncomp_fixed = cbind(CoDa_geo_test_RWD_KL_div_female_ncomp_fixed,
                                             CoDa_geo_test_RWD_JS_div_simple_female_ncomp_fixed, CoDa_geo_test_RWD_JS_div_geo_female_ncomp_fixed)
colnames(CoDa_geo_test_RWD_female_ncomp_fixed) = c("KL_div", "JS_div_simple", "JS_div_geo")


################
# male (simple)
################

## fore_method = "ARIMA"

# order_selection = "eigenvalue ratio"

CoDa_simple_test_ARIMA_KL_div_male =
CoDa_simple_test_ARIMA_JS_div_simple_male = CoDa_simple_test_ARIMA_JS_div_geo_male = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_ARIMA_KL_div_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                             order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_simple_test_ARIMA_JS_div_simple_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                    order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_simple_test_ARIMA_JS_div_geo_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                 order_selection = "eigen_ratio", fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_ARIMA_male = cbind(CoDa_simple_test_ARIMA_KL_div_male,
                                    CoDa_simple_test_ARIMA_JS_div_simple_male, CoDa_simple_test_ARIMA_JS_div_geo_male)
colnames(CoDa_simple_test_ARIMA_male) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_simple_test_ARIMA_KL_div_male_ncomp_fixed =
CoDa_simple_test_ARIMA_JS_div_simple_male_ncomp_fixed = CoDa_simple_test_ARIMA_JS_div_geo_male_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_ARIMA_KL_div_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                         order_selection = "fixed", fore_method = "ARIMA", criterion = "KL_div")
    CoDa_simple_test_ARIMA_JS_div_simple_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                                order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_simple")
    CoDa_simple_test_ARIMA_JS_div_geo_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                             order_selection = "fixed", fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}
CoDa_simple_test_ARIMA_male_ncomp_fixed_results = cbind(CoDa_simple_test_ARIMA_KL_div_male_ncomp_fixed,
                                                        CoDa_simple_test_ARIMA_JS_div_simple_male_ncomp_fixed, CoDa_simple_test_ARIMA_JS_div_geo_male_ncomp_fixed)
colnames(CoDa_simple_test_ARIMA_male_ncomp_fixed_results) = c("KL_div", "JS_div_simple", "JS_div_geo")


## fore_method = "RWD"

# order_selection = "eigenvalue ratio"

CoDa_simple_test_RWD_KL_div_male =
CoDa_simple_test_RWD_JS_div_simple_male = CoDa_simple_test_RWD_JS_div_geo_male = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_RWD_KL_div_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                           order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_simple_test_RWD_JS_div_simple_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                  order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_simple_test_RWD_JS_div_geo_male[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                               order_selection = "eigen_ratio", fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_simple_test_RWD_male = cbind(CoDa_simple_test_RWD_KL_div_male,
                                  CoDa_simple_test_RWD_JS_div_simple_male, CoDa_simple_test_RWD_JS_div_geo_male)
colnames(CoDa_simple_test_RWD_male) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_simple_test_RWD_KL_div_male_ncomp_fixed =
CoDa_simple_test_RWD_JS_div_simple_male_ncomp_fixed = CoDa_simple_test_RWD_JS_div_geo_male_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_simple_test_RWD_KL_div_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                       order_selection = "fixed", fore_method = "RWF_drift", criterion = "KL_div")
    CoDa_simple_test_RWD_JS_div_simple_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                              order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_simple")
    CoDa_simple_test_RWD_JS_div_geo_male_ncomp_fixed[iwk] = CoDa_geom_test(dat = SWE_male_pop, horizon = iwk, weit_choice = "simple",
                                                                           order_selection = "fixed", fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}
CoDa_simple_test_RWD_male_ncomp_fixed = cbind(CoDa_simple_test_RWD_KL_div_male_ncomp_fixed,
                                              CoDa_simple_test_RWD_JS_div_simple_male_ncomp_fixed, CoDa_simple_test_RWD_JS_div_geo_male_ncomp_fixed)
colnames(CoDa_simple_test_RWD_male_ncomp_fixed) = c("KL_div", "JS_div_simple", "JS_div_geo")


#############
# male (geo)
#############

## fore_method = "ARIMA"

# order_selection = "eigenvalue ratio"

CoDa_geo_test_ARIMA_KL_div_male =
CoDa_geo_test_ARIMA_JS_div_simple_male = CoDa_geo_test_ARIMA_JS_div_geo_male = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_ARIMA_KL_div_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_KL_div_male[iwk,1],
                                                          dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                          fore_method = "ARIMA", criterion = "KL_div")

    CoDa_geo_test_ARIMA_JS_div_simple_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_simple_male[iwk,1],
                                                                 dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                                 fore_method = "ARIMA", criterion = "JS_div_simple")

    CoDa_geo_test_ARIMA_JS_div_geo_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_geo_male[iwk,1],
                                                              dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                              fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_ARIMA_male = cbind(CoDa_geo_test_ARIMA_KL_div_male,
                                 CoDa_geo_test_ARIMA_JS_div_simple_male, CoDa_geo_test_ARIMA_JS_div_geo_male)
colnames(CoDa_geo_test_ARIMA_male) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_geo_test_ARIMA_KL_div_male_ncomp_fixed =
CoDa_geo_test_ARIMA_JS_div_simple_male_ncomp_fixed = CoDa_geo_test_ARIMA_JS_div_geo_male_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_ARIMA_KL_div_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_KL_div_male_ncomp_fixed[iwk,1],
                                                                      dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                      fore_method = "ARIMA", criterion = "KL_div")

    CoDa_geo_test_ARIMA_JS_div_simple_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_simple_male_ncomp_fixed[iwk,1],
                                                                             dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                             fore_method = "ARIMA", criterion = "JS_div_simple")

    CoDa_geo_test_ARIMA_JS_div_geo_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_ARIMA_JS_div_geo_male_ncomp_fixed[iwk,1],
                                                                          dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                          fore_method = "ARIMA", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}
CoDa_geo_test_ARIMA_male_ncomp_fixed_results = cbind(CoDa_geo_test_ARIMA_KL_div_male_ncomp_fixed,
                                                     CoDa_geo_test_ARIMA_JS_div_simple_male_ncomp_fixed, CoDa_geo_test_ARIMA_JS_div_geo_male_ncomp_fixed)
colnames(CoDa_geo_test_ARIMA_male_ncomp_fixed_results) = c("KL div", "JS div (simple)", "JS div (geo)")


## fore_method = "RWD"

# order_selection = "eigenvalue ratio"

CoDa_geo_test_RWD_KL_div_male =
CoDa_geo_test_RWD_JS_div_simple_male = CoDa_geo_test_RWD_JS_div_geo_male = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_RWD_KL_div_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_KL_div_male[iwk,1],
                                                        dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                        fore_method = "RWF_drift", criterion = "KL_div")

    CoDa_geo_test_RWD_JS_div_simple_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_simple_male[iwk,1],
                                                               dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                               fore_method = "RWF_drift", criterion = "JS_div_simple")

    CoDa_geo_test_RWD_JS_div_geo_male[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_geo_male[iwk,1],
                                                            dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "eigen_ratio",
                                                            fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_RWD_male = cbind(CoDa_geo_test_RWD_KL_div_male,
                               CoDa_geo_test_RWD_JS_div_simple_male, CoDa_geo_test_RWD_JS_div_geo_male)
colnames(CoDa_geo_test_RWD_male) = c("KL_div", "JS_div_simple", "JS_div_geo")

# order_selection = "fixed"

CoDa_geo_test_RWD_KL_div_male_ncomp_fixed =
CoDa_geo_test_RWD_JS_div_simple_male_ncomp_fixed = CoDa_geo_test_RWD_JS_div_geo_male_ncomp_fixed = vector("numeric", 10)
for(iwk in 1:10)
{
    CoDa_geo_test_RWD_KL_div_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_KL_div_male_ncomp_fixed[iwk,1],
                                                                    dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                    fore_method = "RWF_drift", criterion = "KL_div")

    CoDa_geo_test_RWD_JS_div_simple_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_simple_male_ncomp_fixed[iwk,1],
                                                                           dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                           fore_method = "RWF_drift", criterion = "JS_div_simple")

    CoDa_geo_test_RWD_JS_div_geo_male_ncomp_fixed[iwk] = CoDa_geom_test(lambda_para_optim = CoDa_geo_weight_select_RWD_JS_div_geo_male_ncomp_fixed[iwk,1],
                                                                        dat = SWE_male_pop, horizon = iwk,  weit_choice = "geom", order_selection = "fixed",
                                                                        fore_method = "RWF_drift", criterion = "JS_div_geo")
    print(iwk); rm(iwk)
}

CoDa_geo_test_RWD_male_ncomp_fixed = cbind(CoDa_geo_test_RWD_KL_div_male_ncomp_fixed,
                                           CoDa_geo_test_RWD_JS_div_simple_male_ncomp_fixed,
                                           CoDa_geo_test_RWD_JS_div_geo_male_ncomp_fixed)
colnames(CoDa_geo_test_RWD_male_ncomp_fixed) = c("KL_div", "JS_div_simple", "JS_div_geo")

########################################
# tabulate results via xtable R package
########################################

CoDa_simple_geo_ARIMA_female = cbind(CoDa_simple_test_ARIMA_female, CoDa_geo_test_ARIMA_female,
                                     CoDa_simple_test_ARIMA_female_ncomp_fixed_results, CoDa_geo_test_ARIMA_female_ncomp_fixed_results) * 100

CoDa_simple_geo_RWD_female = cbind(CoDa_simple_test_RWD_female, CoDa_geo_test_RWD_female,
                                   CoDa_simple_test_RWD_female_ncomp_fixed_results, CoDa_geo_test_RWD_female_ncomp_fixed) * 100

CoDa_simple_geo_ARIMA_male = cbind(CoDa_simple_test_ARIMA_male, CoDa_geo_test_ARIMA_male,
                                   CoDa_simple_test_ARIMA_male_ncomp_fixed_results, CoDa_geo_test_ARIMA_male_ncomp_fixed_results) * 100

CoDa_simple_geo_RWD_male = cbind(CoDa_simple_test_RWD_male, CoDa_geo_test_RWD_male,
                                 CoDa_simple_test_RWD_male_ncomp_fixed, CoDa_geo_test_RWD_male_ncomp_fixed) * 100

colnames(CoDa_simple_geo_ARIMA_female) = colnames(CoDa_simple_geo_RWD_female) =
colnames(CoDa_simple_geo_ARIMA_male) = colnames(CoDa_simple_geo_RWD_male) = rep(c("KLD", "JSD(simple)", "JSD(geo)"), 4)

require(xtable)
xtable(rbind(CoDa_simple_geo_ARIMA_female, colMeans(CoDa_simple_geo_ARIMA_female)), digits = 3)
xtable(rbind(CoDa_simple_geo_RWD_female, colMeans(CoDa_simple_geo_RWD_female)), digits = 3)

xtable(rbind(CoDa_simple_geo_ARIMA_male, colMeans(CoDa_simple_geo_ARIMA_male)), digits = 3)
xtable(rbind(CoDa_simple_geo_RWD_male, colMeans(CoDa_simple_geo_RWD_male)), digits = 3)
