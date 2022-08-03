####################################
# CoDa with nonparametric bootstrap
####################################

# dat: first data set (n by p)
# weighting: simple or geom weighting
# geom_weight: geometrically decaying weights
# fh: forecast horizon
# B: number of bootstrap samples
# fmethod: forecasting method

# set a working directory

setwd("~/Dropbox/Todos/cs_geometric/code")

# load R packages

source("load_packages.R")

CoDa_interval_bootstrap <- function(dat, weighting = c("simple", "geom"), geom_weight,
                                    fh, ncomp_selection = c("eigen_ratio", "fixed"),
                                    B = 1000, fmethod = c("ARIMA", "RWF_drift"))
{
    n_year = nrow(dat)
    n_age = ncol(dat)
    year_index = rownames(dat)
    age_index = colnames(dat)

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
    else
    {
      warning("Weighting must be either simple or geometrically decaying.")
    }

    # compute geometric mean

    geom_mean = vector("numeric", n_age)
    geom_mean = rep(0, n_age)
    for(ik in 1:n_year)
    {
      geom_mean = geom_mean + wei[ik] * log(dat[ik,])
    }
    alpha_x = exp(geom_mean)

    # standardization

    h_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
      h_x_t[ik,] = log(dat[ik,]/alpha_x)
    }

    wq = diag(wei)
    h_x_t_weight = wq %*% h_x_t

    # perform an eigen-decomposition to beta_t_u

    eigen_beta = eigen(cov(h_x_t_weight))

    # select ncomp by LRS(JASA, 2020) eigen-value ratio tests

    lambda_val = eigen_beta$values[which(eigen_beta$values > 0)]
    if(ncomp_selection == "eigen_ratio")
    {
      k_max = length(which(lambda_val >= mean(lambda_val)))

      tau = 1/log(max(lambda_val[1], length(lambda_val)))
      eigen_val_ratio = vector("numeric", k_max)
      for(ik in 1:k_max)
      {
        eigen_val_ratio[ik] = lambda_val[ik+1]/lambda_val[ik] * ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 0) + ifelse(lambda_val[ik]/lambda_val[1] < tau, 1, 0)
      }
      ncomp = max(1, which.min(eigen_val_ratio))
      if(exists("ncomp") == "FALSE")
      {
        ncomp = 1
      }
    }
    else if(ncomp_selection == "fixed")
    {
      ncomp = 6
    }
    else
    {
      warning("ncomp selection method is not defined.")
    }

    basis = matrix(eigen_beta$vectors[,1:ncomp], ncol = ncomp)
    scores = h_x_t %*% basis
    recon = basis %*% t(scores)

    resi = t(h_x_t) - recon
    colnames(recon) = rownames(scores) = year_index

    # forecast h-step-ahead PC scores

    olivia = matrix(NA, ncomp, fh)
    if(fmethod == "ARIMA")
    {
      for(ik in 1:ncomp)
      {
        olivia[ik,] = forecast(auto.arima(ts(scores[,ik], start = head(year_index, 1), end = tail(year_index, 1))),
                               h = fh)$mean
      }
    }
    else if(fmethod == "RWF_drift")
    {
      for(ik in 1:ncomp)
      {
        olivia[ik,] = rwf(ts(scores[,ik], start = head(year_index, 1), end = tail(year_index, 1)),
                          h = fh, drift = TRUE)$mean
      }
    }
    else
    {
      warning("fmethod is not listed.")
    }

    # forecast errors based on observed scores

    forerr = matrix(NA, (n_year - ncomp - fh), ncomp)
    for(ik in (fh + 1):(n_year - ncomp))
    {
        k = ik + (ncomp - fh)
        fore = matrix(NA, 1, ncomp)
        if(fmethod == "ARIMA")
        {
            for(j in 1:ncomp)
            {
                fore[,j] = forecast(auto.arima(scores[1:k,j]), h = fh)$mean[fh]
            }
        }
        else if(fmethod == "RWF_drift")
        {
            for(j in 1:ncomp)
            {
                fore[,j] = rwf(scores[1:k,j], h = fh, drift = TRUE)$mean[fh]
            }
        }
        forerr[ik - fh,] = scores[k + fh,] - fore
    }

    # bootstrapping residuals

    q = array(NA, dim = c(n_age, B, fh))
    for(j in 1:fh)
    {
        q[,,j] = resi[,sample(1:n_year, size = B, replace = TRUE)]
    }

    # bootstrapping PC score error

    ny = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
        for(i in 1:ncomp)
        {
            ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE)
        }
    }

    # adding PC score error to PC score forecasts

    oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
    fo = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
        for(i in 1:B)
        {
            fo[,i,j] = oli[,i,j] + ny[,i,j]
        }
    }

    # loading * score forecasts + residuals

    pred = array(NA, dim = c(n_age, B, fh))
    for(j in 1:fh)
    {
        for(i in 1:B)
        {
            pred[,i,j] = basis %*% fo[,i,j] + q[,i,j]
        }
    }

    beta_t_u_boot = pred[,,fh]

    # inverse log-ratio transformation

    f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, B)
    for(ik in 1:B)
    {
        f_x_t_star_fore[,ik] = exp(beta_t_u_boot[,ik])/sum(exp(beta_t_u_boot[,ik]))
        d_x_t_star_fore[,ik] = (f_x_t_star_fore[,ik] * alpha_x)/sum((f_x_t_star_fore[,ik] * alpha_x))
    }
    colnames(d_x_t_star_fore) = 1:B
    rownames(d_x_t_star_fore) = age_index
    PI_80 = apply(d_x_t_star_fore, 1, quantile, c(0.1, 0.9))
    PI_95 = apply(d_x_t_star_fore, 1, quantile, c(0.025, 0.975))
    return(list(PI_80 = PI_80, PI_95 = PI_95))
}

#####################
# CoDa test function
#####################

# interval score

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(mean(score), cover, cpd))
}

####################################################################################################
# finding the optimal geometrically decaying weights for selecting lambdas in the validation sample
####################################################################################################

CoDa_int_geo_weight <- function(lambda_para_optim, dat, weit_choice, horizon,
                                order_selection, fore_method,
                                criterion = c("score_80", "CPD_80", "score_95", "CPD_95"))
{
    n = nrow(dat)

    # computing lower and upper bounds

    train_den_int_80 = train_den_int_95 = array(NA, dim = c(2, ncol(dat), (11 - horizon)))
    for(iw in 1:(11 - horizon))
    {
      dum <- CoDa_interval_bootstrap(dat = dat[1:(n - 21 + iw),], weighting = weit_choice,
                                     geom_weight = lambda_para_optim,
                                     fh = horizon, ncomp_selection = order_selection,
                                     fmethod = fore_method)
      train_den_int_80[,,iw] = dum$PI_80
      train_den_int_95[,,iw] = dum$PI_95
      rm(dum); rm(iw)
    }

    # computing interval scores

    interval_score_80 <- interval_score(holdout = t(dat[(n - 20 + horizon):(n - 10),]),
                                        lb = train_den_int_80[1,,], ub = train_den_int_80[2,,],
                                        alpha = 0.2)
    interval_score_95 <- interval_score(holdout = t(dat[(n - 20 + horizon):(n - 10),]),
                                        lb = train_den_int_95[1,,], ub = train_den_int_95[2,,],
                                        alpha = 0.05)

    if(criterion == "score_80")
    {
      return(interval_score_80[1])
    }
    else if(criterion == "CPD_80")
    {
      return(interval_score_80[3])
    }
    else if(criterion == "score_95")
    {
      return(interval_score_95[1])
    }
    else if(criterion == "CPD_95")
    {
      return(interval_score_95[3])
    }
    else
    {
      warning("criterion is not listed.")
    }
}

#################################
## fore_method = "RWF_drift"
# criterion = "eigenvalue ratio"
#################################

##########################
## criterion = "score_80"
##########################

# female

CoDa_int_geo_weight_RWD_female_para_score_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "score_80")
    CoDa_int_geo_weight_RWD_female_para_score_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_female_para_score_80) = 1:10
colnames(CoDa_int_geo_weight_RWD_female_para_score_80) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_RWD_female_para_score_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "score_80")
    CoDa_int_geo_weight_RWD_female_para_score_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_female_para_score_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_female_para_score_80_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_RWD_male_para_score_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom",horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "score_80")
    CoDa_int_geo_weight_RWD_male_para_score_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_male_para_score_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom",horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "score_80")
    CoDa_int_geo_weight_RWD_male_para_score_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_male_para_score_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_male_para_score_80_ncomp_fixed) = c("Opt_para", "Opt_score_80")

##########################
## criterion = "score_95"
##########################

# female

CoDa_int_geo_weight_RWD_female_para_score_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "score_95")
    CoDa_int_geo_weight_RWD_female_para_score_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_female_para_score_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "score_95")
    CoDa_int_geo_weight_RWD_female_para_score_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_female_para_score_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_female_para_score_95_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_RWD_male_para_score_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "score_95")
    CoDa_int_geo_weight_RWD_male_para_score_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_male_para_score_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "score_95")
    CoDa_int_geo_weight_RWD_male_para_score_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_male_para_score_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_male_para_score_95_ncomp_fixed) = c("opt_para", "opt_obj")


#######################
# criterion = "CPD_80"
#######################

# female

CoDa_int_geo_weight_RWD_female_para_CPD_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "CPD_80")
    CoDa_int_geo_weight_RWD_female_para_CPD_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "CPD_80")
    CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_RWD_male_para_CPD_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "CPD_80")
    CoDa_int_geo_weight_RWD_male_para_CPD_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}


CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "CPD_80")
    CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed) = c("opt_para", "opt_obj")


#######################
# criterion = "CPD_95"
#######################

# female

CoDa_int_geo_weight_RWD_female_para_CPD_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "CPD_95")
    CoDa_int_geo_weight_RWD_female_para_CPD_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "CPD_95")
    CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_RWD_male_para_CPD_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "RWF_drift", criterion = "CPD_95")
    CoDa_int_geo_weight_RWD_male_para_CPD_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}

CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "RWF_drift", criterion = "CPD_95")
    CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed) = c("Opt_para", "Opt_obj")

############################
## fore_method = "RW_drift"
## testing data
############################

## order_selection = "eigen_ratio"

# female

array_int_score_80_RWD = array_int_score_95_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "RWF_drift")
    array_int_score_80_RWD[iwk,] = dum$interval_score_80
    array_int_score_95_RWD[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_RWD) = rownames(array_int_score_95_RWD) = 1:10
colnames(array_int_score_80_RWD) = colnames(array_int_score_95_RWD) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_RWD), 4) # 489.3312 0.9555 0.1555
round(colMeans(array_int_score_95_RWD), 4) # 848.9549 0.9998 0.0498

# male

array_int_score_80_male_RWD = array_int_score_95_male_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "RWF_drift")
    array_int_score_80_male_RWD[iwk,] = dum$interval_score_80
    array_int_score_95_male_RWD[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_male_RWD) = rownames(array_int_score_95_male_RWD) = 1:10
colnames(array_int_score_80_male_RWD) = colnames(array_int_score_95_male_RWD) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_male_RWD), 4) # 542.3092 0.9876 0.1876
round(colMeans(array_int_score_95_male_RWD), 4) # 978.2880 0.9999 0.0499

## order_selection = "fixed"

# female

array_int_score_80_RWD_ncomp_fixed = array_int_score_95_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "RWF_drift")
    array_int_score_80_RWD_ncomp_fixed[iwk,] = dum$interval_score_80
    array_int_score_95_RWD_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_RWD_ncomp_fixed) = rownames(array_int_score_95_RWD_ncomp_fixed) = 1:10
colnames(array_int_score_80_RWD_ncomp_fixed) = colnames(array_int_score_95_RWD_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_RWD_ncomp_fixed), 4) # 462.1234 0.9866 0.1866
round(colMeans(array_int_score_95_RWD_ncomp_fixed), 4) # 829.6810 0.9995 0.0495

# male

array_int_score_80_male_RWD_ncomp_fixed = array_int_score_95_male_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "RWF_drift")
    array_int_score_80_male_RWD_ncomp_fixed[iwk,] = dum$interval_score_80
    array_int_score_95_male_RWD_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_male_RWD_ncomp_fixed) = rownames(array_int_score_95_male_RWD_ncomp_fixed) = 1:10
colnames(array_int_score_80_male_RWD_ncomp_fixed) = colnames(array_int_score_95_male_RWD_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_male_RWD_ncomp_fixed), 4) # 516.5812 0.9942 0.1942
round(colMeans(array_int_score_95_male_RWD_ncomp_fixed), 4) # 939.6206 0.9998 0.0498


#################################
# geometrically decaying weights
#################################

### order_selection = "eigen_ratio"

## female

# criterion = "CPD_80"

geom_female_array_int_CPD_80_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_female_para_CPD_80[iwk,1],
                         fore_method = "RWF_drift")
    geom_female_array_int_CPD_80_RWD[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_female_array_int_CPD_95_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_female_para_CPD_95[iwk,1],
                         fore_method = "RWF_drift")
    geom_female_array_int_CPD_95_RWD[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_female_array_int_CPD_80_RWD) = rownames(geom_female_array_int_CPD_95_RWD) = 1:10
colnames(geom_female_array_int_CPD_80_RWD) = colnames(geom_female_array_int_CPD_95_RWD) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_female_array_int_CPD_80_RWD), 4) #  893.9917 0.8044 0.0585
round(colMeans(geom_female_array_int_CPD_95_RWD), 4) # 1593.7958 0.9117 0.0426

## male

# criterion = "CPD_80"

geom_male_array_int_CPD_80_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_male_para_CPD_80[iwk,1],
                         fore_method = "RWF_drift")
    geom_male_array_int_CPD_80_RWD[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_male_array_int_CPD_95_RWD = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_male_para_CPD_95[iwk, 1],
                         fore_method = "RWF_drift")
    geom_male_array_int_CPD_95_RWD[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_male_array_int_CPD_80_RWD) = rownames(geom_male_array_int_CPD_95_RWD) = 1:10
colnames(geom_male_array_int_CPD_80_RWD) = colnames(geom_male_array_int_CPD_95_RWD) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_male_array_int_CPD_80_RWD), 4) #  968.7310 0.7707 0.0356
round(colMeans(geom_male_array_int_CPD_95_RWD), 4) # 1515.3710 0.9467 0.0059

### order_selection = "fixed"

## female

# criterion = "CPD_80"

geom_female_array_int_CPD_80_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed[iwk,1],
                         fore_method = "RWF_drift")
    geom_female_array_int_CPD_80_RWD_ncomp_fixed[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_female_array_int_CPD_95_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed[iwk,1],
                         fore_method = "RWF_drift")
    geom_female_array_int_CPD_95_RWD_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_female_array_int_CPD_80_RWD_ncomp_fixed) = rownames(geom_female_array_int_CPD_95_RWD_ncomp_fixed) = 1:10
colnames(geom_female_array_int_CPD_80_RWD_ncomp_fixed) = colnames(geom_female_array_int_CPD_95_RWD_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_female_array_int_CPD_80_RWD_ncomp_fixed), 4) #  847.0547 0.7656 0.0445
round(colMeans(geom_female_array_int_CPD_95_RWD_ncomp_fixed), 4) # 1239.8843 0.9369 0.0235

## male

# criterion = "CPD_80"

geom_male_array_int_CPD_80_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed[iwk,1],
                         fore_method = "RWF_drift")
    geom_male_array_int_CPD_80_RWD_ncomp_fixed[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_male_array_int_CPD_95_RWD_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed[iwk,1],
                         fore_method = "RWF_drift")
    geom_male_array_int_CPD_95_RWD_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_male_array_int_CPD_80_RWD_ncomp_fixed) = rownames(geom_male_array_int_CPD_95_RWD_ncomp_fixed) = 1:10
colnames(geom_male_array_int_CPD_80_RWD_ncomp_fixed) = colnames(geom_male_array_int_CPD_95_RWD_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_male_array_int_CPD_80_RWD_ncomp_fixed), 4) #  993.7918 0.7475 0.0613
round(colMeans(geom_male_array_int_CPD_95_RWD_ncomp_fixed), 4) # 1473.0967 0.9423 0.0123

##########
# summary
##########

geom_RWD_array_int_CPD = cbind(array_int_score_80_RWD[,3], array_int_score_95_RWD[,3],
                                 geom_female_array_int_CPD_80_RWD[,3], geom_female_array_int_CPD_95_RWD[,3],
                                 array_int_score_80_RWD_ncomp_fixed[,3], array_int_score_95_RWD_ncomp_fixed[,3],
                                 geom_female_array_int_CPD_80_RWD_ncomp_fixed[,3], geom_female_array_int_CPD_95_RWD_ncomp_fixed[,3],

                                 array_int_score_80_male_RWD[,3], array_int_score_95_male_RWD[,3],
                                 geom_male_array_int_CPD_80_RWD[,3], geom_male_array_int_CPD_95_RWD[,3],
                                 array_int_score_80_male_RWD_ncomp_fixed[,3], array_int_score_95_male_RWD_ncomp_fixed[,3],
                                 geom_male_array_int_CPD_80_RWD_ncomp_fixed[,3], geom_male_array_int_CPD_95_RWD_ncomp_fixed[,3])
rownames(geom_RWD_array_int_CPD) = 1:10
colnames(geom_RWD_array_int_CPD) = c("S_F_CPD_80", "S_F_CPD_95",
                                       "G_F_CPD_80", "G_F_CPD_95",
                                       "S_M_CPD_80", "S_M_CPD_95",
                                       "G_M_CPD_80", "G_M_CPD_95",
                                       "S_F_ncomp_fixed_CPD_80", "S_F_ncomp_fixed_CPD_95",
                                       "G_F_ncomp_fixed_CPD_80", "G_F_ncomp_fixed_CPD_95",
                                       "S_M_ncomp_fixed_CPD_80", "S_M_ncomp_fixed_CPD_95",
                                       "G_M_ncomp_fixed_CPD_80", "G_M_ncomp_fixed_CPD_95")
require(xtable)
xtable(geom_RWD_array_int_CPD[,1:8], digits = 3)
round(colMeans(geom_RWD_array_int_CPD[,1:8]), 3) # 0.156 0.050 0.059 0.043 0.187 0.050 0.045 0.023

xtable(geom_RWD_array_int_CPD[,9:16], digits = 3)
round(colMeans(geom_RWD_array_int_CPD[,9:16]), 3) # 0.188 0.050 0.036 0.006 0.194 0.050 0.061 0.012
