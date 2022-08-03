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
    return(list(fore_count = d_x_t_star_fore, PI_80 = PI_80, PI_95 = PI_95))
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
    # sample size

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

    # computing interval scores when alpha values = 0.2 and 0.05

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
## fore_method = "ARIMA"
# criterion = "eigenvalue ratio"
#################################

##########################
## criterion = "score_80"
##########################

# female

CoDa_int_geo_weight_ARIMA_female_para_score_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "score_80")
    CoDa_int_geo_weight_ARIMA_female_para_score_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_score_80) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_score_80) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_female_para_score_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "score_80")
    CoDa_int_geo_weight_ARIMA_female_para_score_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_score_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_score_80_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_ARIMA_male_para_score_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom",horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "score_80")
    CoDa_int_geo_weight_ARIMA_male_para_score_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_score_80) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_score_80) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_male_para_score_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom",horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "score_80")
    CoDa_int_geo_weight_ARIMA_male_para_score_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_score_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_score_80_ncomp_fixed) = c("Opt_para", "Opt_score_80")

##########################
## criterion = "score_95"
##########################

# female

CoDa_int_geo_weight_ARIMA_female_para_score_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "score_95")
    CoDa_int_geo_weight_ARIMA_female_para_score_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_score_95) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_score_95) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_female_para_score_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "score_95")
    CoDa_int_geo_weight_ARIMA_female_para_score_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_score_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_score_95_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_ARIMA_male_para_score_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "score_95")
    CoDa_int_geo_weight_ARIMA_male_para_score_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_score_95) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_score_95) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_male_para_score_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "score_95")
    CoDa_int_geo_weight_ARIMA_male_para_score_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_score_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_score_95_ncomp_fixed) = c("opt_para", "opt_obj")

#######################
# criterion = "CPD_80"
#######################

# female

CoDa_int_geo_weight_ARIMA_female_para_CPD_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "CPD_80")
    CoDa_int_geo_weight_ARIMA_female_para_CPD_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_CPD_80) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_CPD_80) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "CPD_80")
    CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_ARIMA_male_para_CPD_80 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "CPD_80")
    CoDa_int_geo_weight_ARIMA_male_para_CPD_80[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_CPD_80) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_CPD_80) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "CPD_80")
    CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed) = c("opt_para", "opt_obj")

#####################################
# Selected optimal tuning parameters
#####################################

# ncomp is determined by eigenvalue ratio

CoDa_int_geo_weight_para = cbind(CoDa_int_geo_weight_ARIMA_female_para_CPD_80[,1],
                                 CoDa_int_geo_weight_ARIMA_female_para_CPD_95[,1],
                                 CoDa_int_geo_weight_RWD_female_para_CPD_80[,1],
                                 CoDa_int_geo_weight_RWD_female_para_CPD_95[,1],
                                 CoDa_int_geo_weight_ARIMA_male_para_CPD_80[,1],
                                 CoDa_int_geo_weight_ARIMA_male_para_CPD_95[,1],
                                 CoDa_int_geo_weight_RWD_male_para_CPD_80[,1],
                                 CoDa_int_geo_weight_RWD_male_para_CPD_95[,1])
colnames(CoDa_int_geo_weight_para) = c("ARIMA_F_CPD_80", "ARIMA_F_CPD_95",
                                       "RWD_F_CPD_80", "RWD_F_CPD_95",
                                       "ARIMA_M_CPD_80", "ARIMA_M_CPD_95",
                                       "RWD_M_CPD_80", "RWD_M_CPD_95")
rownames(CoDa_int_geo_weight_para) = 1:10
xtable(CoDa_int_geo_weight_para, digits = 3)

# ncomp = 6

CoDa_int_geo_weight_para_ncomp_fixed = cbind(CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_RWD_female_para_CPD_80_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_RWD_female_para_CPD_95_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_RWD_male_para_CPD_80_ncomp_fixed[,1],
                                             CoDa_int_geo_weight_RWD_male_para_CPD_95_ncomp_fixed[,1])
colnames(CoDa_int_geo_weight_para_ncomp_fixed) = c("ARIMA_F_CPD_80", "ARIMA_F_CPD_95",
                                       "RWD_F_CPD_80", "RWD_F_CPD_95",
                                       "ARIMA_M_CPD_80", "ARIMA_M_CPD_95",
                                       "RWD_M_CPD_80", "RWD_M_CPD_95")
rownames(CoDa_int_geo_weight_para_ncomp_fixed) = 1:10


#######################
# criterion = "CPD_95"
#######################

# female

CoDa_int_geo_weight_ARIMA_female_para_CPD_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "CPD_95")
    CoDa_int_geo_weight_ARIMA_female_para_CPD_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_CPD_95) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_CPD_95) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_female_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "CPD_95")
    CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed) = c("Opt_para", "Opt_obj")

# male

CoDa_int_geo_weight_ARIMA_male_para_CPD_95 = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "eigen_ratio",
                 fore_method = "ARIMA", criterion = "CPD_95")
    CoDa_int_geo_weight_ARIMA_male_para_CPD_95[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_CPD_95) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_CPD_95) = c("Opt_para", "Opt_obj")

CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed = matrix(NA, 10, 2)
for(iwk in 1:10)
{
    dum <- optim(par = 0.1, fn = CoDa_int_geo_weight, method = "L-BFGS-B",
                 lower = 0, upper = 1, dat = SWE_male_pop,
                 weit_choice = "geom", horizon = iwk, order_selection = "fixed",
                 fore_method = "ARIMA", criterion = "CPD_95")
    CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed[iwk,] = c(dum$par, dum$value)
    rm(dum); print(iwk); rm(iwk)
}
rownames(CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed) = 1:10
colnames(CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed) = c("Opt_para", "Opt_obj")


##############################################
# interval forecasts using the testing sample
##############################################

CoDa_int_test <- function(dat, order_selection, horizon, weit_choice, lambda_para_optim, fore_method)
{
    n = nrow(dat)

    # computing lower and upper bounds

    den_int_80 = den_int_95 = array(NA, dim = c(2, ncol(dat), (11 - horizon)))
    for(iw in 1:(11 - horizon))
    {
        dum <- CoDa_interval_bootstrap(dat = dat[1:(n - 11 + iw),], weighting = weit_choice,
                                       geom_weight = lambda_para_optim,
                                       fh = horizon, ncomp_selection = order_selection,
                                       fmethod = fore_method)
        den_int_80[,,iw] = dum$PI_80 * 10^5
        den_int_95[,,iw] = dum$PI_95 * 10^5
        rm(dum); rm(iw)
    }

    # computing interval scores

    interval_score_80 <- interval_score(holdout = t(dat[(n - 10 + horizon):n,]) * 10^5, lb = den_int_80[1,,], ub = den_int_80[2,,], alpha = 0.2)
    interval_score_95 <- interval_score(holdout = t(dat[(n - 10 + horizon):n,]) * 10^5, lb = den_int_95[1,,], ub = den_int_95[2,,], alpha = 0.05)

    return(list(den_int_80 = den_int_80, den_int_95 = den_int_95,
                interval_score_80 = interval_score_80, interval_score_95 = interval_score_95))
}


###################
# simple weighting
###################

### ncomp_selection = "eigen_ratio"

## fore_method = "ARIMA"

# female

array_int_score_80_ARIMA = array_int_score_95_ARIMA = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "ARIMA")
    array_int_score_80_ARIMA[iwk,] = dum$interval_score_80
    array_int_score_95_ARIMA[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_ARIMA) = rownames(array_int_score_95_ARIMA) = 1:10
colnames(array_int_score_80_ARIMA) = colnames(array_int_score_95_ARIMA) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_ARIMA), 4) # 480.5408 0.9361 0.1361
round(colMeans(array_int_score_95_ARIMA), 4) # 746.4935 0.9996 0.0496

# male

array_int_score_80_ARIMA_male = array_int_score_95_ARIMA_male = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "ARIMA")
    array_int_score_80_ARIMA_male[iwk,] = dum$interval_score_80
    array_int_score_95_ARIMA_male[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_ARIMA_male) = rownames(array_int_score_95_ARIMA_male) = 1:10
colnames(array_int_score_80_ARIMA_male) = colnames(array_int_score_95_ARIMA_male) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_ARIMA_male), 4) # 531.4368 0.9273 0.1273
round(colMeans(array_int_score_95_ARIMA_male), 4) # 799.9029 0.9981 0.0481

### ncomp_selection = "fixed"

## fore_method = "ARIMA"

# female

array_int_score_80_ARIMA_ncomp_fixed = array_int_score_95_ARIMA_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "ARIMA")
    array_int_score_80_ARIMA_ncomp_fixed[iwk,] = dum$interval_score_80
    array_int_score_95_ARIMA_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_ARIMA_ncomp_fixed) = rownames(array_int_score_95_ARIMA_ncomp_fixed) = 1:10
colnames(array_int_score_80_ARIMA_ncomp_fixed) = colnames(array_int_score_95_ARIMA_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_ARIMA_ncomp_fixed), 4) # 405.3441 0.9696 0.1696
round(colMeans(array_int_score_95_ARIMA_ncomp_fixed), 4) # 687.2720 0.9988 0.0488

# male

array_int_score_80_ARIMA_male_ncomp_fixed = array_int_score_95_ARIMA_male_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "simple",
                         fore_method = "ARIMA")
    array_int_score_80_ARIMA_male_ncomp_fixed[iwk,] = dum$interval_score_80
    array_int_score_95_ARIMA_male_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(array_int_score_80_ARIMA_male_ncomp_fixed) = rownames(array_int_score_95_ARIMA_male_ncomp_fixed) = 1:10
colnames(array_int_score_80_ARIMA_male_ncomp_fixed) = colnames(array_int_score_95_ARIMA_male_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(array_int_score_80_ARIMA_male_ncomp_fixed), 4) # 426.2014 0.9801 0.1801
round(colMeans(array_int_score_95_ARIMA_male_ncomp_fixed), 4) # 760.0749 0.9995 0.0495


#################################
# geometrically decaying weights
#################################

### order_selection = "eigen_ratio"

## female

# criterion = "CPD_80"

geom_female_array_int_CPD_80 = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_female_para_CPD_80[iwk,1],
                         fore_method = "ARIMA")
    geom_female_array_int_CPD_80[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_female_array_int_CPD_95 = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_female_para_CPD_95[iwk,1],
                         fore_method = "ARIMA")
    geom_female_array_int_CPD_95[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_female_array_int_CPD_80) = rownames(geom_female_array_int_CPD_95) = 1:10
colnames(geom_female_array_int_CPD_80) = colnames(geom_female_array_int_CPD_95) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_female_array_int_CPD_80), 4) # 856.5859  0.7891 0.0398
round(colMeans(geom_female_array_int_CPD_95), 4) # 1206.2826 0.9548 0.0295

## male

# criterion = "CPD_80"

geom_male_array_int_CPD_80 = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_male_para_CPD_80[iwk,1],
                         fore_method = "ARIMA")
    geom_male_array_int_CPD_80[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_male_array_int_CPD_95 = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "eigen_ratio",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_male_para_CPD_95[iwk,1],
                         fore_method = "ARIMA")
    geom_male_array_int_CPD_95[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_male_array_int_CPD_80) = rownames(geom_male_array_int_CPD_95) = 1:10
colnames(geom_male_array_int_CPD_80) = colnames(geom_male_array_int_CPD_95) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_male_array_int_CPD_80), 4) # 846.9539  0.7864 0.0471
round(colMeans(geom_male_array_int_CPD_95), 4) # 1275.8922 0.9319 0.0269

### order_selection = "fixed"

## female

# criterion = "CPD_80"

geom_female_array_int_CPD_80_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_female_para_CPD_80_ncomp_fixed[iwk,1],
                         fore_method = "ARIMA")
    geom_female_array_int_CPD_80_ncomp_fixed[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_female_array_int_CPD_95_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_female_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_female_para_CPD_95_ncomp_fixed[iwk,1],
                         fore_method = "ARIMA")
    geom_female_array_int_CPD_95_ncomp_fixed[iwk,] = dum$interval_score_95
}
rownames(geom_female_array_int_CPD_80_ncomp_fixed) = rownames(geom_female_array_int_CPD_95_ncomp_fixed) = 1:10
colnames(geom_female_array_int_CPD_80_ncomp_fixed) = colnames(geom_female_array_int_CPD_95_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_female_array_int_CPD_80_ncomp_fixed), 4) #  877.9818 0.7688 0.0825
round(colMeans(geom_female_array_int_CPD_95_ncomp_fixed), 4) # 1220.3792 0.9331 0.0246

## male

# criterion = "CPD_80"

geom_male_array_int_CPD_80_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_male_para_CPD_80_ncomp_fixed[iwk,1],
                         fore_method = "ARIMA")
    geom_male_array_int_CPD_80_ncomp_fixed[iwk,] = dum$interval_score_80
    print(iwk); rm(iwk); rm(dum)
}

# criterion = "CPD_95"

geom_male_array_int_CPD_95_ncomp_fixed = matrix(NA, 10, 3)
for(iwk in 1:10)
{
    dum <- CoDa_int_test(dat = SWE_male_pop, order_selection = "fixed",
                         horizon = iwk, weit_choice = "geom",
                         lambda_para_optim = CoDa_int_geo_weight_ARIMA_male_para_CPD_95_ncomp_fixed[iwk,1],
                         fore_method = "ARIMA")
    geom_male_array_int_CPD_95_ncomp_fixed[iwk,] = dum$interval_score_95
    print(iwk); rm(iwk); rm(dum)
}
rownames(geom_male_array_int_CPD_80_ncomp_fixed) = rownames(geom_male_array_int_CPD_95_ncomp_fixed) = 1:10
colnames(geom_male_array_int_CPD_80_ncomp_fixed) = colnames(geom_male_array_int_CPD_95_ncomp_fixed) = c("Interval score", "ECP", "CPD")

round(colMeans(geom_male_array_int_CPD_80_ncomp_fixed), 4) #  907.9574 0.7318 0.0682
round(colMeans(geom_male_array_int_CPD_95_ncomp_fixed), 4) # 1315.6926 0.9352 0.0177

##########
# summary
##########

geom_ARIMA_array_int_CPD = cbind(array_int_score_80_ARIMA[,3], array_int_score_95_ARIMA[,3],
                                 geom_female_array_int_CPD_80[,3], geom_female_array_int_CPD_95[,3],
                                 array_int_score_80_ARIMA_ncomp_fixed[,3], array_int_score_95_ARIMA_ncomp_fixed[,3],
                                 geom_female_array_int_CPD_80_ncomp_fixed[,3], geom_female_array_int_CPD_95_ncomp_fixed[,3],

                                 array_int_score_80_ARIMA_male[,3], array_int_score_95_ARIMA_male[,3],
                                 geom_male_array_int_CPD_80[,3], geom_male_array_int_CPD_95[,3],
                                 array_int_score_80_ARIMA_male_ncomp_fixed[,3], array_int_score_95_ARIMA_male_ncomp_fixed[,3],
                                 geom_male_array_int_CPD_80_ncomp_fixed[,3], geom_male_array_int_CPD_95_ncomp_fixed[,3])
rownames(geom_ARIMA_array_int_CPD) = 1:10
colnames(geom_ARIMA_array_int_CPD) = c("S_F_CPD_80", "S_F_CPD_95",
                                       "G_F_CPD_80", "G_F_CPD_95",
                                       "S_M_CPD_80", "S_M_CPD_95",
                                       "G_M_CPD_80", "G_M_CPD_95",
                                       "S_F_ncomp_fixed_CPD_80", "S_F_ncomp_fixed_CPD_95",
                                       "G_F_ncomp_fixed_CPD_80", "G_F_ncomp_fixed_CPD_95",
                                       "S_M_ncomp_fixed_CPD_80", "S_M_ncomp_fixed_CPD_95",
                                       "G_M_ncomp_fixed_CPD_80", "G_M_ncomp_fixed_CPD_95")
require(xtable)
xtable(geom_ARIMA_array_int_CPD[,1:8], digits = 3)
round(colMeans(geom_ARIMA_array_int_CPD[,1:8]), 3) # 0.136 0.050 0.040 0.030 0.170 0.049 0.083 0.025

xtable(geom_ARIMA_array_int_CPD[,9:16], digits = 3)
round(colMeans(geom_ARIMA_array_int_CPD[,9:16]), 3) # 0.130 0.048 0.047 0.027 0.180 0.050 0.068 0.018
