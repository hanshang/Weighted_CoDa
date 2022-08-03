#############################################################
# Age-specific interval mortality forecast from 2021 to 2070
#############################################################

# female

SWE_female_out_CoDa_boot = SWE_female_out_CoDa_boot_geo = array(NA, dim = c(111, 1000, 50))
for(iwk in 1:50)
{
    SWE_female_out_CoDa_boot[,,iwk] = CoDa_interval_bootstrap(dat = SWE_female_pop, weighting = "simple",
                                                fh = iwk, ncomp_selection = "eigen_ratio", fmethod = "ARIMA")$fore_count * 10^5

    SWE_female_out_CoDa_boot_geo[,,iwk] = CoDa_interval_bootstrap(dat = SWE_female_pop, weighting = "geom",
                                                geom_weight = 0.07, fh = iwk,
                                                ncomp_selection = "eigen_ratio", fmethod = "ARIMA")$fore_count * 10^5
    print(iwk); rm(iwk)
}

# male

SWE_male_out_CoDa_boot = SWE_male_out_CoDa_boot_geo = array(NA, dim = c(111, 1000, 50))
for(iwk in 1:50)
{
    SWE_male_out_CoDa_boot[,,iwk] = CoDa_interval_bootstrap(dat = SWE_male_pop, weighting = "simple",
                                                     fh = iwk, ncomp_selection = "eigen_ratio", fmethod = "RWF_drift")$fore_count * 10^5

    SWE_male_out_CoDa_boot_geo[,,iwk] = CoDa_interval_bootstrap(dat = SWE_male_pop, weighting = "geom",
                                                        geom_weight = 0.11, fh = iwk,
                                                       ncomp_selection = "eigen_ratio", fmethod = "RWF_drift")$fore_count * 10^5
    print(iwk); rm(iwk)
}

######################################################
# Work out lx on the basis of forecast death count dx
######################################################

lx_female_boot = lx_female_boot_geo =
lx_male_boot = lx_male_boot_geo = array(NA, dim = c(111, 1000, 50))
for(iw in 1:1000)
{
    for(ij in 1:50)
    {
        for(ik in 1:111)
        {
            lx_female_boot[ik,iw,ij] = 10^5 - sum(SWE_female_out_CoDa_boot[1:ik,iw,ij])
            lx_female_boot_geo[ik,iw,ij] = 10^5 - sum(SWE_female_out_CoDa_boot_geo[1:ik,iw,ij])

            lx_male_boot[ik,iw,ij] = 10^5 - sum(SWE_male_out_CoDa_boot[1:ik,iw,ij])
            lx_male_boot_geo[ik,iw,ij] = 10^5 - sum(SWE_male_out_CoDa_boot_geo[1:ik,iw,ij])
        }
    }
}

######################################################
# Work out survival probability px based on dx and lx
######################################################

px_female_boot = px_female_boot_geo = px_male_boot = px_male_boot_geo = array(NA, dim = c(111, 1000, 50))
for(iw in 1:1000)
{
    for(ij in 1:50)
    {
        for(ik in 1:111)
        {
            px_female_boot[,iw,ij] = 1 - round(SWE_female_out_CoDa_boot[,iw,ij]/c(10^5, lx_female_boot[1:110,iw,ij]),4)
            px_female_boot_geo[,iw,ij] = 1 - round(SWE_female_out_CoDa_boot_geo[,iw,ij]/c(10^5, lx_female_boot_geo[1:110,iw,ij]),4)

            px_male_boot[,iw,ij] = 1 - round(SWE_male_out_CoDa_boot[,iw,ij]/c(10^5, lx_male_boot[1:110,iw,ij]),4)
            px_male_boot_geo[,iw,ij] = 1 - round(SWE_male_out_CoDa_boot_geo[,iw,ij]/c(10^5, lx_male_boot_geo[1:110,iw,ij]),4)
        }
    }
}

###################################
# interest rate eta = 3% (default)
###################################

annuities_female_boot = annuities_female_boot_geo = annuities_male_boot = annuities_male_boot_geo = array(NA, dim = c(length(maturities), 1000, length(ages)))
for(iu in 1:1000)
{
    for(ij in 1:length(maturities))
    {
        for(iw in 1:length(ages))
        {
            annuities_female_boot[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot[61:110,iu,]), age = ages[iw],
                                                                     maturity = maturities[ij], inRate = 0.03), silent = TRUE)
            annuities_female_boot_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot_geo[61:110,iu,]), age = ages[iw],
                                                                     maturity = maturities[ij], inRate = 0.03), silent = TRUE)

            annuities_male_boot[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot[61:110,iu,]), age = ages[iw],
                                                                     maturity = maturities[ij], inRate = 0.03), silent = TRUE)
            annuities_male_boot_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot_geo[61:110,iu,]), age = ages[iw],
                                                                         maturity = maturities[ij], inRate = 0.03), silent = TRUE)
        }
    }
}

female_boot_annuities = female_boot_annuities_geo = male_boot_annuities = male_boot_annuities_geo = array(NA, dim = c(length(ages), 1000, length(maturities)))
for(iw in 1:1000)
{
    female_boot_annuities[,iw,] = t(matrix(as.numeric(annuities_female_boot[,iw,]), nrow = length(maturities)))
    female_boot_annuities_geo[,iw,] = t(matrix(as.numeric(annuities_female_boot_geo[,iw,]), nrow = length(maturities)))

    male_boot_annuities[,iw,] = t(matrix(as.numeric(annuities_male_boot[,iw,]), nrow = length(maturities)))
    male_boot_annuities_geo[,iw,] = t(matrix(as.numeric(annuities_male_boot_geo[,iw,]), nrow = length(maturities)))
}

female_annuities_lb = round(apply(female_boot_annuities, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_ub = round(apply(female_boot_annuities, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

female_annuities_lb_geo = round(apply(female_boot_annuities_geo, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_ub_geo = round(apply(female_boot_annuities_geo, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_no_Inf = replace(male_boot_annuities, which(male_boot_annuities == "-Inf"), NA)
male_annuities_lb = round(apply(male_boot_annuities_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_ub = round(apply(male_boot_annuities_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_geo_no_Inf = replace(male_boot_annuities_geo, which(male_boot_annuities_geo == "-Inf"), NA)
male_annuities_lb_geo = round(apply(male_boot_annuities_geo_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_ub_geo = round(apply(male_boot_annuities_geo_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

rownames(female_annuities_lb) = rownames(female_annuities_ub) =
rownames(female_annuities_lb_geo) = rownames(female_annuities_ub_geo) =
rownames(male_annuities_lb) = rownames(male_annuities_ub) =
rownames(male_annuities_lb_geo) = rownames(male_annuities_ub_geo) = ages

colnames(female_annuities_lb) = colnames(female_annuities_ub) =
colnames(female_annuities_lb_geo) = colnames(female_annuities_ub_geo) =
colnames(male_annuities_lb) = colnames(male_annuities_ub) =
colnames(male_annuities_lb_geo) = colnames(male_annuities_ub_geo) = maturities

# CoDa

xtable(cbind(female_annuities_lb[,1], female_annuities_ub[,1],
             female_annuities_lb[,2], female_annuities_ub[,2],
             female_annuities_lb[,3], female_annuities_ub[,3],
             female_annuities_lb[,4], female_annuities_ub[,4],
             female_annuities_lb[,5], female_annuities_ub[,5],
             female_annuities_lb[,6], female_annuities_ub[,6]), digits = 3)

xtable(cbind(male_annuities_lb[,1], male_annuities_ub[,1],
             male_annuities_lb[,2], male_annuities_ub[,2],
             male_annuities_lb[,3], male_annuities_ub[,3],
             male_annuities_lb[,4], male_annuities_ub[,4],
             male_annuities_lb[,5], male_annuities_ub[,5],
             male_annuities_lb[,6], male_annuities_ub[,6]), digits = 3)

# Weighted CoDa

xtable(cbind(female_annuities_lb_geo[,1], female_annuities_ub_geo[,1],
             female_annuities_lb_geo[,2], female_annuities_ub_geo[,2],
             female_annuities_lb_geo[,3], female_annuities_ub_geo[,3],
             female_annuities_lb_geo[,4], female_annuities_ub_geo[,4],
             female_annuities_lb_geo[,5], female_annuities_ub_geo[,5],
             female_annuities_lb_geo[,6], female_annuities_ub_geo[,6]), digits = 3)

xtable(cbind(male_annuities_lb_geo[,1], male_annuities_ub_geo[,1],
             male_annuities_lb_geo[,2], male_annuities_ub_geo[,2],
             male_annuities_lb_geo[,3], male_annuities_ub_geo[,3],
             male_annuities_lb_geo[,4], male_annuities_ub_geo[,4],
             male_annuities_lb_geo[,5], male_annuities_ub_geo[,5],
             male_annuities_lb_geo[,6], male_annuities_ub_geo[,6]), digits = 3)

#########################
# interest rate eta = 1%
#########################

annuities_female_boot_eta_0.01 = annuities_female_boot_eta_0.01_geo =
annuities_male_boot_eta_0.01 = annuities_male_boot_eta_0.01_geo = array(NA, dim = c(length(maturities), 1000, length(ages)))
for(iu in 1:1000)
{
    for(ij in 1:length(maturities))
    {
        for(iw in 1:length(ages))
        {
            annuities_female_boot_eta_0.01[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot[61:110,iu,]), age = ages[iw],
                                                                     maturity = maturities[ij], inRate = 0.01), silent = TRUE)
            annuities_female_boot_eta_0.01_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot_geo[61:110,iu,]), age = ages[iw],
                                                                         maturity = maturities[ij], inRate = 0.01), silent = TRUE)

            annuities_male_boot_eta_0.01[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot[61:110,iu,]), age = ages[iw],
                                                                   maturity = maturities[ij], inRate = 0.01), silent = TRUE)
            annuities_male_boot_eta_0.01_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot_geo[61:110,iu,]), age = ages[iw],
                                                                       maturity = maturities[ij], inRate = 0.01), silent = TRUE)
        }
    }
}

female_boot_annuities_eta_0.01 = female_boot_annuities_eta_0.01_geo =
male_boot_annuities_eta_0.01 = male_boot_annuities_eta_0.01_geo = array(NA, dim = c(length(ages), 1000, length(maturities)))
for(iw in 1:1000)
{
    female_boot_annuities_eta_0.01[,iw,] = t(matrix(as.numeric(annuities_female_boot_eta_0.01[,iw,]), nrow = length(maturities)))
    female_boot_annuities_eta_0.01_geo[,iw,] = t(matrix(as.numeric(annuities_female_boot_eta_0.01_geo[,iw,]), nrow = length(maturities)))

    male_boot_annuities_eta_0.01[,iw,] = t(matrix(as.numeric(annuities_male_boot_eta_0.01[,iw,]), nrow = length(maturities)))
    male_boot_annuities_eta_0.01_geo[,iw,] = t(matrix(as.numeric(annuities_male_boot_eta_0.01_geo[,iw,]), nrow = length(maturities)))
}

female_annuities_eta_0.01_lb = round(apply(female_boot_annuities_eta_0.01, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_eta_0.01_ub = round(apply(female_boot_annuities_eta_0.01, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

female_annuities_eta_0.01_lb_geo = round(apply(female_boot_annuities_eta_0.01_geo, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_eta_0.01_ub_geo = round(apply(female_boot_annuities_eta_0.01_geo, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_eta_0.01_no_Inf = replace(male_boot_annuities_eta_0.01, which(male_boot_annuities_eta_0.01 == "-Inf"), NA)
male_annuities_eta_0.01_lb = round(apply(male_boot_annuities_eta_0.01_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_eta_0.01_ub = round(apply(male_boot_annuities_eta_0.01_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_eta_0.01_geo_no_Inf = replace(male_boot_annuities_eta_0.01_geo, which(male_boot_annuities_eta_0.01_geo == "-Inf"), NA)
male_annuities_eta_0.01_lb_geo = round(apply(male_boot_annuities_eta_0.01_geo_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_eta_0.01_ub_geo = round(apply(male_boot_annuities_eta_0.01_geo_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

rownames(female_annuities_eta_0.01_lb) = rownames(female_annuities_eta_0.01_ub) =
rownames(female_annuities_eta_0.01_lb_geo) = rownames(female_annuities_eta_0.01_ub_geo) =
rownames(male_annuities_eta_0.01_lb) = rownames(male_annuities_eta_0.01_ub) =
rownames(male_annuities_eta_0.01_lb_geo) = rownames(male_annuities_eta_0.01_ub_geo) = ages

colnames(female_annuities_eta_0.01_lb) = colnames(female_annuities_eta_0.01_ub) =
colnames(female_annuities_eta_0.01_lb_geo) = colnames(female_annuities_eta_0.01_ub_geo) =
colnames(male_annuities_eta_0.01_lb) = colnames(male_annuities_eta_0.01_ub) =
colnames(male_annuities_eta_0.01_lb_geo) = colnames(male_annuities_eta_0.01_ub_geo) = maturities

# CoDa

xtable(cbind(female_annuities_eta_0.01_lb[,1], female_annuities_eta_0.01_ub[,1],
             female_annuities_eta_0.01_lb[,2], female_annuities_eta_0.01_ub[,2],
             female_annuities_eta_0.01_lb[,3], female_annuities_eta_0.01_ub[,3],
             female_annuities_eta_0.01_lb[,4], female_annuities_eta_0.01_ub[,4],
             female_annuities_eta_0.01_lb[,5], female_annuities_eta_0.01_ub[,5],
             female_annuities_eta_0.01_lb[,6], female_annuities_eta_0.01_ub[,6]), digits = 3)

xtable(cbind(male_annuities_eta_0.01_lb[,1], male_annuities_eta_0.01_ub[,1],
             male_annuities_eta_0.01_lb[,2], male_annuities_eta_0.01_ub[,2],
             male_annuities_eta_0.01_lb[,3], male_annuities_eta_0.01_ub[,3],
             male_annuities_eta_0.01_lb[,4], male_annuities_eta_0.01_ub[,4],
             male_annuities_eta_0.01_lb[,5], male_annuities_eta_0.01_ub[,5],
             male_annuities_eta_0.01_lb[,6], male_annuities_eta_0.01_ub[,6]), digits = 3)

# Weighted CoDa

xtable(cbind(female_annuities_eta_0.01_lb_geo[,1], female_annuities_eta_0.01_ub_geo[,1],
             female_annuities_eta_0.01_lb_geo[,2], female_annuities_eta_0.01_ub_geo[,2],
             female_annuities_eta_0.01_lb_geo[,3], female_annuities_eta_0.01_ub_geo[,3],
             female_annuities_eta_0.01_lb_geo[,4], female_annuities_eta_0.01_ub_geo[,4],
             female_annuities_eta_0.01_lb_geo[,5], female_annuities_eta_0.01_ub_geo[,5],
             female_annuities_eta_0.01_lb_geo[,6], female_annuities_eta_0.01_ub_geo[,6]), digits = 3)

xtable(cbind(male_annuities_eta_0.01_lb_geo[,1], male_annuities_eta_0.01_ub_geo[,1],
             male_annuities_eta_0.01_lb_geo[,2], male_annuities_eta_0.01_ub_geo[,2],
             male_annuities_eta_0.01_lb_geo[,3], male_annuities_eta_0.01_ub_geo[,3],
             male_annuities_eta_0.01_lb_geo[,4], male_annuities_eta_0.01_ub_geo[,4],
             male_annuities_eta_0.01_lb_geo[,5], male_annuities_eta_0.01_ub_geo[,5],
             male_annuities_eta_0.01_lb_geo[,6], male_annuities_eta_0.01_ub_geo[,6]), digits = 3)

#########################
# interest rate eta = 5%
#########################

annuities_female_boot_eta_0.05 = annuities_female_boot_eta_0.05_geo =
annuities_male_boot_eta_0.05 = annuities_male_boot_eta_0.05_geo = array(NA, dim = c(length(maturities), 1000, length(ages)))
for(iu in 1:1000)
{
    for(ij in 1:length(maturities))
    {
        for(iw in 1:length(ages))
        {
            annuities_female_boot_eta_0.05[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot[61:110,iu,]), age = ages[iw],
                                                                              maturity = maturities[ij], inRate = 0.05), silent = TRUE)
            annuities_female_boot_eta_0.05_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_boot_geo[61:110,iu,]), age = ages[iw],
                                                                                  maturity = maturities[ij], inRate = 0.05), silent = TRUE)

            annuities_male_boot_eta_0.05[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot[61:110,iu,]), age = ages[iw],
                                                                            maturity = maturities[ij], inRate = 0.05), silent = TRUE)
            annuities_male_boot_eta_0.05_geo[ij,iu,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_boot_geo[61:110,iu,]), age = ages[iw],
                                                                                maturity = maturities[ij], inRate = 0.05), silent = TRUE)
        }
    }
}

female_boot_annuities_eta_0.05 = female_boot_annuities_eta_0.05_geo =
male_boot_annuities_eta_0.05 = male_boot_annuities_eta_0.05_geo = array(NA, dim = c(length(ages), 1000, length(maturities)))
for(iw in 1:1000)
{
    female_boot_annuities_eta_0.05[,iw,] = t(matrix(as.numeric(annuities_female_boot_eta_0.05[,iw,]), nrow = length(maturities)))
    female_boot_annuities_eta_0.05_geo[,iw,] = t(matrix(as.numeric(annuities_female_boot_eta_0.05_geo[,iw,]), nrow = length(maturities)))

    male_boot_annuities_eta_0.05[,iw,] = t(matrix(as.numeric(annuities_male_boot_eta_0.05[,iw,]), nrow = length(maturities)))
    male_boot_annuities_eta_0.05_geo[,iw,] = t(matrix(as.numeric(annuities_male_boot_eta_0.05_geo[,iw,]), nrow = length(maturities)))
}

female_annuities_eta_0.05_lb = round(apply(female_boot_annuities_eta_0.05, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_eta_0.05_ub = round(apply(female_boot_annuities_eta_0.05, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

female_annuities_eta_0.05_lb_geo = round(apply(female_boot_annuities_eta_0.05_geo, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_eta_0.05_ub_geo = round(apply(female_boot_annuities_eta_0.05_geo, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_eta_0.05_no_Inf = replace(male_boot_annuities_eta_0.05, which(male_boot_annuities_eta_0.05 == "-Inf"), NA)
male_annuities_eta_0.05_lb = round(apply(male_boot_annuities_eta_0.05_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_eta_0.05_ub = round(apply(male_boot_annuities_eta_0.05_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

male_boot_annuities_eta_0.05_geo_no_Inf = replace(male_boot_annuities_eta_0.05_geo, which(male_boot_annuities_eta_0.05_geo == "-Inf"), NA)
male_annuities_eta_0.05_lb_geo = round(apply(male_boot_annuities_eta_0.05_geo_no_Inf, c(1,3), quantile, 0.025, na.rm = TRUE), 3)
male_annuities_eta_0.05_ub_geo = round(apply(male_boot_annuities_eta_0.05_geo_no_Inf, c(1,3), quantile, 0.975, na.rm = TRUE), 3)

rownames(female_annuities_eta_0.05_lb) = rownames(female_annuities_eta_0.05_ub) =
rownames(female_annuities_eta_0.05_lb_geo) = rownames(female_annuities_eta_0.05_ub_geo) =
rownames(male_annuities_eta_0.05_lb) = rownames(male_annuities_eta_0.05_ub) =
rownames(male_annuities_eta_0.05_lb_geo) = rownames(male_annuities_eta_0.05_ub_geo) = ages

colnames(female_annuities_eta_0.05_lb) = colnames(female_annuities_eta_0.05_ub) =
colnames(female_annuities_eta_0.05_lb_geo) = colnames(female_annuities_eta_0.05_ub_geo) =
colnames(male_annuities_eta_0.05_lb) = colnames(male_annuities_eta_0.05_ub) =
colnames(male_annuities_eta_0.05_lb_geo) = colnames(male_annuities_eta_0.05_ub_geo) = maturities

# CoDa

xtable(cbind(female_annuities_eta_0.05_lb[,1], female_annuities_eta_0.05_ub[,1],
             female_annuities_eta_0.05_lb[,2], female_annuities_eta_0.05_ub[,2],
             female_annuities_eta_0.05_lb[,3], female_annuities_eta_0.05_ub[,3],
             female_annuities_eta_0.05_lb[,4], female_annuities_eta_0.05_ub[,4],
             female_annuities_eta_0.05_lb[,5], female_annuities_eta_0.05_ub[,5],
             female_annuities_eta_0.05_lb[,6], female_annuities_eta_0.05_ub[,6]), digits = 3)

xtable(cbind(male_annuities_eta_0.05_lb[,1], male_annuities_eta_0.05_ub[,1],
             male_annuities_eta_0.05_lb[,2], male_annuities_eta_0.05_ub[,2],
             male_annuities_eta_0.05_lb[,3], male_annuities_eta_0.05_ub[,3],
             male_annuities_eta_0.05_lb[,4], male_annuities_eta_0.05_ub[,4],
             male_annuities_eta_0.05_lb[,5], male_annuities_eta_0.05_ub[,5],
             male_annuities_eta_0.05_lb[,6], male_annuities_eta_0.05_ub[,6]), digits = 3)

# Weighted CoDa

xtable(cbind(female_annuities_eta_0.05_lb_geo[,1], female_annuities_eta_0.05_ub_geo[,1],
             female_annuities_eta_0.05_lb_geo[,2], female_annuities_eta_0.05_ub_geo[,2],
             female_annuities_eta_0.05_lb_geo[,3], female_annuities_eta_0.05_ub_geo[,3],
             female_annuities_eta_0.05_lb_geo[,4], female_annuities_eta_0.05_ub_geo[,4],
             female_annuities_eta_0.05_lb_geo[,5], female_annuities_eta_0.05_ub_geo[,5],
             female_annuities_eta_0.05_lb_geo[,6], female_annuities_eta_0.05_ub_geo[,6]), digits = 3)

xtable(cbind(male_annuities_eta_0.05_lb_geo[,1], male_annuities_eta_0.05_ub_geo[,1],
             male_annuities_eta_0.05_lb_geo[,2], male_annuities_eta_0.05_ub_geo[,2],
             male_annuities_eta_0.05_lb_geo[,3], male_annuities_eta_0.05_ub_geo[,3],
             male_annuities_eta_0.05_lb_geo[,4], male_annuities_eta_0.05_ub_geo[,4],
             male_annuities_eta_0.05_lb_geo[,5], male_annuities_eta_0.05_ub_geo[,5],
             male_annuities_eta_0.05_lb_geo[,6], male_annuities_eta_0.05_ub_geo[,6]), digits = 3)
