##############################################################
# Swedish out-of-sample mortality forecasts from 2021 to 2070
##############################################################

## female

# simple weights

SWE_female_out_CoDa_RWD = R_square_fit(dat = SWE_female_pop, weighting = "simple", fh = 50,
                                       ncomp_selection = "eigen_ratio", forecasting_method = "RWF_drift")$fore_count * 10^5

# geometrically decaying weights

SWE_female_out_CoDa_RWD_geo = R_square_fit(dat = SWE_female_pop, weighting = "geom", geom_weight = 0.07,
                                           fh = 50, ncomp_selection = "eigen_ratio", forecasting_method = "RWF_drift")$fore_count * 10^5

## male

# simple weights

SWE_male_out_CoDa_RWD = R_square_fit(dat = SWE_male_pop, weighting = "simple", fh = 50,
                                       ncomp_selection = "eigen_ratio", forecasting_method = "RWF_drift")$fore_count * 10^5

# geometrically decaying weights

SWE_male_out_CoDa_RWD_geo = R_square_fit(dat = SWE_male_pop, weighting = "geom", geom_weight = 0.1,
                                         fh = 50, ncomp_selection = "eigen_ratio", forecasting_method = "RWF_drift")$fore_count * 10^5

colnames(SWE_female_out_CoDa_RWD) = colnames(SWE_male_out_CoDa_RWD) =
colnames(SWE_female_out_CoDa_RWD_geo) = colnames(SWE_male_out_CoDa_RWD_geo) = 2021:2070
rownames(SWE_female_out_CoDa_RWD) = rownames(SWE_male_out_CoDa_RWD) =
rownames(SWE_female_out_CoDa_RWD_geo) = rownames(SWE_male_out_CoDa_RWD_geo) = 0:110

########################################
# Graphic displays (mortality forecast)
########################################

# female

savepdf("Figure_5a", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(SWE_female_pop)*10^5, xname = "Age", yname = "Life-table death count"),
     main = "Swedish female data", col = gray(0.6), ylim = c(0, 25000))
lines(fts(0:110, SWE_female_out_CoDa_RWD))
legend("topright", c("Observed data (1751-2020)", "Forecasts (2021-2070)"), col = c(gray(0.6), "red"),
       lty = c(1,1), cex = 0.8)
dev.off()

savepdf("Figure_5b", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(SWE_female_pop)*10^5, xname = "Age", yname = "Life-table death count"),
     main = "Swedish female data", col = gray(0.6), ylim = c(0, 25000))
lines(fts(0:110, SWE_female_out_CoDa_RWD_geo))
dev.off()

# male

savepdf("Figure_5c", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(SWE_male_pop)*10^5, xname = "Age", yname = "Life-table death count"),
     main = "Swedish male data", col = gray(0.6), ylim = c(0, 25000))
lines(fts(0:110, SWE_male_out_CoDa_RWD))
dev.off()

savepdf("Figure_5d", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(SWE_male_pop)*10^5, xname = "Age", yname = "Life-table death count"),
     main = "Swedish male data", col = gray(0.6), ylim = c(0, 25000))
lines(fts(0:110, SWE_male_out_CoDa_RWD_geo))
dev.off()

############################################
# Work out age-specific survival count (lx)
############################################

# female & male using RWD method

lx_female_RWD = lx_male_RWD = lx_female_RWD_geo = lx_male_RWD_geo = matrix(NA, 111, 50)
for(ij in 1:50)
{
    for(ik in 1:111)
    {
        # simple

        lx_female_RWD[ik,ij] = 10^5 - sum(SWE_female_out_CoDa_RWD[1:ik, ij])
        lx_male_RWD[ik,ij] = 10^5 - sum(SWE_male_out_CoDa_RWD[1:ik, ij])

        # geo

        lx_female_RWD_geo[ik,ij] = 10^5 - sum(SWE_female_out_CoDa_RWD_geo[1:ik,ij])
        lx_male_RWD_geo[ik,ij] = 10^5 - sum(SWE_male_out_CoDa_RWD_geo[1:ik,ij])
    }
}

#############################################
# Work out survival probability (px) from lx
#############################################

px_female_RWD = px_male_RWD = px_female_RWD_geo = px_male_RWD_geo = matrix(NA, 111, 50)
for(ij in 1:50)
{
    for(ik in 1:111)
    {
        px_female_RWD[,ij] = 1 - round(SWE_female_out_CoDa_RWD[,ij]/c(10^5, lx_female_RWD[1:110,ij]), 4)
        px_male_RWD[,ij] = 1 - round(SWE_male_out_CoDa_RWD[,ij]/c(10^5, lx_male_RWD[1:110,ij]), 4)

        px_female_RWD_geo[,ij] = 1 - round(SWE_female_out_CoDa_RWD_geo[,ij]/c(10^5, lx_female_RWD_geo[1:110,ij]), 4)
        px_male_RWD_geo[,ij] = 1 - round(SWE_male_out_CoDa_RWD_geo[,ij]/c(10^5, lx_male_RWD_geo[1:110,ij]), 4)
    }
}

############################################
# define ages and maturities
# empty cell is when age + maturities > 105
############################################

ages = seq(60, 105, by = 5)
maturities = seq(5, 30, by = 5)

# interest rate = 3% (default setting)

annuities_female_RWD = annuities_male_RWD = annuities_female_RWD_geo = annuities_male_RWD_geo = matrix(NA, length(maturities), length(ages))
for(ij in 1:length(maturities))
{
    for(iw in 1:length(ages))
    {
        annuities_female_RWD[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD[61:110,]), age = ages[iw],
                                                         maturity = maturities[ij], inRate = 0.03), silent = TRUE)

        annuities_male_RWD[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD[61:110,]), age = ages[iw],
                                                       maturity = maturities[ij], inRate = 0.03), silent = TRUE)

        annuities_female_RWD_geo[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD_geo[61:110,]), age = ages[iw],
                                                             maturity = maturities[ij], inRate = 0.03), silent = TRUE)

        annuities_male_RWD_geo[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD_geo[61:110,]), age = ages[iw],
                                                           maturity = maturities[ij], inRate = 0.03), silent = TRUE)
    }
}

# female

annuities_female_mat_RWD = t(annuities_female_RWD)
annuities_female_mat_RWD_geo = t(annuities_female_RWD_geo)

# male

annuities_male_mat_RWD = t(annuities_male_RWD)
annuities_male_mat_RWD_geo = t(annuities_male_RWD_geo)

rownames(annuities_female_mat_RWD) = rownames(annuities_male_mat_RWD) =
rownames(annuities_female_mat_RWD_geo) = rownames(annuities_male_mat_RWD_geo) = ages
colnames(annuities_female_mat_RWD) = colnames(annuities_male_mat_RWD) =
colnames(annuities_female_mat_RWD_geo) = colnames(annuities_male_mat_RWD_geo) = maturities

# female

xtable(annuities_female_mat_RWD, digits = 3)
xtable(annuities_female_mat_RWD_geo, digits = 3)

# male

xtable(annuities_male_mat_RWD, digits = 3)
xtable(annuities_male_mat_RWD_geo, digits = 3)

#####################
# interest rate = 1%
#####################

annuities_female_RWD_0.01 = annuities_male_RWD_0.01 =
annuities_female_RWD_geo_0.01 = annuities_male_RWD_geo_0.01 = matrix(NA, length(maturities), length(ages))
for(ij in 1:length(maturities))
{
    for(iw in 1:length(ages))
    {
        annuities_female_RWD_0.01[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD[61:110,]), age = ages[iw],
                                                                  maturity = maturities[ij], inRate = 0.01), silent = TRUE)

        annuities_male_RWD_0.01[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD[61:110,]), age = ages[iw],
                                                                maturity = maturities[ij], inRate = 0.01), silent = TRUE)

        annuities_female_RWD_geo_0.01[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD_geo[61:110,]), age = ages[iw],
                                                                      maturity = maturities[ij], inRate = 0.01), silent = TRUE)

        annuities_male_RWD_geo_0.01[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD_geo[61:110,]), age = ages[iw],
                                                                    maturity = maturities[ij], inRate = 0.01), silent = TRUE)
    }
}

# female

annuities_female_mat_RWD_0.01 = t(annuities_female_RWD_0.01)
annuities_female_mat_RWD_geo_0.01 = t(annuities_female_RWD_geo_0.01)

# male

annuities_male_mat_RWD_0.01 = t(annuities_male_RWD_0.01)
annuities_male_mat_RWD_geo_0.01 = t(annuities_male_RWD_geo_0.01)

rownames(annuities_female_mat_RWD_0.01) = rownames(annuities_male_mat_RWD_0.01) =
rownames(annuities_female_mat_RWD_geo_0.01) = rownames(annuities_male_mat_RWD_geo_0.01) = ages
colnames(annuities_female_mat_RWD_0.01) = colnames(annuities_male_mat_RWD_0.01) =
colnames(annuities_female_mat_RWD_geo_0.01) = colnames(annuities_male_mat_RWD_geo_0.01) = maturities

# female

xtable(annuities_female_mat_RWD_0.01, digits = 3)
xtable(annuities_female_mat_RWD_geo_0.01, digits = 3)

# male

xtable(annuities_male_mat_RWD_0.01, digits = 3)
xtable(annuities_male_mat_RWD_geo_0.01, digits = 3)

#####################
# interest rate = 5%
#####################

annuities_female_RWD_0.05 = annuities_male_RWD_0.05 =
annuities_female_RWD_geo_0.05 = annuities_male_RWD_geo_0.05 = matrix(NA, length(maturities), length(ages))
for(ij in 1:length(maturities))
{
    for(iw in 1:length(ages))
    {
        annuities_female_RWD_0.05[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD[61:110,]), age = ages[iw],
                                                                  maturity = maturities[ij], inRate = 0.05), silent = TRUE)

        annuities_male_RWD_0.05[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD[61:110,]), age = ages[iw],
                                                                maturity = maturities[ij], inRate = 0.05), silent = TRUE)

        annuities_female_RWD_geo_0.05[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female_RWD_geo[61:110,]), age = ages[iw],
                                                                      maturity = maturities[ij], inRate = 0.05), silent = TRUE)

        annuities_male_RWD_geo_0.05[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male_RWD_geo[61:110,]), age = ages[iw],
                                                                    maturity = maturities[ij], inRate = 0.05), silent = TRUE)
    }
}

# female

annuities_female_mat_RWD_0.05 = t(annuities_female_RWD_0.05)
annuities_female_mat_RWD_geo_0.05 = t(annuities_female_RWD_geo_0.05)

# male

annuities_male_mat_RWD_0.05 = t(annuities_male_RWD_0.05)
annuities_male_mat_RWD_geo_0.05 = t(annuities_male_RWD_geo_0.05)

rownames(annuities_female_mat_RWD_0.05) = rownames(annuities_male_mat_RWD_0.05) =
rownames(annuities_female_mat_RWD_geo_0.05) = rownames(annuities_male_mat_RWD_geo_0.05) = ages
colnames(annuities_female_mat_RWD_0.05) = colnames(annuities_male_mat_RWD_0.05) =
colnames(annuities_female_mat_RWD_geo_0.05) = colnames(annuities_male_mat_RWD_geo_0.05) = maturities

# female

xtable(annuities_female_mat_RWD_0.05, digits = 3)
xtable(annuities_female_mat_RWD_geo_0.05, digits = 3)

# male

xtable(annuities_male_mat_RWD_0.05, digits = 3)
xtable(annuities_male_mat_RWD_geo_0.05, digits = 3)
