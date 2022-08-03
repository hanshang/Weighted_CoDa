# load R packages

source("load_packages.R")

##########################
# set a working directory
##########################

setwd("~/Dropbox/Todos/cs_geometric/data")

female_qx = t(matrix(read.table("SWE_lt_female_death.txt", header = TRUE)[,4], 111, 270))
male_qx   = t(matrix(read.table("SWE_lt_male_death.txt", header = TRUE)[,4], 111, 270))
female_e0 = matrix(read.table("SWE_lt_female_death.txt", header = TRUE)[,"ex"],111,270)[1,]
male_e0 = matrix(read.table("SWE_lt_male_death.txt", header = TRUE)[,"ex"],111,270)[1,]

n_col = ncol(female_qx)
n_row = nrow(female_qx)
female_pop = male_pop = matrix(NA, n_row, n_col)
for(ij in 1:n_row)
{
    start_pop_female = start_pop_male = 1
    for(ik in 1:n_col)
    {
        female_pop[ij,ik] = female_qx[ij,ik] * start_pop_female
        start_pop_female = start_pop_female - female_pop[ij,ik]

        male_pop[ij,ik] = male_qx[ij,ik] * start_pop_male
        start_pop_male = start_pop_male - male_pop[ij,ik]
        rm(ik)
    }
    print(ij); rm(ij)
}
colnames(female_pop) = colnames(male_pop) = 0:110
rownames(female_pop) = rownames(male_pop) = 1751:2020

# replace zero by small value

SWE_female_pop = replace(female_pop, which(female_pop == 0), 10^-5)
SWE_male_pop   = replace(male_pop,   which(male_pop == 0), 10^-5)

all(round(rowSums(SWE_female_pop), 4) == 1) # TRUE
all(round(rowSums(SWE_male_pop), 4) == 1) # TRUE

# graphical displays

savepdf("Fig_2c", width = 12, height = 10, toplines = 0.5)
plot(1751:2020, female_e0, type = "l", xlab = "Year", ylab = "e(0)")
dev.off()

savepdf("Fig_2d", width = 12, height = 10, toplines = 0.5)
plot(1751:2020, male_e0, type = "l", xlab = "Year", ylab = "e(0)")
dev.off()

savepdf("Fig_2_3_a", width = 12, height = 10, toplines = 0.8)
plot(1751:2020, wei_kappa_0.05, type="l", xlab = "Year", ylab = "Geometrically decaying weights")
lines(1751:2020, wei_kappa_0.95, col = 2, lty = 2)
legend("topleft", c(expression(paste(kappa,"=0.05",sep="")), expression(paste(kappa,"=0.95",sep=""))),
       col = c(1, 2), lty = c(1, 2), cex = 0.8)
dev.off()
