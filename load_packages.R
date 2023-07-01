#################
# load R package
#################

# Package names

packages <- c("psych", "ftsa", "meboot", "pracma", "reldist", "flexmix", "demography", "xtable")

# Install packages not yet installed

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading

invisible(lapply(packages, library, character.only = TRUE))
