#####################################################
##        ___    ____     __                       ##
##       /   |  / / /__  / /_  __________  _____   ##
##      / /| | / / / _ \/ / / / / ___/ _ \/ ___/   ##
##     / ___ |/ / /  __/ / /_/ (__  )  __/ /       ##
##    /_/  |_/_/_/\___/_/\__, /____/\___/_/        ##
##                      /____/                     ##
#####################################################


# Packages
library(tidyverse)
library(readxl)
library(here)

# Data loading
data <- read_xlsx(here("APCgenotypesAnonym.xlsx"))


# Selecting a SNP
a <- menu(colnames(data[,-1]), title="Please, select the SNP of your interest:") + 1
SNP <- data[, a]
print(paste(colnames(SNP), "selected"))

# Testing for Hardy-Weinberg equilibrium
HWTable <- table(SNP)

suma <- 2 * (HWTable[[1]] + HWTable[[2]] + HWTable[[3]])
p <- (HWTable[[1]] * 2 + HWTable[[2]])/suma
q <- (HWTable[[3]] * 2 + HWTable[[2]])/suma


HWTable <-  rbind(HWTable, c(p^2 , 2*p*q , q^2)* sum(HWTable))
HWTable



