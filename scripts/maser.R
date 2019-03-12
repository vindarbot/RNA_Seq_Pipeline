library(ggplot2)
library(tidyverse)
library(maser)

dir = setwd("~/Desktop/Data")


path2 <- system.file("extdata", file.path("MATS_output"),
                    package = "maser")


path <- "rMATS_hon4HS/"

hon4 <- maser(path, c("ColHS","hon4_HS"), ftype="JCEC")

head(summary(hon4, type = "SE")[, 1:8])

