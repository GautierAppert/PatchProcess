#! /usr/bin/Rscript

options(guiToolkit="RGtk2")

library(PatchProcess34)

# setwd("../../..")
# PatchProcessWidget()

filename = gfile("Enter filename", type = "open")
print(filename)
