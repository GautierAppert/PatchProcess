#! /usr/bin/Rscript

#####################################################
#         Olivier Catoni & Gautier Appert           #
#         GNU LESSER GENERAL PUBLIC LICENSE         # 
#####################################################


options(guiToolkit="RGtk2")

library(PatchProcess34)

# setwd("../../..")
# PatchProcessWidget()

filename = gfile("Enter filename", type = "open")
print(filename)
