#! /usr/bin/Rscript

#####################################################
#         Olivier Catoni & Gautier Appert           #
#         GNU LESSER GENERAL PUBLIC LICENSE         # 
#####################################################


options(guiToolkit="RGtk2")

library(PatchProcess34)

# don't forget to change the path to your local repository
setwd("/home/gautier/Dropbox/THESE/PatchProcess34/images/")
PatchProcessWidget()

if (FALSE) { 
myList = list(beta = 1.2, gamma = 9.3)
testFunction(myList)
str(myList)
}
