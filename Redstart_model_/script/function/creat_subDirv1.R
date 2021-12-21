#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
# B.Droz -- February 2017
##########################################

creat.subDir <- function (mainDir,subDir)
{
  if (file.exists(subDir)){
    #setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    #setwd(file.path(mainDir, subDir))
  }
}
################################################