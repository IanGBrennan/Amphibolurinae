library(corrplot)


############################################################################

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

load("Data/Amphibolurinae_Data.RData")

##############################################################################

LSR.reorder <- all.LSR.traits[,c("HeadWidth", "HeadDepth", "PosSkull", "SnoutEye", "EyeDiameter",
                                 "BodyWidth", "Interlimb", "PelvicHeight", "PelvicWidth",
                                 "UpperArm", "LowerArm", "Hand", "UpperLeg", "LowerLeg", "Foot",
                                 "TailWidth", "TailLength",
                                 "Neck", "Size")]


corr.small <- as.matrix(cor(LSR.reorder))
#corrplot::corrplot.mixed(corr.small, lower="ellipse", upper="number")
corrplot(corr.small)
