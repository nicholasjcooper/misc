# DESCRIPTION: This script conducts a paired t-test analysis for either a case-control or sibling dataset, where there may be
# trios, or quads as well as pairs. It uses 4 alternative approaches for dealing with the uneven cells.
# You should look for converging evidence between the approaches, not just selecting one with the p-values that
# suit the hypothesis. There can be as many columns of paired measurements as you like, they must simply have their
# column headings lists in pairs in the cc.pairs/sb.pairs list(s) below. These were originally TREG/CD4 levels,
# but could be any similar measurement.
# AUTHOR: Nick Cooper
# DATE: February 10, 2015
# FOR: Tony Cutler
library(reader) # use package 'reader'; if this fails, run: install.packages("reader")


########### BEGIN USER MODIFIED INPUT PARAMETERS SECTION ##############

## Sibling dataset location and format
sb.fn <- "/home/ncooper/TonyPairsTriosETC/SibDataTony.csv"
# first column should be 'DAY' with day IDs/numbers; or set value of 'catvar' if the label is different, e.g catvar <- "PAIR"
# subsequent columns should be paired, e.g, unaffected vs affected, for some measurement
# need to list these pairs of column names below as list named 'sb.pairs'
# set the value of 'sb.units.label' below to give the y axis label
# if using siblings, set 'sibling.mode' to TRUE below

## Case control dataset location and format
cc.fn <- "/home/ncooper/TonyPairsTriosETC/CCDataTony.csv"
# first column should be 'PAIR' with pair ids; or set value of 'catvar' if the label is different, e.g catvar <- "DAY"
# subsequent columns should be paired, e.g, case vs control, for some measurement
# need to list these pairs of column names below as list named 'cc.pairs'
# set the value of 'units.label' below to give the y axis label

# location for output plots
plot.dir <- "/home/ncooper/TonyPairsTriosETC/"

sibling.mode <- FALSE # TRUE = analyse sibling dataset, FALSE = analyse Case-Control dataset


## Define column names for measurement pairs to analyse ##

cc.pairs <- list(c("CTRL_TREG","T1D_TREG"),
                 c("CTRL_MEM","T1D_MEM"),
                 c("CTRL_NAIVE","T1D_NAIVE"))

sb.pairs <- list(c("UNAFF_TREG","AFF_TREG"),
                 c("UNAFF_MEM","AFF_MEM"),
                 c("UNAFF_NAIVE","AFF_NAIVE"))

# define units for y-axis labels for each analysis pair respectively
units.label <- c("T_reg","Memory CD4+ T_eff","Naive CD4+ T_eff")

####################### END USER MODIFY SECTION ########################



## NOW WILL APPLY SETTINGS BASED ON WHETHER TO ANALYSE A SIBLING OR CASE-CONTROL DATASET ##
if(sibling.mode) {
  SB <- reader(sb.fn)
  cc.pairs <- sb.pairs; CC <- SB; catvar <- "DAY" 
  units.label <- sb.units.label
} else {
  CC <- reader(cc.fn)
  catvar <- "PAIR"
}

# names for the 4 analysis methods #
analysis.methods <- c("remove last obs.","average matched pairs","duplicate singletons","bootstrap")

## DEFINE MEAN REPLACEMENT FUNCTION ##
mean.rpl <- function(x,def=1) {
  # replace missing values in x with the mean
  if(length(x)>0) { if(length(x)==1 & all(is.na(x))) { x <- def } else { x[is.na(x)] <- mean(x,na.rm=T) }  }
  return(x)
}
######################################



### MAIN PROGRAM ###

RESULTS <- matrix(numeric(),nrow=4,ncol=length(cc.pairs)) # initialise table for summary of p values for each approach
pdf(cat.path(plot.dir,"PairedAnalysisPlots",suf=simple.date(),ext="pdf"),width=12) # initiate pdf for plots to save to

#......................................................................................#
#       ANALYSE WITH 4 ALTERNATIVE APPROACHES FOR DEALING WITH TRIOS/empty-cells       #
#......................................................................................#

#############################################################################
### 1. for each of the trios, remove the 3rd ID, paired t-test 
#############################################################################

CC_1 <- na.omit(CC)

np <- length(cc.pairs)
n.per.row  <- 3  # use 3 plots per row
plot.dims <- c(((np%/%n.per.row)+if(np%%n.per.row!=0) { 1 } else { 0 }),if(np>(n.per.row-1)) { n.per.row } else { np }) 
par(mfrow=plot.dims)

dif <- vector("list",np)
for (i in 1:np) {
  # calculate the difference score for each pair
  dif[[i]] <- CC_1[[ cc.pairs[[i]][2] ]] - CC_1[[ cc.pairs[[i]][1] ]]
}
glob.ylim <- range(sapply(dif,range,na.rm=T),na.rm=T)
pair.labs <- CC_1[[catvar]]; am <- analysis.methods[1]
source("/home/ncooper/TonyPairsTriosETC/tonyPairedPlot.R")
RESULTS[1,] <- result


#############################################################################
### 2. for each of the trios, average the doubles scores, do paired t-test 
#############################################################################

CC_2 <- CC

dif <- vector("list",np)
for (i in 1:np) {
  # calculate the difference score for each pair, averaging when multiples (trios)
  dif[[i]] <- tapply(CC_2[[ cc.pairs[[i]][2] ]],factor(CC_2[[catvar]]),mean,na.rm=T) - tapply(CC_2[[ cc.pairs[[i]][1] ]],factor(CC_2[[catvar]]),mean,na.rm=T)
}
am <- analysis.methods[2]
#glob.ylim <- range(sapply(dif,range,na.rm=T),na.rm=T) # comment out to use the same limits as for part #1
# pair labels as above, as missing from before are now amalgamated into 1 pair group
source("/home/ncooper/TonyPairsTriosETC/tonyPairedPlot.R")
RESULTS[2,] <- result


#############################################################################
### 3. for each of the trios, duplicate the singleton scores, do paired t-test 
#############################################################################

CC_3 <- CC

dif <- vector("list",np)
for (i in 1:np) {
  # calculate the difference score for each pair, averaging when multiples (trios)
  dif[[i]] <- unlist(tapply(CC_3[[ cc.pairs[[i]][2] ]],factor(CC_3[[catvar]]),mean.rpl)) - unlist(tapply(CC_3[[ cc.pairs[[i]][1] ]],factor(CC_3[[catvar]]),mean.rpl))
}
#glob.ylim <- range(sapply(dif,range,na.rm=T),na.rm=T) # comment out to use the same limits as for part #1
pair.labs <- CC_3[[catvar]]; am <- analysis.methods[3]
source("/home/ncooper/TonyPairsTriosETC/tonyPairedPlot.R")
RESULTS[3,] <- result

dev.off() # no plots for bootstrap, so finish up the pdf here


#############################################################################
### 4. bootstapping in 2 strata (pairs, trios), use weight 3/4 for trios
#############################################################################

CC_4 <- CC
# replace missing with neighbours
if(TRUE){
  for (i in 1:np) {
    # averaging when multiples (trios)
    CC_4[,cc.pairs[[i]][1]] <- unlist(tapply(CC_4[[ cc.pairs[[i]][1] ]],factor(CC_4[[catvar]]),mean.rpl)) 
    CC_4[,cc.pairs[[i]][2]] <- unlist(tapply(CC_4[[ cc.pairs[[i]][2] ]],factor(CC_4[[catvar]]),mean.rpl)) 
  }
}
dups <- dup.pairs(CC_4[[catvar]]) # which are duplicated 'PAIR' IDs
ndup <- length(dups) # dup count
nsing <- nrow(CC_4)-ndup  #singleton count
CC_4 <- rbind(CC_4[-dups,],CC_4[dups,]) # organise dataset with singletons first, then multiples (e.g, dups/trios)
w <- c(rep(1, nsing), rep(3/4, ndup)); w <- w/sum(w)

# custom function to calculate the mean difference of each pair of columns (excluding 1st column)
avrg <- function(d, weight=w) {
  n <- (ncol(d)-1)/2; ans <- numeric(n) # assumes 1 header column
  for (i in 1:n) { ans[i] <- sum(apply(d[,((2*i)+0):((2*i)+1)],1,diff)*weight,na.rm=T) }
  return(ans) 
}

# sample 10,000 times using bootstrap
library(boot)
bstp <- boot(CC_4, avrg, R=10000, weights = w, stype = "w", strata = c(rep(1,nsing),rep(2,ndup)))
#source("/home/ncooper/TonyPairsTriosETC/tonyPairedPlot.R")
m = bstp$t0
s = sqrt(diag(var(bstp$t)))
print(bstp)

# calculate the p-values based on bootstrap mean and Std error
for (i in 1:np) {
  t.n = m[i]/s[i]
  pval = 2*pt(abs(t.n), 50, lower.tail=FALSE)
  cat(paste(paste(cc.pairs[[i]],collapse=" vs "),"; p =",round(pval,4)),"\n")
  RESULTS[4,i] <- round(pval,5)
}


#############################################################################
###             FINAL SUMMARY OF P-VALUES ACROSS 4 APPROACHES             ###
#############################################################################

Header("Overall Results")
rownames(RESULTS) <- analysis.methods 
colnames(RESULTS) <- sapply(cc.pairs,function(X) { paste(X,collapse=" vs ") })
print(RESULTS)

#############################################################################
