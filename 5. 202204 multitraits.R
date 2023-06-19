#### standardization of selected data ####
#### By Jie Hu, 202204

####### put data together for hypervolume analysis ########
## percentage coverage in organic fields: 
## percentage coverage in conventional fields
## Farmers' perception: 
## Weed root mycobiota "all fungi" and each phylum diversity, diversity in FUNGuild data
## Weed ability to transfer root mycobiota, diversity in FUNGuild data
## Growth competition

floristic.crop.s1.o
perception.s2.o<-subset(perception.s2,perception.s2$type!="conventional")
richness.weed
richness.phylum
richness.fun

richness.share.fun
com.ind2

m1<-aggregate(floristic.crop.s1.o,list(floristic.crop.s1.o$ind),FUN = mean)
m2<-aggregate(perception.s2.o,list(perception.s2.o$weed),FUN = mean)
m3<-aggregate(richness.fun,list(richness.fun$weed),FUN = mean)
m4<-aggregate(richness.share.fun,list(richness.share.fun$weed.neighbour),FUN = mean)
m5<-aggregate(com.ind2,list(com.ind2$weed.neighbour),FUN = mean)

write.table(m1,file="m1.txt",sep="\t")
write.table(m2,file="m2.txt",sep="\t")
write.table(m3,file="m3.txt",sep="\t")
write.table(m4,file="m4.txt",sep="\t")
write.table(m5,file="m5.txt",sep="\t")

###### This script will standardize data from 0 to 1 ####
library(devtools)  # for the install_github function
library(rJava)    # if doesn't work reinstall rJava from http://www.java.com/en/download/manual.jsp (Windows Offline (64-bit))
library(reshape)
library(multifunc) # install_github("multifunc", username="jebyrnes")
library(ggplot2)
library(grid)
library(plyr)
library(Hmisc)    # for correlations
library(corrplot) # for correlations
library(nlme)
library(lme4)
library(boot)     # for the inv.logit

##### create a function to standardize each function between 0 and 1 as in Soliveres ####
##### adapted from the getStdAndMeanFunctions ####
standardize01 <- function(afun, minValue=min(afun, na.rm=T), maxValue=max(afun, na.rm=T)){
  (afun - minValue) / (maxValue - minValue)
}

getStdFunction <- function (data, vars, standardizeFunction = standardize01)
{
  ret <- colwise(standardizeFunction)(data[, which(names(data) %in% vars)])
  names(ret) <- paste(names(ret), "_std", sep = "")
  sumFunction <- rowSums(ret, na.rm=T)
  numberFunction <- apply(ret, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]))
  ret$sumFunction <- sumFunction
  ret$numberFunction <- numberFunction
  ret$meanFunction <- ret$sumFunction/ret$numberFunction
  return(ret)
}

# the getFuncsMaxed function does not allow to calculate the number of functions that pass a threshold if there are NAs.
# Adapted from the getFuncsMaxed which contains the getFuncMaxed function to calculate even when there are NAs. The proportion is adjusted as well.
getFuncMaxed_1 <- function (adf, vars = NA, thresh = 0.7, proportion = F, prepend = "Diversity", 
                            maxN = 1) 
{
  if (is.na(vars)[1]) 
    stop("You need to specify some response variable names")
  vars <- whichVars(adf, vars)
  getMaxValue <- function(x) {
    l <- length(x)
    mean(sort(x, na.last = F)[l:(l - maxN + 1)], na.rm = T)
  }  
  funcMaxed <- colwise(function(x) x >= thresh * getMaxValue(x))(adf[, which(names(adf) %in% vars)])
  sumFuncMaxed <- rowSums(colwise(function(x) x >= thresh * getMaxValue(x))(adf[, which(names(adf) %in% vars)]), na.rm = T) # just add , na.rm=T here
  if (proportion) 
    sumFuncMaxed <- sumFuncMaxed/apply(MFdata_std[,which(names(MFdata_std) %in% allVars_std)], MARGIN = 1, FUN = function(x) length(x[!is.na(x)])) # just add the apply function for the denominator to get the exact number of functions that do not have NAs
  ret <- data.frame(cbind(adf[, which(names(adf) %in% prepend)], sumFuncMaxed))
  names(ret) <- c(names(adf)[which(names(adf) %in% prepend)], "sumFuncMaxed")
  ret$nFunc <- apply(adf[, which(names(adf) %in% vars)], MARGIN = 1, FUN = function(x) length(x[!is.na(x)]))
  ret$propFuncMaxed <- ret$sumFuncMaxed / ret$nFunc
  ret$percentFuncMaxed <- 100*ret$propFuncMaxed
  ret <- data.frame(cbind(ret, funcMaxed))
  ret
}

getFuncsMaxed_1 <- function (adf, vars = NA, threshmin = 0.05, threshmax = 0.99, 
                             threshstep = 0.01, proportion = F, prepend = "Diversity", 
                             maxN = 1) 
{
  ret_dummy <- data.frame(thresholds = seq(threshmin, threshmax, 
                                           threshstep))
  ret <- ddply(ret_dummy, .variables = "thresholds", function(x) {
    getFuncMaxed_1(adf, vars = vars, thresh = x[1, 1], proportion = proportion, 
                   maxN = maxN, prepend = c(prepend))  # specify that getFuncMaxed_1 is used instead of getFuncMaxed, otherwise it overwrite with the original
  })
  ret$thresh <- as.numeric(ret$thresh)
  ret
}

#' load the Rfunctions

q <- function(...) {
  sapply(match.call()[-1], deparse)
}

########

#### Which parameters we used for this calculation ####
multi.traits<-read.table("mtotal.txt", sep="\t", header=T, row.names=1) ## read traits
head(multi.traits)

allVars.b<- q(coverage,perception,otu.richness.weed,otu.richness.asco.r,otu.richness.basi.r,otu.richness.chyt.r,otu.richness.glom.r,otu.richness.zygo.r,otu.richness.path,otu.richness.sym,otu.richness.sap,share.num,richness.path,richness.sym,richness.sap,cd)

allVars_std.b<-q(coverage_std,perception_std,otu.richness.weed_std,otu.richness.asco.r_std,otu.richness.basi.r_std,otu.richness.chyt.r_std,otu.richness.glom.r_std,otu.richness.zygo.r_std,otu.richness.path_std,otu.richness.sym_std,otu.richness.sap_std,share.num_std,richness.path_std,richness.sym_std,richness.sap_std,cd_std)

####' calculate the standardized parameters ####
multi_std.b <- cbind(multi.traits, getStdFunction(multi.traits, allVars.b))
head(multi_std.b)

multi_std.b$meanFunction

####' Check the correlation #####
cor_data <- multi_std.b[, which(names(multi_std.b) %in% allVars_std.b)]
res2 <- rcorr(as.matrix(cor_data))
res2
res2$r
res2$P

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

corrplot(res2$r, type="upper",order="hclust")

#make a cluster analysis to chose with functions to group

df_short <- multi_std.b[c("coverage_std","perception_std","otu.richness.weed_std","otu.richness.sym_std",
                          "share.num_std","cd_std")]

short_m<-rcorr(as.matrix(df_short))

corrplot(short_m$r,type="upper",order="hclust",method='number')

d <- dist(t(df_short), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) 

#Reshape for plotting all the functions
head(multi_std.b)
names(multi_std.b)
allVars_std.c<-c("coverage_std","perception_std","cd_std","otu.richness.weed_std","otu.richness.sym_std",
             "share.num_std")

allVars.short<- q(coverage,perception,cd,otu.richness.weed,otu.richness.sym,share.num_std)
allVars_std.short<- q(coverage_std,perception_std,cd_std,otu.richness.weed_std,otu.richness.sym_std,share.num_std)

####' calculate the standardized parameters ####
multi_std.short <- cbind(multi.traits, getStdFunction(multi.traits, allVars.short))
head(multi_std.short)

multi_std.short$meanFunction

df_short <- multi_std.b[c("perception_std","coverage_std","cd_std","otu.richness.weed_std","otu.richness.glom.r_std",
                          "share.num_std")]

#MeanForPlotting.b<-cbind(multi.data,cc)








