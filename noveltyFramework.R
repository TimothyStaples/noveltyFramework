# ========================================================================= ####
# Title:        A conceptual framework for measuring ecological novelty     ####    
# Description:  This script performs analyses for the quantitative case study
# Author:       Timothy Staples
# Date Created: 2022-01-01
# Last Modified: 2024-06-11
# Version:      1.0
# License:      MIT License
# ========================================================================= #### 
# SET DIRECTORY ####
setwd("/Users/uqtstapl/Library/CloudStorage/Dropbox/Tim/Post-doc/Research projects/PaleoNovelty/prodCode")

# READ IN SUPPLEMENTARY FUNCTIONS ####
sapply(list.files("./functions", full.names = TRUE), source)

# DEPENDENCIES ####
package.loader(c("rworldmap", "vegan", "sf", "shape", "lme4"))

# Dependencies:
# - rworldmap: world map shape files
# - vegan: ecological analysis and dissimilarity indices
# - sf: spatial data capability
# - shape: supplementary plotting tools
# - lme4: mixed-effects modelling

# Notes:
# Neotoma download and taxonomic synonymization are conducted in a separate script
# (neotomaProcessing.R). The output of the processing script has been pre-loaded
# into the "rawdata" subfolder of this repository for ease of reproduction.

# ========================================================================= ####

# 1. NEOTOMA DATA IMPORT ####

# read in processed Neotoma data
neo <- readRDS("./rawdata/processedRecords.rds")

# only pollen (exclude spores etc)
neo = neo[neo$elementtype == "pollen",]

# Collate coordinates and geographical data using world map and lat/long coords
# restrict data to Nth America and pollen samples within -50-5000 ybp
neo <- neo[complete.cases(neo[,c("long","lat")]),]
world <- getMap("high")
site.df <- neo[!duplicated(neo$siteid), c("age", "siteid", "sitename", "long", "lat", "elev")]
site.coords.raw <- site.df[,c("long", "lat")]
coordinates(site.coords.raw) <- c("long", "lat")
proj4string(site.coords.raw) <- proj4string(world)
site.df.save <- site.df
site.df <- cbind(site.df, site.coords.raw %over% world)
nthAmSite <- site.df[site.df$REGION == "North America",]
neo <- droplevels(neo[neo$siteid %in% nthAmSite$siteid,])
neo <- neo[neo$age <= 5000 & neo$age > -50 & !is.na(neo$age),]

#               Identify high density sample region ####

# North american plot
plot(nthAmSite$lat ~ nthAmSite$long)
# target region
rect(xleft=-97.5, xright=-90, ybottom=41, ytop=49, border="red")

# subset site data to target region
subSite <- droplevels(nthAmSite[nthAmSite$long <= -90 & nthAmSite$long >= -97.5 &
                                 nthAmSite$lat <= 49 & nthAmSite$lat >= 41,])
plot(subSite$lat ~ subSite$long)

#               Sub-sample data ####

# subset sample data to just sites in target region
neoSub <- droplevels(neo[neo$siteid %in% subSite$siteid,])

# add in geographical metadata to sample dataframe
neoSub <- merge(neoSub, subSite[,-c(1,3:6)],
             by.x="siteid", by.y="siteid", all.x=TRUE, all.y=FALSE, sort=FALSE)

# cut pollen ages into 500 year bins (beginning at -50 ybp which = 2000 AD) 
neoSub$timeBin <- cutCenter(neoSub$age, 
                            breaks=seq(-50,5500,500),
                            include.lowest = TRUE)

# create site x taxa x timebin array (3D array)
neoArr <- tapply(neoSub$value,
                 list(neoSub$siteid, neoSub$family, neoSub$timeBin), sum, na.rm=TRUE)
neoArr[is.na(neoArr)] = 0
dim(neoArr) # dimensions are ordered as site, taxa, time

# calculate the number of samples within each bin
neoSampleSize = with(neoSub[!duplicated(neoSub[,c("siteid", "sampleid")]),],
                                        table(siteid, timeBin))

# remove aquatic families
neoArr <- neoArr[,!dimnames(neoArr)[[2]] %in% c("Potamogetonaceae", "Nymphaeaceae", 
                                                "Typhaceae", "Cabombaceae", 
                                                "Alismataceae", "Haloragaceae"),]

# Use only the 20 most abundant families. First create an array where each
# site x time bin combination is relativized
neoArrProp = array(apply(neoArr, 3, function(x){prop.table(x, margin=1)}),
                   dim=dim(neoArr), dimnames=dimnames(neoArr))

# remove taxa with abundance < the 20th highest mean
taxaMean = apply(neoArrProp[,,1:2], 2, mean, na.rm=TRUE)
neoArr = neoArr[,taxaMean >= sort(taxaMean, decreasing=TRUE)[20],]

#               Calculate novelty ####

# now we can make novelty assessments across modern vs pre-modern time slices

# subset array to time series with slices in the modern (-50-450 ybp) and pre-modern
# (450-950 ybp). Look for sites with data in both time bins

neoArrSub = neoArr[rowSums(neoArr[,,1] > 0) & rowSums(neoArr[,,2] >0),,1:2]
neoArrRoot = sqrt(neoArrSub)

novGrid <- data.frame(siteid = rownames(neoArrSub))
novGrid <- merge(novGrid, subSite[,c("siteid", "long", "lat")],
                 by.x="siteid", by.y="siteid", all.x=TRUE, all.y=FALSE, sort=FALSE)

# novelty calculated four ways:

# how the past differs from the present (no-analog)
novGrid = cbind(novGrid, 
                do.call('rbind', lapply(1:dim(neoArrRoot)[1], function(n){
                  x = novelty(neoArrRoot[n,,2] / sum(neoArrRoot[n,,2]), 
                              prop.table(neoArrRoot[,,1], 1), method="bray")
                  
                  return(data.frame(noAnalog = x,
                                    noAnalogWhich = names(x)))
                })))
                
# how the present differs from the past (time's arrow novelty)
novGrid = cbind(novGrid, 
                do.call("rbind", lapply(1:dim(neoArrRoot)[1], function(n){
                  x = novelty(neoArrRoot[n,,1] / sum(neoArrRoot[n,,1]), 
                              prop.table(neoArrRoot[,,2], 1), method="bray")
                  return(data.frame(timeArrow = x,
                                    timeArrowWhich = names(x)))
                })))

# how the past differs from the past (past comparison)
novGrid = cbind(novGrid,
                do.call("rbind", lapply(1:dim(neoArrRoot)[1], function(n){
 x = novelty(neoArrRoot[n,,2] / sum(neoArrRoot[n,,2]), prop.table(neoArrRoot[,,2], 1), method="bray", nSize=2)[2]
 return(data.frame(pastComp = x,
                   pastCompWhich = names(x)))
 })))

# how the present differs from the present (present comparison)
novGrid = cbind(novGrid,
                do.call("rbind", lapply(1:dim(neoArrRoot)[1], function(n){
  x = novelty(neoArrRoot[n,,1] / sum(neoArrRoot[n,,1]), prop.table(neoArrRoot[,,1], 1), method="bray", nSize=2)[2]
  return(data.frame(presentComp = x,
                    presentCompWhich = names(x)))
  })))

# add in sampling size
novGrid$pastSamples = neoSampleSize[match(novGrid$siteid, rownames(neoSampleSize)), 2]
novGrid$presentSamples = neoSampleSize[match(novGrid$siteid, rownames(neoSampleSize)), 1]

# re-order columns so novelty columns are all together
novGrid = novGrid[,c("siteid", "long", "lat", 
                     "noAnalog", "timeArrow", "pastComp", "presentComp",
                     "noAnalogWhich", "timeArrowWhich", "pastCompWhich", "presentCompWhich",
                     "pastSamples", "presentSamples")]

# does sample size correlate with novelty measures?
cor(cbind(sapply(novGrid[,4:7], logit),
          sapply(novGrid[,12:13], log)), 
          use="complete.obs")

# FIGURES ####
#         novelty PCA ####

# novel correlations (logit transformed)
cor(sapply(novGrid[,4:7], logit), use="complete.obs")

# re-order data based on overall novelty sum (lowest to highest)
novGrid = novGrid[order(rowSums(novGrid[,4:7])),]

# principal component analysis
a = princomp(sapply(novGrid[,4:7], logit))
aVar = a$sdev / sum(a$sdev)
#a = princomp(novGrid[,4:7]) # does not logit-transform dissimilarities

# highlight points with at least one novelty value over threshold
novelCat = rowSums(novGrid[,4:7] > 0.15) > 0

# save novelty data
novGrid$fig3Number = ""
novGrid$fig3Number[novelCat] = sum(novelCat):1
write.csv(novGrid, "./outputs/novelResults.csv")

pdf("./plots/noveltyPCA.pdf", height=3.5, width=7, useDingbats = FALSE)
split.screen(rbind(c(0.1,0.49,0.15,0.99),
                   c(0.6,0.99,0.15,0.99)))
vectorScale = 1

col1 = rgb(colorRamp(c("white", "grey85", "black"), bias=2)(unitScale(a$scores[,1]))/255)
col2 = rgb(colorRamp(c("purple1", "orange"), bias=2)(unitScale(a$scores[,2]))/255)
col3 = rgb(colorRamp(c("red", "white", "blue"))(unitScale(a$scores[,3], custMin= -max(abs(a$scores[,3])), custMax=max(abs(a$scores[,3]))))/255)
col4 = rgb(colorRamp(c("yellow", "white", "darkgreen"))(unitScale(a$scores[,4], custMin= -max(abs(a$scores[,4])), custMax=max(abs(a$scores[,4]))))/255)

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(a$scores[,1:2], asp=1, type="n", axes=FALSE, xlab="", ylab="")
axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1.25, text=paste0("PC1 (", sprintf("%.2f", aVar[1]*100), "%)"), las=0)
axis(side=2); mtext(side=2, line=2, text=paste0("PC2 (", sprintf("%.2f", aVar[2]*100), "%)"), las=0)
Arrows(x0=0, y0=0, col="grey50",
       x1=a$loadings[,1]*vectorScale, 
       y1=a$loadings[,2]*vectorScale, 
       arr.type="triangle", arr.length=0.15, arr.width=0.15, lwd=2)
text(x=a$loadings[,1]*vectorScale, y=a$loadings[,2]*vectorScale, 
     labels=c("No-analog", "Time's arrow", "Past comp", "Present comp"),
     adj=0, col="grey50", pos=4)

points(a$scores[,1:2], pch=c(21,22)[as.factor(novelCat)], lwd=0.5,
       bg=mapply(c1=col1, c2=col2, colorMix), cex=ifelse(novelCat, 1.5, 1))
text(a$scores[novelCat,1:2], labels=sum(novelCat==TRUE):1, cex=0.575, col="white")

text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
     labels="(A)", font=1, adj=0)
box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(a$scores[,3:4], asp=1, type="n", axes=FALSE, xlab="", ylab="", ylim=c(-0.8,0.8), xlim=c(-1,1))
axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1.25, text=paste0("PC3 (", sprintf("%.2f", aVar[3]*100), "%)"), las=0)
axis(side=2); mtext(side=2, line=2, text=paste0("PC4 (", sprintf("%.2f", aVar[4]*100), "%)"), las=0)
Arrows(x0=0, y0=0, col="grey50",
       x1=a$loadings[,3]*vectorScale, 
       y1=a$loadings[,4]*vectorScale, 
       arr.type="triangle", arr.length=0.15, arr.width=0.15, lwd=2)
text(x=a$loadings[,3]*vectorScale, y=a$loadings[,4]*vectorScale, 
     labels=c("No-analog", "Time's arrow", "Past comp", "Present comp"),
     adj=0, col="grey50", pos=c(3,1,2,3))
points(a$scores[,3:4], pch=c(21,22)[as.factor(novelCat)],  lwd=0.5,
       bg=mapply(c1=col3, c2=col4, colorMix), cex=ifelse(novelCat, 1.5, 1))
text(a$scores[novelCat,3:4], labels=sum(novelCat==TRUE):1, cex=0.575, 
     col="black")
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
     labels="(B)", font=1, adj=0)
box()
close.screen(2)
dev.off()

#         study map ####

usOutline = read_sf("./rawdata/shapeFiles/cb_2018_us_nation_20m/cb_2018_us_nation_20m.shp")
usStates = read_sf("./rawdata/shapeFiles/cb_2018_us_state_5m/cb_2018_us_state_5m.shp")
worldmap = countriesLow

pdf("./plots/studyMap.pdf", height=4.5, width=3.75, useDingbats = FALSE)
split.screen(rbind(c(0.15,0.95,0.05,0.99),
                   c(0.75,0.95,0.05,0.2)))
vectorScale = 1

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(NULL, xlim=c(-97.5,-90), ylim=c(41,49), xlab="", ylab="", axes=FALSE)
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     col="grey90")
axis(side=1, at=seq(-100,-90,2), labels=parse(text=paste0(abs(rev(seq(90,100,2))), "*degree~W")), mgp=c(3,0.2,0))
axis(side=2, mgp=c(3,0.5,0),
     at=seq(40,50,2), labels=parse(text=paste0(seq(40,50,2), "*degree~N")))

subStates = st_crop(usStates, xmin=-102, xmax=-85, ymin=35, ymax=55)
plot(subStates, max.plot=1, col="white", bg="grey80", border="black", add=TRUE)

points(novGrid$lat ~ novGrid$long, pch=21, lwd=0.5,
       bg="grey50")
box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), mgp=c(3,0,0), tcl=-0.25)
plot(worldmap, max.plot=1, xlim=c(-125,-65), ylim=c(35,45), col="grey80", border="grey80", bg="white")
plot(usOutline, max.plot=1, add=TRUE, col="grey50")
rect(xleft=-97.5, xright=-90, ybottom=41, ytop=49, border="red", lwd=2)
box()
close.screen(2)
close.screen(all.screens=TRUE)

dev.off()

#         study map with novelty overlaid ####

pdf("./plots/studyMapNovelty.pdf", height=8.35, width=6, useDingbats = FALSE)
split.screen(rbind(c(0.1,0.545,0.52,0.99),
                   c(0.545,0.99,0.52,0.99),
                   c(0.1,0.545,0.05,0.52),
                   c(0.545,0.99,0.05,0.52),
                   c(0.375,0.545,0.47,0.55),
                   c(0.12,0.32,0.55,0.58),
                   c(0.565,0.765,0.55,0.58),
                   c(0.12,0.32,0.08,0.11),
                   c(0.565,0.765,0.08,0.11)))

novMat = novGrid[,c("noAnalog", "timeArrow", "pastComp", "presentComp")]

sapply(1:4, function(n){
  
screen(n)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(NULL, xlim=c(-97.5,-90), ylim=c(41,49.75), xlab="", ylab="", axes=FALSE)
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     col="grey90")

subStates = st_crop(usStates, xmin=-102, xmax=-85, ymin=35, ymax=55)
plot(subStates, max.plot=1, col="white", bg="grey80", border="black", add=TRUE)

novCol = c("darkblue", "red", "darkgreen", "purple")[n]
novRamp = colorRampPalette(c("white", novCol))
points(novGrid$lat ~ novGrid$long, cex=1.25,
       pch=c(21,22)[as.factor(novMat[,n] > 0.15)], 
       lwd=0.5,
       bg=novRamp(4)[cut(ifelse(novMat[,n] > 0.2, 0.2, novMat[,n]), breaks=seq(0,0.20,0.05))])
box()
text(novGrid$lat ~ novGrid$long,
     labels=novGrid$fig3Number, col="white", cex=0.5)


if(n %in% c(3:4)){
  axis(side=1, at=seq(-100,-90,2), labels=parse(text=paste0(abs(rev(seq(90,100,2))), "*degree~W")), mgp=c(3,0.2,0))
  axis(side=3, at=seq(-100,-90,2), labels=NA, tcl=0.25)
}
if(n %in% c(1,3)){
  axis(side=2, mgp=c(3,0.5,0), at=seq(40,50,2), labels=parse(text=paste0(seq(40,50,2), "*degree~N")))
  axis(side=4, mgp=c(3,0.5,0), at=seq(40,50,2), labels=NA, tcl=0.25)
}

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj = 0, labels=paste0("(", LETTERS[n], ")"), font=2)
text(x=relative.axis.point(0.11, "x"),
     y=relative.axis.point(0.95, "y"),
     adj = 0, labels=c("No-analog", "Time's arrow", "Past comparison", "Present comparison")[n])
rect(xleft=par("usr")[1], xright=relative.axis.point(0.5,"x"), ybottom=par("usr")[3], ytop=relative.axis.point(0.1,"y"),
     col="white", border=NA)

close.screen(n)


screen(5+n)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.1,0), las=1)
plot.new()
gridX = seq(par('usr')[1], par("usr")[2], len=5)
rect(xleft=gridX[-length(gridX)], xright=gridX[-1],
     ybottom=par("usr")[3], ytop=par("usr")[4],
     col=novRamp(4))
axis(side=1, at=gridX[c(1,3,5)], labels=seq(0,0.2,0.1))
box()
close.screen(5+n)

})

screen(5)
par(mar=c(0,0,0,0), mgp=c(3,0,0), tcl=-0.25)
plot(worldmap, max.plot=1, xlim=c(-125,-65), ylim=c(35,45), col="grey80", border="grey80", bg="white")
plot(usOutline, max.plot=1, add=TRUE, col="grey50")
rect(xleft=-97.5, xright=-90, ybottom=41, ytop=49, border="red", lwd=2)
box()
close.screen(5)
close.screen(all.screens=TRUE)
dev.off()

#         nMDS plot ####

# run ordination on relativized counts
neoArrSubProp = array(apply(neoArrSub, 3, function(x){prop.table(x, margin=1)}),
                   dim=dim(neoArrSub), dimnames=dimnames(neoArrSub))

# collapse 3D array to site-taxa matrix
neoNMDSArr = rbind(neoArrSubProp[,,1], neoArrSubProp[,,2])

set.seed(001155)
neoOrd <- metaMDS(neoNMDSArr)

pdf("./plots/newOrd.pdf", height=8, width=8, useDingbats = FALSE)
par(mfrow=c(2,2), oma=c(3,3,0.5,0.5), mar=c(0,0,0,0), mgp=c(3,0.5,0), ps=10, tcl=-0.25, las=1)

sapply(1:4, function(n){

plot(neoOrd$points, type="n", asp=1, xlab="", ylab="", xaxt="n", yaxt="n")
  
if(n %in% c(3:4)){
axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1, text="nMDS1")
} else {axis(side=1, labels=NA)}
  
if(n %in% c(1,3)){
axis(side=2); mtext(side=2, line=2, las=0, text="nMDS2")
} else {axis(side=2, labels=NA)}

novCol = c("darkblue", "red", "darkgreen", "purple")[n]
novRamp = colorRampPalette(c("white", novCol))(5)

points(neoOrd$species, pch=16, lwd=2, col="grey")
text(neoOrd$species, labels=gsub("aceae|inosae|ositae", "", rownames(neoOrd$species)), pos=2, col="grey")

presentPoints = neoOrd$points[1:dim(neoArrSubProp)[1],]
pastPoints = neoOrd$points[(dim(neoArrSubProp)[1]+1):nrow(neoOrd$points),]
novGridOrder = novGrid[match(rownames(presentPoints), novGrid$siteid),]

novScaling = function(x){c(0.5,0.8,1,1.5, 1.5)[cut(x, breaks=c(0,0.05,0.1,0.15,0.2,1))]}

# no-analog points
if(n==1){
  
  segments(x0 = pastPoints[,1],
           y0 = pastPoints[,2],
           x1 = presentPoints[match(novGridOrder$noAnalogWhich, rownames(presentPoints)), 1],
           y1 = presentPoints[match(novGridOrder$noAnalogWhich, rownames(presentPoints)), 2],
           lty="31")
  
  points(presentPoints, pch=4, col="black", cex=0.5)

  points(pastPoints, 
         pch=21, cex=novScaling(novGridOrder$noAnalog), 
         bg=novRamp[cut(novGridOrder$noAnalog, breaks=seq(0,0.25,0.05))])
  
}

# time's arrow points
if(n==2){
  
  segments(x0= presentPoints[,1],
           y0= presentPoints[,2],
           x1 = pastPoints[match(novGridOrder$timeArrowWhich, rownames(pastPoints)), 1],
           y1 = pastPoints[match(novGridOrder$timeArrowWhich, rownames(pastPoints)), 2],
           lty="31")
  
  points(pastPoints, pch=4, col="black", cex=0.5)
  
  
  points(presentPoints, 
         pch=21, cex=novScaling(novGridOrder$timeArrow), 
         bg=novRamp[cut(novGridOrder$timeArrow, breaks=seq(0,0.25,0.05))])
  
}

# past comparison points
if(n==3){
  
  segments(x0 = pastPoints[,1],
           y0 = pastPoints[,2],
           x1 = pastPoints[match(novGridOrder$pastCompWhich, rownames(pastPoints)), 1],
           y1 = pastPoints[match(novGridOrder$pastCompWhich, rownames(pastPoints)), 2],
           lty="31")
  
  points(pastPoints, 
         pch=21, cex=novScaling(novGridOrder$pastComp), 
         bg=novRamp[cut(novGridOrder$pastComp, breaks=seq(0,0.25,0.05))])
  
  
  legend(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(0.9, "y"),
         pch=c(rep(21,4),4), pt.cex=novScaling(seq(0.025, 0.175, 0.05)),
         pt.bg=c(colorRampPalette(c("white", "grey20"))(4), NA),
         legend = c("< 0.05", "0.05 - 0.10","0.10 - 0.15","> 0.15", "Ref set"),
         title = "Novelty")
  
}

if(n==4){
  
  segments(x0 = presentPoints[,1],
           y0 = presentPoints[,2],
           x1 = presentPoints[match(novGridOrder$presentCompWhich, rownames(presentPoints)), 1],
           y1 = presentPoints[match(novGridOrder$presentCompWhich, rownames(presentPoints)), 2],
           lty="31")
  
  points(presentPoints, 
         pch=21, cex=novScaling(novGridOrder$presentComp), 
         bg=novRamp[cut(novGridOrder$presentComp, breaks=seq(0,0.25,0.05))])
  
}

text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
     labels=paste0("(",LETTERS[n],")"), font=2, adj=0)
text(x=relative.axis.point(0.07, "x"), y=relative.axis.point(0.95, "y"),
     labels=c("No-analog", "Time's arrow", "Past comparison", "Present comparison")[n], 
     adj=0)

text(x=relative.axis.point(0.035, "x"), y=relative.axis.point(0.035, "y"),
     labels=paste0("Stress = ", sprintf("%.3f", neoOrd$stress)), adj=0)

})
dev.off()

# DIMENSIONALITY TEST ####

# randomly order taxa so we can draw them in sequence and calculate how novelty
# changes with increased dimensionality.
taxaSamp = t(replicate(999, sample(1:20, 20)))

# for each random order, calculate novelty successively across increasing
# numbers of taxa.
taxaDimNov = lapply(1:nrow(taxaSamp), function(topN){

  x = taxaSamp[topN,]
  
  # reorder taxa matrices
  dimArr = neoArrSub[,x,]
  
  # run no-analog test for successive taxa, return novelty matrix
  dimNovMat = sapply(2:dim(dimArr)[2], function(taxaN){

    # replace sites with zero taxa as NAs
    taxaMat = dimArr[,1:taxaN,]
    dimAnalog = rep(NA,nrow(dimArr))
    noZeros = which(rowSums(taxaMat[,,2])>0)
    
    refSet = taxaMat[rowSums(taxaMat[,,1]) > 0,,1]
    
    dimAnalog[noZeros] = sapply(noZeros, function(n){
      novelty(taxaMat[n,,2] / sum(taxaMat[n,,2]), prop.table(refSet, 1), method="bray")
    })
    
    return(dimAnalog)
    
  })
  
})

# longform for modelling
taxaDimLong = do.call("rbind", 
                lapply(1:nrow(taxaDimNov[[1]]), function(n){
  print(n)
  # pull out site results from each rep
  siteReps = t(sapply(taxaDimNov, function(x){x[n,]}))
  
  siteDf = data.frame(novelty = as.vector(siteReps),
                      noveltyDev = as.vector(siteReps - siteReps[,ncol(siteReps)]),
                      site = dimnames(neoArrSub)[[1]][n],
                      iter = rep(1:nrow(siteReps), ncol(siteReps)),
                      taxaN = rep(1:ncol(siteReps), each = nrow(siteReps)))
                        
  return(siteDf)
  }))

# mixed LM
dimM = lmer(noveltyDev ~ log(taxaN) + (1|iter) + (1|site), data=taxaDimLong)
summary(dimM)

pdf("./plots/taxaDimTest.pdf", height=5, width=7.5, useDingbats=FALSE)
par(mar=c(3,3.5,0.5,0.5), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))
boxplot(tapply(taxaDimLong$noveltyDev, 
               list(taxaDimLong$site,
                    taxaDimLong$taxaN), 
               median, na.rm=TRUE), xaxt="n", col="grey80")
abline(h=0, lty="31")
axis(side=1, mgp=c(3,0.1,0), at=1:19, labels=2:20); mtext(side=1, line=1.25, text="Number of dimensions (taxa)")
mtext(side=2, line=2.5, las=0, text="Deviation from known no-analog novelty")

lines(y=predict(dimM, newdata=data.frame(taxaN=2:20), re.form=NA), x=1:19, col="red", lwd=3)

dev.off()

# GENUS-LEVEL NOVELTY ####

# create time-bin x grid x taxa array
neoArrGen <- tapply(neoSub$value,
                 list(neoSub$siteid, neoSub$genus, neoSub$timeBin), sum, na.rm=TRUE)
neoArrGen[is.na(neoArrGen)] = 0

# remove genera aquatic families
aquGen = unique(neoSub$genus[neoSub$family %in% c("Potamogetonaceae", "Nymphaeaceae", 
                                                "Typhaceae", "Cabombaceae", 
                                                "Alismataceae", "Haloragaceae")])

neoArrGen <- neoArrGen[,!dimnames(neoArrGen)[[2]] %in% aquGen,]

# remove extremely rare taxa - 50 most abundant genera
neoArrGenProp = array(apply(neoArrGen, 3, function(x){prop.table(x, margin=1)}),
                   dim=dim(neoArrGen), dimnames=dimnames(neoArrGen))

genMean = apply(neoArrGenProp, 2, mean, na.rm=TRUE)
neoArrGen = neoArrGen[,genMean >= sort(genMean, decreasing=TRUE)[50],]

#               Calculate novelty ####

# now we can make novelty assessments across modern vs pre-modern time slices

# subset array to time series with slices in the modern (-50-450 ybp) and pre-modern
# (450-950 ybp)

neoArrGenSub = neoArrGen[rowSums(neoArrGen[,,1] > 0) & rowSums(neoArrGen[,,2] >0),,1:2]
neoArrGenRoot = sqrt(neoArrGenSub)

novGridGen <- data.frame(siteid = rownames(neoArrGenSub))
novGridGen <- merge(novGridGen, subSite[,c("siteid", "long", "lat")],
                 by.x="siteid", by.y="siteid", all.x=TRUE, all.y=FALSE, sort=FALSE)

# novelty calculated four ways:


# how the past differs from the present (no-analog)
novGridGen = cbind(novGridGen, 
                do.call('rbind', lapply(1:dim(neoArrGenRoot)[1], function(n){
                  x = novelty(neoArrGenRoot[n,,2] / sum(neoArrGenRoot[n,,2]), 
                              prop.table(neoArrGenRoot[,,1], 1), method="bray")
                  
                  return(data.frame(noAnalog = x,
                                    noAnalogWhich = names(x)))
                })))

# how the present differs from the past (time's arrow novelty)
novGridGen = cbind(novGridGen, 
                do.call("rbind", lapply(1:dim(neoArrGenRoot)[1], function(n){
                  x = novelty(neoArrGenRoot[n,,1] / sum(neoArrGenRoot[n,,1]), 
                              prop.table(neoArrGenRoot[,,2], 1), method="bray")
                  return(data.frame(timeArrow = x,
                                    timeArrowWhich = names(x)))
                })))

# how the past differs from the past (past comparison)
novGridGen = cbind(novGridGen,
                do.call("rbind", lapply(1:dim(neoArrGenRoot)[1], function(n){
                  x = novelty(neoArrGenRoot[n,,2] / sum(neoArrGenRoot[n,,2]), prop.table(neoArrGenRoot[,,2], 1), method="bray", nSize=2)[2]
                  return(data.frame(pastComp = x,
                                    pastCompWhich = names(x)))
                })))

# how the present differs from the present (present comparison)
novGridGen = cbind(novGridGen,
                do.call("rbind", lapply(1:dim(neoArrGenRoot)[1], function(n){
                  x = novelty(neoArrGenRoot[n,,1] / sum(neoArrGenRoot[n,,1]), prop.table(neoArrGenRoot[,,1], 1), method="bray", nSize=2)[2]
                  return(data.frame(presentComp = x,
                                    presentCompWhich = names(x)))
                })))

# re-order columns so novelty columns are all together
novGridGen = novGridGen[,c("siteid", "long", "lat", 
                     "noAnalog", "timeArrow", "pastComp", "presentComp",
                     "noAnalogWhich", "timeArrowWhich", "pastCompWhich", "presentCompWhich")]

# does sample size correlate with novelty measures?
cor(sapply(novGridGen[,4:7], logit), 
    use="complete.obs")

#         novelty PCA ####

novGridGen = novGridGen[order(rowSums(novGridGen[,4:7])),]
a = princomp(sapply(novGridGen[,4:7], logit))
aVar = a$sdev / sum(a$sdev)
#a = princomp(novGridGen[,4:7]) # does not logit-transform dissimilarities

# highlight points with at least one novelty value > 0.2
novelCat = rowSums(novGridGen[,4:7] > 0.2) > 0

novGridGen$fig3Number = ""
novGridGen$fig3Number[novelCat] = sum(novelCat):1
write.csv(novGridGen, "./outputs/novelResultsGenus.csv")

pdf("./plots/noveltyPCAGenus.pdf", height=3.5, width=7, useDingbats = FALSE)
split.screen(rbind(c(0.1,0.49,0.15,0.99),
                   c(0.6,0.99,0.15,0.99)))
vectorScale = 1

col1 = rgb(colorRamp(c("white", "grey85", "black"), bias=2)(unitScale(a$scores[,1]))/255)
col2 = rgb(colorRamp(c("purple1", "orange"), bias=2)(unitScale(a$scores[,2]))/255)
col3 = rgb(colorRamp(c("red", "white", "blue"))(unitScale(a$scores[,3], custMin= -max(abs(a$scores[,3])), custMax=max(abs(a$scores[,3]))))/255)
col4 = rgb(colorRamp(c("yellow", "white", "darkgreen"))(unitScale(a$scores[,4], custMin= -max(abs(a$scores[,4])), custMax=max(abs(a$scores[,4]))))/255)

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(a$scores[,1:2], asp=1, type="n", axes=FALSE, xlab="", ylab="")
axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1.25, text=paste0("PC1 (", sprintf("%.2f", aVar[1]*100), "%)"), las=0)
axis(side=2); mtext(side=2, line=2, text=paste0("PC2 (", sprintf("%.2f", aVar[2]*100), "%)"), las=0)
Arrows(x0=0, y0=0, col="grey50",
       x1=a$loadings[,1]*vectorScale, 
       y1=a$loadings[,2]*vectorScale, 
       arr.type="triangle", arr.length=0.15, arr.width=0.15, lwd=2)
text(x=a$loadings[,1]*vectorScale, y=a$loadings[,2]*vectorScale, 
     labels=c("No-analog", "Time's arrow", "Past comp", "Present comp"),
     adj=0, col="grey50", pos=4)

points(a$scores[,1:2], pch=c(21,22)[as.factor(novelCat)], lwd=0.5,
       bg=mapply(c1=col1, c2=col2, colorMix), cex=ifelse(novelCat, 1.5, 1))
text(a$scores[novelCat,1:2], labels=sum(novelCat==TRUE):1, cex=0.575, col="white")

text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
     labels="(A)", font=1, adj=0)
box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(a$scores[,3:4], asp=1, type="n", axes=FALSE, xlab="", ylab="", ylim=c(-0.8,0.8))
axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1.25, text=paste0("PC3 (", sprintf("%.2f", aVar[3]*100), "%)"), las=0)
axis(side=2); mtext(side=2, line=2, text=paste0("PC4 (", sprintf("%.2f", aVar[4]*100), "%)"), las=0)
Arrows(x0=0, y0=0, col="grey50",
       x1=a$loadings[,3]*vectorScale, 
       y1=a$loadings[,4]*vectorScale, 
       arr.type="triangle", arr.length=0.15, arr.width=0.15, lwd=2)
text(x=a$loadings[,3]*vectorScale, y=a$loadings[,4]*vectorScale, 
     labels=c("No-analog", "Time's arrow", "Past comp", "Present comp"),
     adj=0, col="grey50", pos=c(3,1,2,3))
points(a$scores[,3:4], pch=c(21,22)[as.factor(novelCat)],  lwd=0.5,
       bg=mapply(c1=col3, c2=col4, colorMix), cex=ifelse(novelCat, 1.5, 1))
text(a$scores[novelCat,3:4], labels=sum(novelCat==TRUE):1, cex=0.575, 
     col="black")
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
     labels="(B)", font=1, adj=0)
box()
close.screen(2)
dev.off()

#         study map with novelty overlaid ####

pdf("./plots/studyMapNoveltyGenus.pdf", height=8.35, width=6, useDingbats = FALSE)
split.screen(rbind(c(0.1,0.545,0.52,0.99),
                   c(0.545,0.99,0.52,0.99),
                   c(0.1,0.545,0.05,0.52),
                   c(0.545,0.99,0.05,0.52),
                   c(0.375,0.545,0.47,0.55),
                   c(0.12,0.32,0.55,0.58),
                   c(0.565,0.765,0.55,0.58),
                   c(0.12,0.32,0.08,0.11),
                   c(0.565,0.765,0.08,0.11)))

vectorScale = 1
novMat = novGridGen[,c("noAnalog", "timeArrow", "pastComp", "presentComp")]

sapply(1:4, function(n){
  
  screen(n)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(NULL, xlim=c(-97.5,-90), ylim=c(41,49.75), xlab="", ylab="", axes=FALSE)
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
       col="grey90")
  
  subStates = st_crop(usStates, xmin=-102, xmax=-85, ymin=35, ymax=55)
  plot(subStates, max.plot=1, col="white", bg="grey80", border="black", add=TRUE)
  
  novCol = c("darkblue", "red", "darkgreen", "purple")[n]
  novRamp = colorRampPalette(c("white", novCol))
  points(novGridGen$lat ~ novGridGen$long, 
         pch=c(21,22)[as.factor(novMat[,n] > 0.2)], 
         lwd=0.5, cex=1.25,
         bg=novRamp(5)[cut(ifelse(novMat[,n] > 0.25, 0.25, novMat[,n]), breaks=seq(0,0.25,0.05))])
  box()
  
  text(novGridGen$lat ~ novGridGen$long,
       labels=novGridGen$fig3Number, col="white", cex=0.5)
  
  if(n %in% c(3:4)){
    axis(side=1, at=seq(-100,-90,2), labels=parse(text=paste0(abs(rev(seq(90,100,2))), "*degree~W")), mgp=c(3,0.2,0))
    axis(side=3, at=seq(-100,-90,2), labels=NA, tcl=0.25)
  }
  if(n %in% c(1,3)){
    axis(side=2, mgp=c(3,0.5,0), at=seq(40,50,2), labels=parse(text=paste0(seq(40,50,2), "*degree~N")))
    axis(side=4, mgp=c(3,0.5,0), at=seq(40,50,2), labels=NA, tcl=0.25)
  }
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       adj = 0, labels=paste0("(", LETTERS[n], ")"), font=2)
  text(x=relative.axis.point(0.11, "x"),
       y=relative.axis.point(0.95, "y"),
       adj = 0, labels=c("No-analog", "Time's arrow", "Past comparison", "Present comparison")[n])
  rect(xleft=par("usr")[1], xright=relative.axis.point(0.5,"x"), ybottom=par("usr")[3], ytop=relative.axis.point(0.1,"y"),
       col="white", border=NA)
  
  close.screen(n)
  
  
  screen(5+n)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.1,0), las=1)
  plot.new()
  gridX = seq(par('usr')[1], par("usr")[2], len=6)
  rect(xleft=gridX[-length(gridX)], xright=gridX[-1],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=novRamp(5))
  axis(side=1, at=gridX[c(1,3,5,7)], labels=seq(0,0.3,0.1))
  box()
  close.screen(5+n)
  
})

screen(5)
par(mar=c(0,0,0,0), mgp=c(3,0,0), tcl=-0.25)
plot(worldmap, max.plot=1, xlim=c(-125,-65), ylim=c(35,45), col="grey80", border="grey80", bg="white")
plot(usOutline, max.plot=1, add=TRUE, col="grey50")
rect(xleft=-97.5, xright=-90, ybottom=41, ytop=49, border="red", lwd=2)
box()
close.screen(5)
close.screen(all.screens=TRUE)
dev.off()

#         nMDS plot ####

# run ordination on relativized counts
neoArrGenSubProp = array(apply(neoArrGenSub, 3, function(x){prop.table(x, margin=1)}),
                      dim=dim(neoArrGenSub), dimnames=dimnames(neoArrGenSub))

# collapse 3D array to site-taxa matrix
neoNMDSArr = rbind(neoArrGenSubProp[,,1], # modern
                   neoArrGenSubProp[,,2]) # pre-modern

set.seed(001155)
neoOrd <- metaMDS(neoNMDSArr)

pdf("./plots/newOrdGenus.pdf", height=8, width=8, useDingbats = FALSE)
par(mfrow=c(2,2), oma=c(3,3,0.5,0.5), mar=c(0,0,0,0), mgp=c(3,0.5,0), ps=10, tcl=-0.25, las=1)

sapply(1:4, function(n){
  
  plot(neoOrd$points, type="n", asp=1, xlab="", ylab="", xaxt="n", yaxt="n",
       xlim=c(-1, 1.2))
  
  if(n %in% c(3:4)){
    axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1, text="nMDS1")
  } else {axis(side=1, labels=NA)}
  
  if(n %in% c(1,3)){
    axis(side=2); mtext(side=2, line=2, las=0, text="nMDS2")
  } else {axis(side=2, labels=NA)}
  
  novCol = c("darkblue", "red", "darkgreen", "purple")[n]
  novRamp = colorRampPalette(c("white", novCol))(6)
  
  points(neoOrd$species, pch=16, lwd=2, col="grey")
  text(neoOrd$species, labels=gsub("aceae|inosae|ositae", "", rownames(neoOrd$species)), pos=2, col="grey")
  
  presentPoints = neoOrd$points[1:dim(neoArrGenSubProp)[1],]
  pastPoints = neoOrd$points[(dim(neoArrGenSubProp)[1]+1):nrow(neoOrd$points),]
  novGridOrder = novGridGen[match(rownames(neoArrGenSubProp), novGridGen$siteid),]
  
  novScaling = function(x){c(0.4,0.6,0.9,1.2,1.5,2)[cut(x, breaks=c(0,0.05,0.1,0.15,0.2,0.25,1))]}
  
  # no-analog points
  if(n==1){
    
    segments(x0 = pastPoints[,1],
             y0 = pastPoints[,2],
             x1 = presentPoints[match(novGridOrder$noAnalogWhich, rownames(pastPoints)), 1],
             y1 = presentPoints[match(novGridOrder$noAnalogWhich, rownames(pastPoints)), 2])
    
    points(presentPoints, pch=4, col="black", cex=0.5)
    
    points(pastPoints, 
           pch=21, cex=novScaling(novGridOrder$noAnalog), 
           bg=novRamp[cut(novGridOrder$noAnalog, breaks=c(0,0.05,0.1,0.15,0.2,0.25,1))])
    
  }
  
  # time's arrow points
  if(n==2){
    
    segments(x0= presentPoints[,1],
             y0= presentPoints[,2],
             x1 = pastPoints[match(novGridOrder$timeArrowWhich, rownames(pastPoints)), 1],
             y1 = pastPoints[match(novGridOrder$timeArrowWhich, rownames(pastPoints)), 2],
             lty="31")
    
    points(pastPoints, pch=4, col="black", cex=0.5)
    
    
    points(presentPoints, 
           pch=21, cex=novScaling(novGridOrder$timeArrow), 
           bg=novRamp[cut(novGridOrder$timeArrow, breaks=c(0,0.05,0.1,0.15,0.2,0.25,1))])
    
  }
  
  # past comparison points
  if(n==3){
    
    segments(x0 = pastPoints[,1],
             y0 = pastPoints[,2],
             x1 = pastPoints[match(novGridOrder$pastCompWhich, rownames(pastPoints)), 1],
             y1 = pastPoints[match(novGridOrder$pastCompWhich, rownames(pastPoints)), 2],
             lty="31")
    
    points(pastPoints, 
           pch=21, cex=novScaling(novGridOrder$pastComp), 
           bg=novRamp[cut(novGridOrder$pastComp, breaks=c(0,0.05,0.1,0.15,0.2,0.25,1))])
    
    
    legend(x=relative.axis.point(0.02, "x"),
           y=relative.axis.point(0.9, "y"),
           pch=c(rep(21,4),4), pt.cex=novScaling(seq(0.025, 0.175, 0.05)),
           pt.bg=c(colorRampPalette(c("white", "grey20"))(5), NA),
           legend = c("< 0.05", "0.05 - 0.10","0.10 - 0.15","0.15 - 0.20", "> 0.20", "Ref set"),
           title = "Novelty")
    
  }
  
  if(n==4){
    
    segments(x0 = presentPoints[,1],
             y0 = presentPoints[,2],
             x1 = presentPoints[match(novGridOrder$presentCompWhich, rownames(presentPoints)), 1],
             y1 = presentPoints[match(novGridOrder$presentCompWhich, rownames(presentPoints)), 2],
             lty="31")
    
    points(presentPoints, 
           pch=21, cex=novScaling(novGridOrder$presentComp), 
           bg=novRamp[cut(novGridOrder$presentComp, breaks=c(0,0.05,0.1,0.15,0.2,0.25,1))])
    
  }
  
  text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
       labels=paste0("(",LETTERS[n],")"), font=2, adj=0)
  text(x=relative.axis.point(0.07, "x"), y=relative.axis.point(0.95, "y"),
       labels=c("No-analog", "Time's arrow", "Past comparison", "Present comparison")[n], 
       adj=0)
  
  text(x=relative.axis.point(0.035, "x"), y=relative.axis.point(0.035, "y"),
       labels=paste0("Stress = ", sprintf("%.3f", neoOrd$stress)), adj=0)
  
})
dev.off()

# GENUS VS FAMILY NOVELTY PLOT ####

# would be good to know what proportion of pollen in site was
# ID-able to genus
topSub = neoSub[neoSub$timeBin <= 700 & !is.na(neoSub$timeBin), ]
genProp = sapply(split(topSub, f=topSub$siteid), 
                 function(x){
                   sum(x$value[!is.na(x$genus)]) / sum(x$value)
                 })
genProp = data.frame(siteid = names(genProp),
                     genProp = genProp)

novGrid = merge(novGrid, genProp, by.x="siteid", by.y="siteid",
                all.x=TRUE, all.y=FALSE, sort=FALSE)

genRamp = colorRamp(rev(c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF")))

familyNov = novGrid[,4:7]
genusNov = novGridGen[,4:7]

pdf("./plots/genFamComp.pdf", height=8, width=8)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.1,0.545,0.545,0.99),
                   c(0.545,0.99,0.545,0.99),
                   c(0.1,0.545,0.1,0.545),
                   c(0.545,0.99,0.1,0.545),
                   c(0.775,0.975,0.62,0.65)))

sapply(1:ncol(familyNov), function(n){
  
  screen(n)
  plot(genusNov[,n] ~ familyNov[,n], type="n", axes=FALSE, xlab="", ylab="",
       xlim=c(0,0.23), ylim=c(0,0.3))
  
  if(n %in% c(3,4)){axis(side=1, mgp=c(3,0.1,0)); mtext(side=1, line=1.25, text="Family novelty")
  } else {axis(side=1, labels=NA)}
  if(n %in% c(1,3)){axis(side=2); mtext(side=2, line=2, las=0, text="Genus novelty")
  } else {axis(side=2, labels=NA)}
  
  abline(a=0,b=1,lty="31",col="grey")
  
  novCor = cor(cbind(familyNov[,n], genusNov[,n]))[1,2]
  
  abline(h=0.2,col="red")
  segments(x0=0.15, x1=0.15 , y0=ifelse(n==2, relative.axis.point(0.2, "y"),par("usr")[3]), 
           y1=relative.axis.point(0.9,"y"), col="red")
  
  points(genusNov[,n] ~ familyNov[,n], pch=21, 
         bg=rgb(genRamp(novGrid$genProp)/255))
  text(x=relative.axis.point(0.02, "x"), y = relative.axis.point(0.95, "y"), 
       labels=paste0("(",LETTERS[n],")"), font=2, adj=0)
  text(x=relative.axis.point(0.1, "x"), y = relative.axis.point(0.95, "y"), 
       labels=paste0(c("No-analog", "Time's arrow", "Past comparison", "Present comparison")[n],
                     " (R = ", sprintf("%.3f", novCor), ")"), adj=0)
  box()
  close.screen(n)
})

screen(5)
par(mgp=c(3,0,0), tcl=-0.25)
image(x=seq(0,1,len=200), y=c(0,1),
      z=matrix(1:200, ncol=1), col=rgb(genRamp(seq(0,1,len=200))/255),
      axes=FALSE, xlab="", ylab="", useRaster=TRUE)
axis(side=1)
mtext(side=1, text="Pollen fraction identifiable\nto genus", line=1.75)
box()
close.screen(5)
close.screen(all.screens = TRUE)
dev.off()

# Summary stats ####

novGrid$coordID = paste0(novGrid$long,":",novGrid$lat)
coordDists <- expand.grid(novGrid$coordID, novGrid$coordID)
distMat <- gc.dist(novGrid$lat[match(coordDists$Var1, novGrid$coordID)],
                   novGrid$long[match(coordDists$Var1, novGrid$coordID)],
                   novGrid$lat[match(coordDists$Var2, novGrid$coordID)],
                   novGrid$long[match(coordDists$Var2, novGrid$coordID)])
summary(tapply(distMat, coordDists[,1], function(x){max(x)}))

# rank families by total % count and then export
famArr = neoArr[,,c(1,2)]
famArr = prop.table(famArr, 1)

famTable = apply(famArr, 2, mean, na.rm=TRUE)
famTable = data.frame(family = names(famTable),
                      meanProp = famTable)
famTable = famTable[order(famTable$meanProp, decreasing=TRUE),]

# also merge in the ordination coords (why not)
ordfam = data.frame(family = rownames(neoOrd$species),
                    MDS1 = neoOrd$species[,1],
                    MDS2 = neoOrd$species[,2])

famTable = merge(famTable, ordfam, by.x="family", by.y="family",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

famTable$meanProp = sprintf("%.2f", famTable$meanProp * 100)
famTable = famTable[,c("family", "meanProp", "MDS1", "MDS2")]
rownames(famTable) = NULL
write.csv(famTable, "./outputs/familyTable.csv")


genArr = neoArrGen[,,c(1,2)]
genArr = prop.table(genArr, 1)

genTable = apply(genArr, 2, mean, na.rm=TRUE)
genTable = data.frame(genus = names(genTable),
                      meanProp = genTable)
# merge in family names
genTable = merge(genTable, neoSub[!duplicated(neoSub$genus),c("genus", "family")],
                 by.x="genus", by.y="genus")

# also merge in the ordination coords (why not)
ordGen = data.frame(family = rownames(neoOrd$species),
                    MDS1 = neoOrd$species[,1],
                    MDS2 = neoOrd$species[,2])

genTable = merge(genTable, ordGen, by.x="genus", by.y="family",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

genTable = genTable[order(genTable$meanProp, decreasing=TRUE),]

genTable$meanProp = sprintf("%.2f", genTable$meanProp * 100)
genTable = genTable[,c("family","genus", "meanProp", "MDS1", "MDS2")]
rownames(genTable) = NULL

write.csv(genTable, "./outputs/genusTable.csv")
