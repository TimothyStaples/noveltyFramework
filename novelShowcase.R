rm(list=ls())

library(tidyr)
library(vegan)
library(data.table)
library(rworldmap)
library(maptools)
library(abind)
library(viridisLite)
library(hilldiv)
library(VennDiagram)
library(DHARMa)
library(betareg)
library(raster)
library(performance)
library(shape)
library(prettymapr)

setwd("/Users/uqtstapl/Library/CloudStorage/Dropbox/Tim/Post-doc/Research projects/PaleoNovelty/Manifesto/GEB/novelFramework")

sapply(list.files("./functions", full.names = TRUE), source)
cutCenter <- function(x){
  rowMeans(cbind(as.numeric( sub("\\((.+),.*", "\\1", x)),
                 as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x) )))
}
unitScale <- function(x){
  (x - min(x)) / (max(x)-min(x))
}

neoSub = read.csv("./rawdata/neoComp.csv")
neoSub = neoSub[,!colnames(neoSub) %in% c("X.2", "X.1", "X", "site.1")]

subSite = read.csv("./rawdata/neoSite.csv")
siteCoords = subSite[,c("site", "long", "lat")]
subSite = cbind(subSite[,!colnames(subSite) %in% colnames(neoSub)], site = subSite$site)

#               Sub-sampling ####

neoSub <- merge(neoSub, subSite,
             by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE)

neoSub$timeBin <- cutCenter(cut(neoSub$age, breaks=seq(300,2300,500)))

# create time-bin x grid x taxa array
neoArr <- tapply(neoSub$count,
                 list(neoSub$site, neoSub$family, neoSub$timeBin), sum, na.rm=TRUE)
neoArr[is.na(neoArr)] = 0

# remove aquatic families
neoArr <- neoArr[,!dimnames(neoArr)[[2]] %in% c("Potamogetonaceae", "Nymphaeaceae", "Typhaceae", "Cabombaceae", "Alismataceae",
                                                "Haloragaceae"),]

# remove extremely rare taxa
taxaMean = apply(neoArr, 2, mean)

neoArr = neoArr[,taxaMean >= sort(taxaMean, decreasing=TRUE)[20],]

# get an appropriate subsampling size
binSize <- round(apply(neoArr, c(1,3), sum))
summary(as.vector(binSize)>100)

# rand pollen draws
iters = 999
sampN = 100
randLevs <- colnames(neoArr)

# create random draws of pollen
neoRand <- apply(neoArr, c(1,3), function(x){
  
  if(sum(x)==0){return(matrix(0, nrow=iters, ncol=length(randLevs),
                       dimnames=list(NULL, randLevs)))}
  
  x <- x[x>0]
  xLevs <- names(x)
  longx <- unlist(sapply(1:length(x), function(n){rep(names(x)[n], x[n])}))
  
  # sample with replacement under threshold
  if(sum(x) < sampN){sampReplace = TRUE} else {sampReplace=FALSE}
  
  reps <- replicate(iters, sample(1:length(longx), sampN, replace=sampReplace))
  repMat <- matrix(longx[reps], nrow=sampN)
  
  repArr <- t(apply(repMat, 2, function(x){table(factor(x, levels=randLevs))}))
  
}, simplify=FALSE)

# this was removing sites (rows) without complete coverage. It's now going to fill empties
# with zeros
#neoRand <- neoRand[apply(neoRand, 1, function(x){sum(sapply(x, is.null))}) == 0,]

neoRandArr <- abind(apply(neoRand, 1, function(x){
  a <- abind(x, along=3)
  return(a)
}, simplify=FALSE), along=4)

# dims are: reps, taxa, time, site

#               Calculate novelty ####

# now we can make novelty assessments across modern vs pre-modern time slices
# for each random sample, look for similarities

targArr <- neoRandArr[,,1,]
refArr <- neoRandArr[,,-1,]

refArr <- apply(refArr, c(1,3,4), function(x){x})

# now novelty is across all space for the past, so aggregate ref array into matrix
# for each random iteration
neoNov <- lapply(1:dim(targArr)[1], function(n){
  print(n)
  targSub <- t(targArr[n,,])
  refSub <- refArr[,n,,]
  refSub <- do.call('rbind', lapply(1:dim(refSub)[2], function(n1){
    x <- refSub[,n1,]
    x <- t(x)
    rownames(x) = paste0(rownames(x), ":",dimnames(refSub)[[2]][n]) 
    return(x)
  }))
  
  # remove empty rows
  targSub <- targSub[rowSums(targSub) > 0, ]
  refSub = refSub[rowSums(refSub) > 0, ]
  
  tempNov <- novelty(sqrt(targSub), sqrt(refSub), method="bray", nSize=1)
  return(list(tempNov, rownames(targSub)))
})

neoNovArr = t(do.call("rbind", lapply(neoNov, function(x){x[[1]]})))

novGrid <- data.frame(id = rownames(neoNovArr))
novGrid <- merge(novGrid, siteCoords,
                 by.x="id", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE)

novGrid$meanNov <- apply(neoNovArr, 1, median)
novGrid$stdNov <- apply(neoNovArr, 1, sd)

# richness - species count
neoH0 <- t(apply(targArr, c(1,3), function(x){hill_div(x, 0)}))
neoH0 = neoH0[match(novGrid$id,rownames(neoH0)),]
novGrid$meanH0 <- apply(neoH0, 1, median)
novGrid$stdH0 <- apply(neoH0, 1, sd)

neoH1 <- t(apply(targArr, c(1,3), function(x){hill_div(x, 1)}))
neoH1 = neoH1[match(novGrid$id,rownames(neoH1)),]
novGrid$meanH1 <- apply(neoH1, 1, median)
novGrid$stdH1 <- apply(neoH1, 1, sd)

# target - target novelty
targNov <- sapply(1:dim(targArr)[1], function(n){
  print(n)
  targSub <- t(targArr[n,,])
  
  targSub <- targSub[rowSums(targSub) > 0, ]
  
  novelty(sqrt(targSub), 
                   sqrt(targSub), method="bray", nSize=2)[,-1]
  
})
rownames(targNov) == novGrid$id
novGrid$meanNovT = apply(targNov, 1, median)
novGrid$stdNovT = apply(targNov, 1, sd)

#               calculate position in distribution ####

# outliers
meanH0M <- glm(meanH0 ~ 1, data=novGrid, family=poisson)
plot(simulateResiduals(meanH0M))
H0res = residuals(meanH0M)
novGrid$H0q = pnorm(H0res, mean=mean(H0res), sd=sd(H0res))
novGrid$meanH0Out = residuals(meanH0M) > quantile(residuals(meanH0M), 0.75)
plot(novGrid$meanH0 ~ novGrid$H0q, col=ifelse(novGrid$meanH0Out, "red", "black"))

meanH1M <- glm(meanH1 ~ 1, data=novGrid, family=Gamma)
plot(simulateResiduals(meanH1M))
H1res = residuals(meanH1M)
novGrid$H1q = pnorm(H1res, mean=mean(H1res), sd=sd(H1res))
novGrid$meanH1Out = residuals(meanH1M) > quantile(residuals(meanH1M), 0.75)
plot(novGrid$meanH1 ~ novGrid$H1q, col=ifelse(novGrid$meanH1Out, "red", "black"))

meanNovM <- betareg(meanNov ~ 1, data=novGrid)
plot(simulateResidualsBeta(meanNovM))
novRes = residuals(meanNovM)
novGrid$Novq = pnorm(novRes, mean=mean(novRes), sd=sd(novRes))
novGrid$meanNovOut = residuals(meanNovM) > quantile(residuals(meanNovM), 0.75)
plot(novGrid$meanNov ~ novGrid$Novq, col=ifelse(novGrid$meanNovOut, "red", "black"))

meanNovTM <- betareg(meanNovT ~ 1, data=novGrid)
plot(simulateResidualsBeta(meanNovTM))
novTRes = residuals(meanNovTM)
novGrid$NovTq = pnorm(novTRes, mean=mean(novTRes), sd=sd(novTRes))
novGrid$meanNovTOut = residuals(meanNovTM) > quantile(residuals(meanNovTM), 0.9)
plot(novGrid$meanNovT ~ novGrid$NovTq, col=ifelse(novGrid$meanNovTOut, "red", "black"))

# summary stats
novCor <- cor(novGrid[,c("meanH0", "meanH1", "meanNov", "meanNovT")])
novCor <- novCor[upper.tri(novCor)]
novOut <- apply(t(combn(c("meanH0Out", "meanH1Out", "meanNovOut", "meanNovTOut"), 2)), 1, 
                function(comb){
                  print(comb)
                  a <- rowSums(novGrid[,comb])
                  sum(a==2) / sum(a>0)
                })
write.csv(data.frame(cor = novCor, out = novOut), "./outputs/metricCorrelation.csv")

# New richness vs Novelty only plot ####

novCor <- as.vector(novCor[upper.tri(novCor)])
novGrid$colCount <- rowSums(novGrid[, c("meanH0Out", "meanNovTOut")])
circRad=0.25

# colour order is: meanH0Out meanH1Out meanNovOut meanNovTOut
catCols = c("#3F6B17", "#F89238")

novGrid$sum <- rowSums(novGrid[, c("H0q", "H1q", "Novq", "NovTq")])

library(mapsf)
library(sf)
library(plotrix)
library(rworldmap)
worldmap = countriesLow
usOutline = read_sf("./shapefiles/cb_2018_us_nation_20m.shp")
usStates = read_sf("./shapefiles/cb_2018_us_state_5m.shp")
usUrban = read_sf("./shapefiles/tl_2020_us_uac20.shp")
usPA = read_sf("./shapefiles/PADUS3_0Combined_Region3.shp")
usPAMarine = read_sf("./shapefiles/PADUS3_0Marine.shp")
usPA1 = read_sf("./shapefiles/PADUS3_0Combined_Region4.shp")
usPA2 = read_sf("./shapefiles/PADUS3_0Combined_Region5.shp")
usPA = st_transform(usPA, crs(usStates))
usPAm <- st_transform(usPAMarine, crs(usStates))
usPA1 = st_transform(usPA1, crs(usStates))
usPA2 = st_transform(usPA2, crs(usStates))

# which sample points are in protected areas?
novSf <- st_as_sf(novGrid, coords = c("long","lat"))
novSf = st_set_crs(novSf, crs(usPA)) 

sf_use_s2(FALSE)
novProtected = st_intersection(novSf,usPA)
novProtected1 = st_intersection(novSf,usPA1)
novProtected2 = st_intersection(novSf,usPA2)
novProtectedM = st_intersection(novSf,usPAm)

novGrid$protected = ifelse(novGrid$id %in% novProtected$id |
                           novGrid$id %in% novProtected1$id | 
                           novGrid$id %in% novProtected2$id |
                           novGrid$id %in% novProtectedM$id,
                           TRUE, FALSE)

usBig <- usUrban[usUrban$ALAND20 >= exp(18),]

pdf("./plots/USSubshowcaseDIV.pdf", height=5.25, width=4.75, useDingbats = FALSE)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

split.screen(rbind(c(0.1,0.99,0.05,0.99),
                   c(0.62,0.95,0.95,0.98),
                   c(0.1,0.3,0.85,0.99)))

screen(1)
plot(NULL, xlim=c(-97.5,-90), ylim=c(41,49), xlab="", ylab="", axes=FALSE)
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     col="grey90")

subStates = st_crop(usStates, xmin=par("usr")[1], xmax=par("usr")[2], ymin=par("usr")[3], ymax=par("usr")[4])
subUrban = st_crop(usBig, xmin=par("usr")[1], xmax=par("usr")[2], ymin=par("usr")[3], ymax=par("usr")[4])
#subPA = st_crop(usPA, xmin=par("usr")[1], xmax=par("usr")[2], ymin=par("usr")[3], ymax=par("usr")[4])

mf_base(subStates, max.plot=1, col="white", bg="grey80", border="grey", add=TRUE)
mf_base(subUrban, max.plot=1, col=rgb(0.5,0.5,0.5,0.5), bg="black", lwd=0.75, add=TRUE)

urbanCent = st_centroid(subUrban)
text(st_coordinates(urbanCent), labels=substr(urbanCent$NAME20, 1, regexpr(",|-", urbanCent$NAME20)-1),
     pos=4, cex=0.75)

axis(side=1, at=seq(-100,-90,2), labels=parse(text=paste0(abs(rev(seq(90,100,2))), "*degree~W")), mgp=c(3,0.2,0))
axis(side=2, mgp=c(3,0.5,0),
     at=seq(40,50,2), labels=parse(text=paste0(seq(40,50,2), "*degree~N")))

summary(novGrid$NovTq)

# sort points by low to high novelty
novGrid <- novGrid[order(novGrid$NovTq, decreasing=TRUE),]

# zero colours
points(novGrid$lat ~ novGrid$long, 
       cex=ifelse(novGrid$NovTq == 0, 0.6, 0.6 + novGrid$NovTq * 2), 
       pch=21,
       bg=rgb(colorRamp(viridis(5))(novGrid$NovTq)/255))

# add numbers to top 9 points
with(novGrid[novGrid$NovTq > 0.9,],
     text(lat ~ long, labels="N", cex=0.8))

box()

addnortharrow(pos="topright", scale=0.6)
close.screen(1)

screen(2)
par(mgp=c(3,0,0), tcl=-0.25)
image(x=seq(0,1,len=200), y=c(0,1),
      z=matrix(1:200, ncol=1), col=rgb(colorRamp(viridis(5))(seq(0,1,len=200))/255),
      axes=FALSE, xlab="", ylab="", useRaster=TRUE)
axis(side=1)
mtext(side=1, text="Novelty", line=0.75)
box()
close.screen(2)

screen(3)
par(mgp=c(3,0,0), tcl=-0.25)
plot(worldmap, max.plot=1, xlim=c(-125,-65), ylim=c(35,45), col="grey80", border="grey80", bg="white")
plot(usOutline, max.plot=1, add=TRUE, col="grey50")
rect(xleft=-97.5, xright=-90, ybottom=41, ytop=49, border="red", lwd=2)
box()
close.screen(3)
dev.off()

#               nMDS plot ####

neoTargArr <- neoRandArr[,,1,]
neoTargMat <- t(apply(neoTargArr, c(2,3), function(x){mean(x)}))

neoTargMat <- neoTargMat[rowSums(neoTargMat) > 0, colSums(neoTargMat) >0]

neoOrd <- metaMDS(neoTargMat)

neoFit <- ordisurf(neoOrd ~ meanH0, novGrid[,c("meanH0","meanH1")],
                   plot=FALSE)

ordPoints <- as.data.frame(neoOrd$points)
ordPoints$id <- rownames(ordPoints)

novGrid <- novGrid[match(rownames(ordPoints), novGrid$id),]

ordPoints <- merge(ordPoints, 
                   novGrid[,c("id", "meanH0Out", "meanH1Out", "meanNovOut", "meanNovTOut", "colCount",
                              "H0q", "H1q", "Novq", "NovTq")],
                   by.x="id", by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)

pdf("./plots/showcaseOrd.pdf", height=4.5, width=4.75, useDingbats = FALSE)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.15,0.99,0.1,0.99),
                   c(0.84,0.88,0.125,0.45)))

screen(1)

plot(ordPoints$MDS2 ~ ordPoints$MDS1, type="n", asp=1,
     axes=FALSE, xlab="", ylab="", ylim=c(-0.9,0.65))

axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.1, text="nMDS 1")
axis(side=2)
mtext(side=2, line=2, text="nMDS 2", las=0)

neoSpMag <- sqrt(rowSums(neoOrd$species^2))
targMag = 0
Arrows(x0=0, y0=0,
       x1=neoOrd$species[neoSpMag>targMag,1],
       y1=neoOrd$species[neoSpMag>targMag,2], 
       lwd=1, arr.type="triangle", arr.width=0.15, arr.length=0.15,
       col="red")
text(x=neoOrd$species[neoSpMag>targMag,1], y=neoOrd$species[neoSpMag>targMag,2],
     col="red", labels=gsub("aceae", "", rownames(neoOrd$species)[neoSpMag>targMag]),
     cex=0.5)

ordData <- expand.grid(x1=seq(par("usr")[1], par("usr")[2], len=200),
                       x2=seq(par("usr")[3], par("usr")[4], len=200))
neoPred <- predict(neoFit, newdata=ordData)
contour(x=seq(par("usr")[1], par("usr")[2], len=200), 
        y=seq(par("usr")[3], par("usr")[4], len=200), 
        z=matrix(neoPred, nrow=200, ncol=200), add=TRUE, col="grey40",
        levels=7:12)

ordPoints <- ordPoints[order(ordPoints$NovTq, decreasing=TRUE),]

points(ordPoints$MDS2 ~ ordPoints$MDS1, 
       cex=ifelse(ordPoints$NovTq == 0, 0.6, 0.6 + ordPoints$NovTq * 2), 
       pch=21,
       bg=rgb(colorRamp(viridis(5))(ordPoints$NovTq)/255))

# add numbers to top 9 points
with(ordPoints[ordPoints$NovTq > 0.9,],
     text(MDS2 ~ MDS1, labels="N", cex=0.8))

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.05, "y"),
     adj=0, labels=paste0("Stress = ", round(neoOrd$stress,3)))

box()
close.screen(1)

screen(2)
par(mgp=c(3,0.5,0), tcl=-0.25)
image(y=seq(0,1,len=200), x=c(0,1),
      z=matrix(1:200, nrow=1), col=rgb(colorRamp(viridis(5))(seq(0,1,len=200))/255),
      axes=FALSE, xlab="", ylab="", useRaster=TRUE)
axis(side=4)
mtext(side=4, text="Novelty quantile", line=1.5, las=0)
box()
close.screen(2)

dev.off()


# genus novelty vs family novelty plot ####

genGrid = read.csv("./rawdata/novGridGenus.csv")
genGrid = genGrid[,c("id", "NovTq")]
colnames(genGrid) = c("id", "NovTqGen")
combGrid = merge(novGrid, genGrid, by.x="id", by.y="id",
                 all.x=TRUE, all.y=TRUE, sort=FALSE)

# would be good to know what proportion of pollen in site was
# ID-able to genus
topSub = neoSub[neoSub$timeBin == 550 & !is.na(neoSub$timeBin), ]
genProp = sapply(split(topSub, f=topSub$site), 
       function(x){
         sum(x$count[!is.na(x$genus)]) / sum(x$count)
       })
genProp = data.frame(id = names(genProp),
                     genProp = genProp)

combGrid = merge(combGrid, genProp, by.x="id", by.y="id",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

genRamp = colorRamp(rev(viridis(5, option="plasma")))

pdf("./plots/genFamComp.pdf", height=4.5, width=4.5)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.125,0.99,0.125,0.99),
                   c(0.6,0.95,0.27,0.32)))

screen(1)
plot(combGrid$NovTq ~ combGrid$NovTqGen, type="n", xaxt="n")
axis(side=1, mgp=c(3,0.1,0))
mtext(side=1, line=1.25, text="Family novelty")
mtext(side=2, line=1.75, las=0, text="Genus novelty")
abline(a=0,b=1,lty="31",col="grey")

novCor = cor(combGrid[,c("NovTq", "NovTqGen")])[1,2]

abline(h=0.9,col="red")
segments(x0=0.9, x1=0.9, y0=0.15, y1=par('usr')[4],col="red")

points(combGrid$NovTq ~ combGrid$NovTqGen, pch=21, 
       bg=rgb(genRamp(combGrid$genProp)/255))
text(x=0.1, y=0, labels=paste0("R = ", round(novCor, 3)), adj=0)
close.screen(1)

screen(2)
par(mgp=c(3,0,0), tcl=-0.25)
image(x=seq(0,1,len=200), y=c(0,1),
      z=matrix(1:200, ncol=1), col=rgb(genRamp(seq(0,1,len=200))/255),
      axes=FALSE, xlab="", ylab="", useRaster=TRUE)
axis(side=1)
mtext(side=1, text="Pollen fraction\nidentifiable to genus", line=1.75)
box()
close.screen(2)
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
genArr = neoArr[,,1]
genArr = genArr[rowSums(genArr) > 0 ,]
genArr = prop.table(genArr, 1)

genTable = apply(genArr, 2, mean)
genTable = data.frame(family = names(genTable),
                      meanProp = genTable)
genTable = genTable[order(genTable$meanProp, decreasing=TRUE),]

# also merge in the ordination coords (why not)
ordGen = data.frame(family = rownames(neoOrd$species),
                    MDS1 = neoOrd$species[,1],
                    MDS2 = neoOrd$species[,2])

genTable = merge(genTable, ordGen, by.x="family", by.y="family",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

genTable$meanProp = sprintf("%.2f", genTable$meanProp * 100)
genTable = genTable[,c("family", "meanProp", "MDS1", "MDS2")]
rownames(genTable) = NULL
write.csv(genTable, "./outputs/familyTable.csv")

