#testing 

library(tidyverse)
library(geomorph)

setwd("C:/Users/hrlexd/Dropbox/PlantAndFood (1)/NativeBees/Bee_pictures/Sam_april/Forewing")

#bring in TPS data
bees<-readland.tps('TPS_forewingMay2025_dig.TPS',specID = 'ID')
?readland.tps
bees[,,1]

#bring in links and classification data
species_class <- factor(c("Imat", "Imat", "Imat", "Imat","Imat","Mont","Mont","Mont","Mont","Mont"))
bee_links<-read.table('links_all_landmarks.txt',sep=',',header=F)
as.matrix(bee_links)
?gpagen
#do a generalixed procrustes analysis on TPS
beesgpa<-gpagen(bees)

#visualise points before procrustes
plot(bees)
#visualise after procrustes
plot(beesgpa)
#run PCA analysis
PCA<-gm.prcomp(beesgpa$coords)
#plot with species coluurs
plot(PCA,col=species_class)
#see summary
summary(PCA)

PCA<-gm.prcomp(beesgpa$coords)
bee.pca.plot<-plot(PCA)
#then can do a representative deformation grid from each species from the PCA
picknplot.shape(bee.pca.plot)

#or can do average from each species using this (I think)
#shape change Immatatus compared to average
ref<-mshape(beesgpa$coords)

gpl.mn<-mshape(beesgpa$coords[,,1:5])
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)
#shape change Monticola compared to average
gpl.mn<-mshape(beesgpa$coords[,,6:10])
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#check if there is a relationship between size and shape
gdf<-geomorph.data.frame(beesgpa,species=species_class)
fit.size<-procD.lm(coords~log(Csize), data = gdf, print.progress=F)
summary(fit.size)
anova(fit.size)
#no sig relationship between size and shape
#visiualise
plotAllometry(fit.size,size=gdf$Csize, logsz=TRUE,method='PredLine',col=species_class)
#same thing here:
plot(fit.size,type='regression',reg.type="PredLine", predictor =log(gdf$Csize))
#first principal component scores fitted values - PredLine
#regression scores-standardised regression shapes 
plotAllometry(fit.size,size=gdf$Csize, logsz=TRUE,method='RegScore',col=species_class)
#partial least squares for relationship between shape and size
PLS<-two.b.pls(log(gdf$Csize),gdf$coords,print.progress=F)
PLS
#plot common allometric component
plot(PLS,col=species_class)
#first one should be same as PLS
plotAllometry(fit.size,size=gdf$Csize,logsz = T,method='CAC',col=species_class)
#can do a size shape PCA too
plotAllometry(fit.size,size=gdf$Csize,logsz = T,method='size.shape')
#can save the shape of those with picknplot too
fit.size<-procD.lm(coords~log(Csize), iter=0,data = gdf, print.progress=F)
PA<-plotAllometry(fit.size,size=gdf$Csize,logsz = T,method="PredLine",pch=19)
picknplot.shape(PA)

#look at the influence of shape, size and species
fit.common<-procD.lm(coords~log(Csize)+species,data=gdf)
summary(fit.common)

#relationship between size and species
fit1<-procD.lm(log(gdf$Csize)~species_class)
summary(fit1)
anova(fit1)
#relationship between shape and species
fit2<-procD.lm(gdf$coords~species_class)
summary(fit2)
PW<-pairwise(fit2,groups=species_class)
summary(PW)

#dos shape differ between species, while accounting for shape co varying with size
fit3<-procD.lm(gdf$coords ~log(gdf$Csize)*species_class)
summary(fit3)
PW<-pairwise(fit3,groups=species_class)
summary(PW)

#can do one where it would be
fit.common<-procD.lm(coords~log(Csize)+sex/species,data=gdf)
fit.null<-procD.lm(coords~log(Csize)+species,data=gdf)
#to see pairwise differences between means fo sex
PW<-pairwise(fit.common,fit.null,groups=Sex)
summary(PW,test.type='dist',confidence=0.95)
#to see pairwise differences between means fo species
PW<-pairwise(fit.common,fit.null,groups=Species)
summary(PW,test.type='dist',confidence=0.95)

morphol.disparity(coords~Csize,groups=NULL,data=gdf)
#outputs procrusted variances, pairwise absolute differences and significances
morphol.disparity(coords~1,groups=~species_class,data=gdf)

morphol.disparity(coords~species_class,groups=species_class,data=gdf)



two.d.array(beesgpa$coords)

#removing weird landmark
bees
bees_removefirst <- bees[-1,,] 
bees_removefirst
beesgpa<-gpagen(bees_removefirst)







#plot points manually
#same as plot(beesgpa)
consensus<-apply(beesgpa$coords,c(1,2),mean)
plot(consensus, asp=1,type='n')
for(i in 1:length(beesgpa$coords[,,1]))
    points(beesgpa$coords[,,i])
points(consensus,col='Red',cex=2,pch=20)
ref<-mshape(beesgpa$coords)





gdf<-geomorph.data.frame(beesgpa,species=species_class)
fit.size<-procD.lm(coords~log(Csize), data = gdf, print.progress=F)
summary(fit.size)
plotAllometry(fit.size,size=gdf$Csize,logsz=TRUE,method='PredLine',col=species_class)
morphol.disparity(coords~1,groups=species_class,data=gdf)

fit1<-procD.lm(gdf$Csize~species_class)
summary(fit1)
fit2<-procD.lm(gdf$coords~species_class)
summary(fit2)


gdf$species

data("plethodon")
attrspeciesattributes(plethodon)
plethodon$links
