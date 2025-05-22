#testing 

library(tidyverse)
library(geomorph)

setwd("C:/Users/hrlexd/Dropbox/PlantAndFood (1)/NativeBees/Bee_pictures/Sam_april/Forewing")

#bring in TPS data
#bees<-readland.tps('TPS_forewingMay2025_dig.TPS',specID = 'imageID')
bees<-readland.tps('../../Forewing_file_tps_dig_15_05_2025_1.TPS',specID = 'imageID')

?readland.tps
bees[,,1]
bees
bees[1,,]

typeof(bees)
Sample_names<-labels(bees)[[3]]
Sample_names
Sample_names<-gsub('.*\\\\','',Sample_names)
species_list<-gsub('_.*','',Sample_names)
species_list
unique(species_list)
species_list<-gsub('L.Kanap','L.kanap',species_list)

#bring in links and classification data
#species_class <- factor(c("Imat", "Imat", "Imat", "Imat","Imat","Mont","Mont","Mont","Mont","Mont"))
bee_links<-read.table('links_all_landmarks.txt',sep=',',header=F)
bee_links<-read.table('links_all_landmarks_remove1_sam.txt',sep=',',header=F)
as.matrix(bee_links)
bees
bees_removefirst <- bees[-1,,] 
bees_removefirst

#do a generalixed procrustes analysis on TPS
beesgpa<-gpagen(bees)
beesgpa<-gpagen(bees_removefirst)


#beesgpa<-gpagen(bees)

#visualise points before procrustes
plot(bees)
#visualise after procrustes
plot(beesgpa)
#quite a bit of variation in that end point

#run PCA analysis
PCA<-gm.prcomp(beesgpa$coords)
#plot with species coluurs
?plot
plot(PCA,pch = 19)
#plot(PCA,col=species_class, pch = 19)
plot(PCA,col=as.factor(species_list), pch = 19)
legend('bottomright',legend=unique(as.factor(species_list)),pch=19,bty='n',col=1:9)
PCA$shapes
#ggplot with legend
PCaxis<-as.data.frame(PCA$x[,1:2])
library(RColorBrewer)
ggplot(PCaxis,aes(x=Comp1,y=Comp2,colour = species_list))+
  geom_point(size=4)+
  theme_bw()+
  scale_color_brewer(palette='Paired')+
  labs(colour='Species',x='PC2',y='PC2',title="PCA")


bees_leioproctus <- bees[,,-c(11:15,31:35)] 
beesleio<-gpagen(bees_leioproctus)
Sample_names_leio<-labels(bees_leioproctus)[[3]]
Sample_names_leio
Sample_names_leio<-gsub('.*\\\\','',Sample_names_leio)
species_list_leio<-gsub('_.*','',Sample_names_leio)
species_list_leio
unique(species_list_leio)
species_list_leio<-gsub('L.Kanap','L.kanap',species_list_leio)





#see summary
summary(PCA)

PCA<-gm.prcomp(beesgpa$coords)

#shape change species compared to average
ref<-mshape(beesgpa$coords)
species_list

#nesocoletes (note average is average across all Leioproctus sp. just grouped for veiwing)
#paahua
gpl.mn<-mshape(beesgpa$coords[,,31:35])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#fulv
gpl.mn<-mshape(beesgpa$coords[,,c(11:15)])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#leioproctus leioproctus
#mont
gpl.mn<-mshape(beesgpa$coords[,,2:6])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#bolt
gpl.mn<-mshape(beesgpa$coords[,,c(1,7:10)])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)


#huakiwi
gpl.mn<-mshape(beesgpa$coords[,,16:20])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#imat
gpl.mn<-mshape(beesgpa$coords[,,21:25])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#kanap
gpl.mn<-mshape(beesgpa$coords[,,26:30])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#pango
gpl.mn<-mshape(beesgpa$coords[,,36:40])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#vest
gpl.mn<-mshape(beesgpa$coords[,,41:45])
#plotRefToTarget(M1=ref,M2=gpl.mn,mag=2)
plotRefToTarget(M1=ref,M2=gpl.mn,mag=2,links=bee_links)

#check if there is a relationship between size and shape
gdf<-geomorph.data.frame(beesgpa,species=species_list)
fit.size<-procD.lm(coords~log(Csize), data = gdf, print.progress=F)
summary(fit.size)
anova(fit.size)
#sig relationship between size and shape
#visiualise
plotAllometry(fit.size,size=gdf$Csize, logsz=TRUE,method='PredLine',col=as.factor(species_list),pch=19)
#same thing here:
#plot(fit.size,type='regression',reg.type="PredLine", predictor =log(gdf$Csize))
#first principal component scores fitted values - PredLine
#regression scores-standardised regression shapes 
plotAllometry(fit.size,size=gdf$Csize, logsz=TRUE,method='RegScore',col=as.factor(species_list),pch=19)
#partial least squares for relationship between shape and size
PLS<-two.b.pls(log(gdf$Csize),gdf$coords,print.progress=F)
PLS
#plot common allometric component
#plot(PLS,col=species_class)
#first one should be same as PLS
#plotAllometry(fit.size,size=gdf$Csize,logsz = T,method='CAC',col=species_class)
#can do a size shape PCA too
#plotAllometry(fit.size,size=gdf$Csize,logsz = T,method='size.shape')
#can save the shape of those with picknplot too
#fit.size<-procD.lm(coords~log(Csize), iter=0,data = gdf, print.progress=F)
#PA<-plotAllometry(fit.size,size=gdf$Csize,logsz = T,method="PredLine",pch=19)
#picknplot.shape(PA)

#look at the influence of shape, size and species
fit.common<-procD.lm(coords~log(Csize)+species,data=gdf)
summary(fit.common)

#relationship between size and species
fit1<-procD.lm(log(gdf$Csize)~species_list)
summary(fit1)
anova(fit1)
#relationship between shape and species
fit2<-procD.lm(gdf$coords~species_list)
summary(fit2)
PW<-pairwise(fit2,groups=species_list)
summary(PW)

#dos shape differ between species, while accounting for shape co varying with size
fit3<-procD.lm(gdf$coords ~log(gdf$Csize)*species_list)
summary(fit3)
PW<-pairwise(fit3,groups=species_list)
summary(PW)


morphol.disparity(coords~Csize,groups=NULL,data=gdf)
#outputs procrusted variances, pairwise absolute differences and significances
#test whether different species are significantly different in shape to overall variation of dataset
morphol.disparity(coords~1,groups=~species_list,data=gdf)
#tests whether or which species are significantly different to each other
morphol.disparity(coords~species_list,groups=~species_list,data=gdf)

morphol.disparity(coords~log(Csize)+species_list,groups=species_list,data=gdf)

morphol.disparity(coords~log(Csize),groups=~species_list,data=gdf)

#when we do some males we can look at sex/species

#can do one where it would be
#fit.common<-procD.lm(coords~log(Csize)+sex/species,data=gdf)
#fit.null<-procD.lm(coords~log(Csize)+species,data=gdf)
#to see pairwise differences between means fo sex
#PW<-pairwise(fit.common,fit.null,groups=Sex)
#summary(PW,test.type='dist',confidence=0.95)
#to see pairwise differences between means fo species
#PW<-pairwise(fit.common,fit.null,groups=Species)
#summary(PW,test.type='dist',confidence=0.95)




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
