########################################################
########################################################
################ Weed paper two: section 3 #############

match<-read.table("w2s3.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
taxo.g<-read.table("labtaxa.guilds_weed.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
ngs.r<-read.table("otu.neighbourhood.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata
ngs.m<-as.data.frame((ngs.r[,-1]>1)*1)
com.ind<-read.table("sample_information.txt", header=T, row.names=1, sep="\t") 

wm1<-rownames(ngs.m) %in% match$code.sequence.weed
wm2<-rownames(ngs.m) %in% match$code.sequence.wheat

ngs.weed.m<-ngs.m[wm1,]
ngs.wheat.m<-ngs.m[wm2,]

ngs.weed.m1<-ngs.weed.m[order(row.names(ngs.weed.m)),]
ngs.wheat.m1<-ngs.wheat.m[order(row.names(ngs.wheat.m)),]

matrix.shared<-matrix(NA,nrow=75,ncol=283)
for (i in 1:75)
{
  matrix.shared[i,]<-as.numeric(colSums(rbind(ngs.weed.m1[i,],ngs.wheat.m1[i,]))==2)
}

matrix.shared<-as.data.frame(matrix.shared)
colnames(matrix.shared)<-colnames(ngs.m)

test<-as.data.frame(colSums(matrix.shared))
test1<-subset(test,test$`colSums(matrix.shared)`>1)

share.num<-rowSums(matrix.shared)

##### Each phylum analysis, richness ####

share.a<-subset(t(fun.shared),taxo.g$Kingdom=="Ascomycota")
share.b<-subset(t(fun.shared),taxo.g$Kingdom=="Basidiomycota")
share.c<-subset(t(fun.shared),taxo.g$Kingdom=="Chytridiomycota")
share.g<-subset(t(fun.shared),taxo.g$Kingdom=="Glomeromycota")
share.z<-subset(t(fun.shared),taxo.g$Kingdom=="Zygomycota")

library(vegan)
###### Ascomycota
shannon.a<- diversity(t(share.a),index="shannon") 
richness.a<-rowSums(t(share.a)>0)
evenness.a<- shannon.a/log(richness.a) 

###### Basidiomycota
shannon.b<- diversity(t(share.b),index="shannon") 
richness.b<-rowSums(t(share.b)>0)
evenness.b<- shannon.b/log(richness.b) 

###### Chytridiomycota
shannon.c<- diversity(t(share.c),index="shannon") 
richness.c<-rowSums(t(share.c)>0)
evenness.c<- shannon.c/log(richness.c) 

###### Glomeromycota
shannon.g<- diversity(t(share.g),index="shannon") 
richness.g<-rowSums(t(share.g)>0)
evenness.g<- shannon.g/log(richness.gh) 

###### Zygomycota
shannon.z<- diversity(t(share.z),index="shannon") 
richness.z<-rowSums(t(share.z)>0)
evenness.z<- shannon.z/log(richness.z) 

richness.share.phy<-cbind(richness.share.fun,richness.a,shannon.a,evenness.a,
                          richness.b,shannon.b,evenness.b,
                          richness.c,shannon.c,evenness.c)

##### Each phylum analysis, abundance ####



##### FUNGuild analysis ####
fun.shared<-matrix.shared[,colnames(matrix.shared)%in% rownames(test1)]
taxo.gf<-taxo.g[rownames(taxo.g)%in% rownames(test1),]

share.path<-subset(t(fun.shared),taxo.gf$trophicMode =="Pathotroph")
share.sym<-subset(t(fun.shared),taxo.gf$trophicMode=="Symbiotroph")
share.sap<-subset(t(fun.shared),taxo.gf$trophicMode=="Saprotroph")

######Pathotroph
shannon.path<- diversity(t(share.path),index="shannon") 
richness.path<-rowSums(t(share.path)>0)
evenness.path<- shannon.path/log(richness.path) 

######Symbiotroph
shannon.sym<- diversity(t(share.sym),index="shannon") 
richness.sym<-rowSums(t(share.sym)>0)
evenness.sym<- shannon.sym/log(richness.sym) 

######Saprotroph
shannon.sap<- diversity(t(share.sap),index="shannon") 
richness.sap<-rowSums(t(share.sap)>0)
evenness.sap<-shannon.sap/log(richness.sap) 

richness.share.fun<-cbind(match,share.num,richness.path,shannon.path,evenness.path,
                    richness.sym,shannon.sym,evenness.sym,
                    richness.sap,shannon.sap,evenness.sap)

#####
mod<-lm(richness.path~code,data=richness.share.fun)
summary(mod)
Anova(mod)

mod<-lm(richness.sym~code,data=richness.share.fun)
summary(mod)
Anova(mod)

mod<-lm(richness.sap~code,data=richness.share.fun)
summary(mod)
Anova(mod)

###### calculate weed-wheat competition ####
ctrl.1<-subset(com.ind,com.ind$weed.neighbour=="wheat")
mean(ctrl.1$Wheat.weight)

ctrl.2<-subset(com.ind,com.ind$weed.neighbour=="control")
mean(ctrl.2$Wheat.weight)

com.ind1<-subset(com.ind,com.ind$weed.neighbour!="control" & com.ind$weed.neighbour!="wheat")

com2<-(mean(ctrl.2$Wheat.weight)-com.ind1$Wheat.weight)/mean(ctrl.2$Wheat.weight)

com.ind1$com1<-(mean(ctrl.1$Wheat.weight)-com.ind1$Wheat.weight)/mean(ctrl.1$Wheat.weight)
com.ind1$com2<-(mean(ctrl.2$Wheat.weight)-com.ind1$Wheat.weight)/mean(ctrl.2$Wheat.weight)







