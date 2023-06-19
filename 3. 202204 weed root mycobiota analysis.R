
################################################################
######NGS data analysis for weed root mycobiota in fields#######

library(ecodist)
library(vegan)
library(ggplot2)
library(gridExtra)

##### loading sequence data and sample information ####
ngs.weed<-read.table("ngs.weed.txt", header=T, row.names=1, sep="\t")  
env.weed<-read.table("env.txt", sep="\t",header=T,row.names=1)        ## sample properties, we need this
taxo<-read.table("taxo.txt", sep="\t",header=T,row.names=1)        ## sample properties, we need this

#####
ngs.pweed<-aggregate(t(ngs.weed), list(taxo$phylum), FUN = sum)    ## otu reads inside phylum,based on normalized matrix
rownames(ngs.pweed)<-ngs.pweed[,1]  
ngs.ptweed<-t(ngs.pweed[,-1])

sum(rowSums(ngs.pweed[,-1]))

#### weed roots 
ngs.pweed$seq<-rowSums(ngs.pweed[,-1])
ngs.pweed$prop<-rowSums(ngs.pweed[,-1])/sum(rowSums(ngs.pweed[,-1]))
ngs.pweed$number<-rowSums(ngs.pweed[,-1])
ngs.pweed$phylum<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Multi-affiliation","GZygomycota")
ngs.pweed$richness<-c(159,66,44,23,3,61)
ngs.pweed[order(ngs.pweed$phylum),]
ngs.pweed$variable<-c("Phyla","Phyla","Phyla","Phyla","Phyla","Phyla")

p.pw= ggplot(ngs.pweed, aes(x=variable,y=prop, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Weed root mycobiota in field")+
  ylab("Percentage of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.pw

p.rw= ggplot(ngs.pweed, aes(x=variable,y=richness, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Weed root mycobiota in field")+
  ylab("Richness of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.rw

pdf("Figure Sx summary.pdf",height=4,width=6)
grid.arrange(p.pw,p.rw,ncol=2)
dev.off() 

#### Calculate diversity index ####
##### Microbial diversity of weed ####
#### overview, calculate value for method description ####
rowSums(ngs.weed)
mean(rowSums(ngs.weed))
min(rowSums(ngs.weed))
max(rowSums(ngs.weed))

library(vegan)
rarecurve(ngs.weed,step=100,xlim=c(0,8200),ylim=c(0,180),xlab="Sequence depth", ylab="Number of sequence clusters in weed root mycobiota",col="blue",label=FALSE)

otu.shannon.weed<- diversity(ngs.weed, index = "shannon") 
otu.richness.weed<-rowSums(ngs.weed>0)
otu.evenness.weed<- otu.shannon.weed/log(otu.richness.weed)        ## Pielou's evenness
otu.abundance.weed<-rowSums(ngs.weed)

richness.weed<-cbind(env.weed,otu.richness.weed,otu.shannon.weed,otu.evenness.weed,otu.abundance.weed)

#### For weed root mycobiota
lm.r<-glmer.nb(otu.richness.root~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

lm.r<-lmer(otu.evenness.root~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

###### PCoA analysis####
ngs.log<-log2(ngs.weed+1)

ngs.bray.all<-vegdist(ngs.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.all<-pco(ngs.bray.all,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.all$vectors[,2]~pcoaVS.all$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=as.factor(env.weed$weed),
     axes = TRUE, main = "PCoA (ecodist) on Soil and Root")

pcoa.all<-cbind(pcoaVS.all$vectors[,2],pcoaVS.all$vectors[,1],env.weed)


pdf("PCoA of weed root mycobiota in field.pdf",height=5,width=7)
ggplot(data=pcoa.all,aes(x=pcoaVS.all$vectors[, 1],y=pcoaVS.all$vectors[, 2]))+
  xlim(-0.25,0.25)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)),size=3)+
  scale_color_manual(values=c("Lamium purpureum"="#800080",
                               "Matricaria sp"="#008000","Poa annua"="#FF0000",
                               "Papaver rhoeas"="#44A8DB","Poa trivialis"="#FB8072",
                               "Trifolium repens"="#F962C3","Veronica persica"="#A52A2A",
                               "Vicia sativa"="#0000FF"))+
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=weed))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on weed root mycobiota")
dev.off()


#### statistics to check the effect of weed species on root mycobiota community composition ####
adonis.all<-adonis(ngs.log~weed,env.weed,by=NULL, method="bray") 
adonis.all
rda(ngs.log)

###### For each phylum ####
ngs.asco<-subset(t(ngs.weed),taxo$phylum=="Ascomycota")
ngs.basi<-subset(t(ngs.weed),taxo$phylum=="Basidiomycota")
ngs.chyt<-subset(t(ngs.weed),taxo$phylum=="Chytridiomycota")
ngs.glom<-subset(t(ngs.weed),taxo$phylum=="Glomeromycota")
ngs.zygo<-subset(t(ngs.weed),taxo$phylum=="Zygomycota")

#####
##### calculate diversity index inside each phylum for weed root ####
library(vegan)
#otu.shannon.asco.r<- diversity(t(ngs.asco),index="shannon") 
otu.richness.asco.r<-rowSums(t(ngs.asco)>0)
#otu.evenness.asco.r<- otu.shannon.asco.r/log(otu.richness.asco.r)  ## Pielou's evenness
otu.abundance.asco.r<-rowSums(t(ngs.asco))

########Basidiomycota
#otu.shannon.basi.r<- diversity(t(ngs.basi),index="shannon") 
otu.richness.basi.r<-rowSums(t(ngs.basi)>0)
#otu.evenness.basi.r<- otu.shannon.basi.r/log(otu.richness.basi.r)  ## Pielou's evenness
otu.abundance.basi.r<-rowSums(t(ngs.basi))

##### Chytridiomycota
#otu.shannon.chyt.r<- diversity(t(ngs.chyt),index="shannon") 
otu.richness.chyt.r<-rowSums(t(ngs.chyt)>0)
#otu.evenness.chyt.r<- otu.shannon.chyt.r/log(otu.richness.chyt.r)  ## Pielou's evenness
otu.abundance.chyt.r<-rowSums(t(ngs.chyt))

######Glomeromycota
#otu.shannon.glom.r<- diversity(t(ngs.glom),index="shannon") 
otu.richness.glom.r<-rowSums(t(ngs.glom)>0)
#otu.evenness.glom.r<- otu.shannon.glom.r/log(otu.richness.glom.r)  ## Pielou's evenness
otu.abundance.glom.r<-rowSums(t(ngs.glom))

######Zygomycota
#otu.shannon.zygo.r<- diversity(t(ngs.zygo),index="shannon") 
otu.richness.zygo.r<-rowSums(t(ngs.zygo)>0)
#otu.evenness.zygo.r<- otu.shannon.zygo.r/log(otu.richness.zygo.r)  ## Pielou's evenness
otu.abundance.zygo.r<-rowSums(t(ngs.zygo))

richness.phylum<-cbind(richness.weed,otu.richness.asco.r,otu.abundance.asco.r,
                       otu.richness.basi.r,otu.abundance.basi.r,
                       otu.richness.chyt.r,otu.abundance.chyt.r,
                       otu.richness.glom.r,otu.abundance.glom.r,
                       otu.richness.zygo.r,otu.abundance.zygo.r)

### Statistics for diversity inside each phylum ####
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

lm.r<-glmer.nb(otu.richness.weed~weed+(1|field),richness.phylum)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

lm.r<-glmer.nb(otu.abundance.weed~weed+(1|field),richness.phylum)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)


### effect for Ascomycota####
#lm.asc<-lmer(otu.richness.asco~log2(weed)*log2(otu.richness.soil)+(1|field)+(1|quadrats),richness.phylum)
lm.asc<-glmer.nb(otu.richness.asco.r~weed+(1|field),richness.phylum)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.abundance.asco.r~weed+(1|field),richness.phylum)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

#### effect for Basidiomycota ####
lm.bas<-glmer.nb(otu.richness.basi.r~weed+(1|field),richness.phylum)
#lm.bas<-glmer.nb(otu.richness.basi~weed+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.abundance.basi.r~weed+(1|field),richness.phylum)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

#### effect for Chytridiomycota####
lm.chy<-glmer.nb(otu.richness.chyt.r~weed+(1|field),richness.phylum)
#lm.chy<-glmer.nb(otu.richness.chyt~weed+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.abundance.chyt.r~weed+(1|field),richness.phylum)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

#### effect for Glomeromycota ####
lm.glom<-glmer.nb(otu.richness.glom.r~weed+(1|field),richness.phylum)
#lm.glom<-glmer.nb(otu.richness.glom~weed+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.glom)
AIC(lm.glom)
r.squaredGLMM(lm.glom)

lm.glom<-lmer(otu.abundance.glom.r~weed+(1|field),richness.phylum)
Anova(lm.glom)
AIC(lm.glom)
r.squaredGLMM(lm.glom)

####effect on Zygomycota ####
lm.zygo<-glmer.nb(otu.richness.zygo.r~weed+(1|field),richness.phylum)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.abundance.zygo.r~weed+(1|field),richness.phylum)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

###### FUNGuild dataset#####
taxo.fun<-read.table("taxo.guilds_field_weeds.txt", sep="\t",header=T,row.names=1)        ## sample properties, we need this
fun<-rownames(taxo) %in% rownames(taxo.fun)
ngs.fun<-ngs.weed[,fun]

ngs.path<-subset(t(ngs.fun),taxo.fun$trophicMode =="Pathotroph")
ngs.sym<-subset(t(ngs.fun),taxo.fun$trophicMode=="Symbiotroph")
ngs.sap<-subset(t(ngs.fun),taxo.fun$trophicMode=="Saprotroph")

######Pathotroph
#otu.shannon.path<- diversity(t(ngs.path),index="shannon") 
otu.richness.path<-rowSums(t(ngs.path)>0)
#otu.evenness.path<- otu.shannon.path/log(otu.richness.path) 
otu.abundance.path<-rowSums(t(ngs.path))

######Symbiotroph
#otu.shannon.sym<- diversity(t(ngs.sym),index="shannon") 
otu.richness.sym<-rowSums(t(ngs.sym)>0)
#otu.evenness.sym<- otu.shannon.sym/log(otu.richness.sym) 
otu.abundance.sym<-rowSums(t(ngs.sym))

######Saprotroph
#otu.shannon.sap<- diversity(t(ngs.sap),index="shannon") 
otu.richness.sap<-rowSums(t(ngs.sap)>0)
#otu.evenness.sap<- otu.shannon.sap/log(otu.richness.sap) 
otu.abundance.sap<-rowSums(t(ngs.sap))

richness.fun<-cbind(richness.phylum,otu.richness.path,otu.abundance.path,
                       otu.richness.sym,otu.abundance.sym,
                       otu.richness.sap,otu.abundance.sap)
####effect on Pathotroph ####
lm.path<-glmer.nb(otu.richness.path~weed+(1|field),richness.fun)
Anova(lm.path)
AIC(lm.path)
r.squaredGLMM(lm.path)

lm.path<-lmer(otu.abundance.path~weed+(1|field),richness.fun)
Anova(lm.path)
AIC(lm.path)
r.squaredGLMM(lm.path)

####effect on Symbiotroph ####
lm.sym<-glmer.nb(otu.richness.sym~weed+(1|field),richness.fun)
Anova(lm.sym)
AIC(lm.sym)
r.squaredGLMM(lm.sym)

lm.sym<-lmer(otu.abundance.sym~weed+(1|field),richness.fun)
Anova(lm.sym)
AIC(lm.sym)
r.squaredGLMM(lm.sym)

####effect on Saprotroph ####
lm.sap<-glmer.nb(otu.richness.sap~weed+(1|field),richness.fun)
Anova(lm.sap)
AIC(lm.sap)
r.squaredGLMM(lm.sap)

lm.sap<-lmer(otu.abundance.sap~weed+(1|field),richness.fun)
Anova(lm.sap)
AIC(lm.sap)
r.squaredGLMM(lm.sap)


