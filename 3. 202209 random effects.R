###Statistics for making Tables, Jie Hu and Cendrine Mony
#### loading the packages 
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

###### At individual level####
richness.tp<-subset(richness,richness$plant=="Trifolium"&richness$position!="cor")
### richness
### global fungi
###for richness try both
#?lmer and glmer.nb ####### wait a bit 
#mod<-glmer.nb(otu.richness~code*corridor+(1|pool.code), data=richness.tp) 
#mod<-lmer(otu.richness~code*corridor+(1|pool.code), data=richness.tp)
mod<-glmer.nb(otu.richness~code*corridor+(1|pool.code), data=richness.tp)
#mod<-glm.nb(otu.richness~code*corridor, data=richness.tp)

summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.tp)
lsmeans(mod,pairwise~code,by="corridor",data=richness.tp)

#### Ascomycota
#write.table(richness.with.phylum.nomb,file="Each phylum richness at individual level.txt",sep="\t")
#mod<-glmer.nb(otu.richness.asco~code*corridor+(1|pool.code),richness.with.phylum.nomb)
mod<-lmer(otu.richness.asco~code*corridor+(1|pool.code),richness.with.phylum.nomb)
mod<-glmer.nb(otu.richness.asco~code*corridor+(1|pool.code),richness.with.phylum.nomb)
mod<-glm.nb(otu.richness.asco~code*corridor,richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

#### Corridor effect for Basidiomycota
#initCtrl=list(limit=1000)
#mod<-glmer.nb(otu.richness.Basidio~code*corridor+(1|pool.code),data=richness.with.phylum.nomb)
mod<-glmer.nb(otu.richness.Basidio~code*corridor+(1|pool.code),data=richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

#### Corridor effect for Chytridiomycota
mod<-lmer(otu.richness.Chytridio~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

#### Corridor effect for Glomeromycota
mod<-lmer(otu.richness.Glomero~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

#### Corridor effect on Zygomycota
mod<-lmer(otu.richness.Zygo~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

####Evenness####
### global fungi
#mod<-glmer.nb(otu.evenness~code*corridor+(1|pool.code), data=richness.tp)
mod<-lmer(otu.evenness~code*corridor+(1|pool.code), data=richness.tp)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.tp)
lsmeans(mod,pairwise~code,by="corridor",data=richness.tp)

#### Ascomycota
mod<-lmer(otu.evenness.asco~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

#### Basidiomycota
mod<-lmer(otu.evenness.Basidio~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

####Chytridiomycota
mod<-lmer(otu.evenness.Chytridio~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

####Glomeromycota
mod<-lmer(otu.evenness.Glomero~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

####Zygomycota
mod<-lmer(otu.evenness.Zygo~code*corridor+(1|pool.code),richness.with.phylum.nomb)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="code",data=richness.with.phylum.nomb)
lsmeans(mod,pairwise~code,by="corridor",data=richness.with.phylum.nomb)

####Patch level####
########## Rand effect model for global fungi at patch level 
####Richness
mod<-lmer(otu.richness~Group.4*Group.6+(1|Group.7),data=richness.tpp)
mod<-glmer.nb(otu.richness~Group.4*Group.6+(1|Group.7),data=richness.tpp)
mod<-glm.nb(otu.richness~Group.4*Group.6,data=richness.tpp)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.tpp)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.tpp)

write.table(richness.tpp,file="Global richness at patch level.txt",sep="\t")
#### Ascomycota
write.table(richness.pool.phylum.tri,file="Each phylum richness at patch level.txt",sep="\t")
mod<-lmer(otu.richness.asco.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
mod<-glmer.nb(otu.richness.asco.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
mod<-glm.nb(otu.richness.asco.p~Group.4*Group.6,data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

#### Basidiomycota
mod<-lmer(otu.richness.Basidio.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Chytridiomycota
mod<-lmer(otu.richness.Chytridio.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Glomeromycota
mod<-lmer(otu.richness.Glomero.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Zygomycota
mod<-lmer(otu.richness.Zygo.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Evenness
mod<-lmer(otu.evenness~Group.4*Group.6+(1|Group.7),data=richness.tpp)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.tpp)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.tpp)

#### Ascomycota
mod<-lmer(otu.evenness.asco.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

#### Basidiomycota
mod<-lmer(otu.evenness.Basidio.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Chytridiomycota
mod<-lmer(otu.evenness.Chytridio.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Glomeromycota
mod<-lmer(otu.evenness.Glomero.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

####Zygomycota
mod<-lmer(otu.evenness.Zygo.p~Group.4*Group.6+(1|Group.7),data=richness.pool.phylum.tri)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.6,by="Group.4",data=richness.pool.phylum.tri)
lsmeans(mod,pairwise~Group.4,by="Group.6",data=richness.pool.phylum.tri)

#####landscape level####
########## Random effect model for global fungi at landscape level 
####Richness
####Global fungi
library(MASS)
mod<-glm.nb(otu.richness.pp~campaign*corridor,data=richness.tpmp)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

write.table(richness.tpmp,file="Global richness at landscape level.txt",sep="\t")
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="campaign",data=richness.tpmp)
lsmeans(mod,pairwise~campaign,by="corridor",data=richness.tpmp)

#### Ascomycota
write.table(richness.pp.phylum,file="Each phylum richness at landscape level.txt",sep="\t")
mod<-glm.nb(otu.richness.asco.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

#### Basidiomycota
mod<-glm(otu.richness.Basidio.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Chytridiomycota
mod<-glm(otu.richness.Chytridio.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Glomeromycota
mod<-glm(otu.richness.Glomero.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Zygomycota
mod<-glm(otu.richness.Zygo.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Evenness
#### Global fungi
mod<-glm(otu.evenness.pp~campaign*corridor,data=richness.tpmp)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~corridor,by="campaign",data=richness.tpmp)
lsmeans(mod,pairwise~campaign,by="corridor",data=richness.tpmp)

#### Ascomycota
mod<-glm(otu.evenness.asco.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

#### Basidiomycota
mod<-glm(otu.evenness.Basidio.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Chytridiomycota
mod<-glm(otu.evenness.Chytridio.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Glomeromycota
mod<-glm(otu.evenness.Glomero.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

####Zygomycota
mod<-glm(otu.evenness.Zygo.pp~Group.2*Group.3,data=richness.pp.phylum)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))
####group comparison after model corrections
lsmeans(mod,pairwise~Group.3,by="Group.2",data=richness.pp.phylum)
lsmeans(mod,pairwise~Group.2,by="Group.3",data=richness.pp.phylum)

#### Random effect for global fungal Bray-Curtis dissimilarity at individual scale
bray.i<-read.table("Bray curtis distance at individual level.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata
mod<-lmer(distance~code*corridor+(1|mesocosm),data=bray.i)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~code,by="corridor",data=bray.i)
lsmeans(mod,pairwise~corridor,by="code",data=bray.i)

#####
dis.phylum.all1<-read.table("Bray curtis distance at individual level All phylum1.txt", header=T, sep="\t")  ## otu abundance rawdata

mod<-lmer(Ascomycota~time*corridor+(1|mesocosm),data=dis.phylum.all1)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~corridor,by="time",data=dis.phylum.all1)


mod<-lmer(Basidiomycota~time*corridor+(1|mesocosm),data=dis.phylum.all1)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~corridor,by="time",data=dis.phylum.all1)

mod<-lmer(Chytridiomycota~time*corridor+(1|mesocosm),data=dis.phylum.all1)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~corridor,by="time",data=dis.phylum.all1)


mod<-lmer(Glomeromycota~time*corridor+(1|mesocosm),data=dis.phylum.all1)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~corridor,by="time",data=dis.phylum.all1)


mod<-lmer(Zygomycota~time*corridor+(1|mesocosm),data=dis.phylum.all1)
summary(mod)
Anova(mod)
rsq(mod)
hist(residuals(mod))

lsmeans(mod,pairwise~corridor,by="time",data=dis.phylum.all1)
