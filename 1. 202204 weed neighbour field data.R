##############################################
###### Data format, 202203, Jie Hu ####
#I would rather put a figure with floristic coverage 
#in total field, 
#in organic fields, 
#in conventional fields, 
#in total border, 
#in hedgerow borders, 
#in grassy borders

floristic.edge<-read.table("floristic_edge.txt", header=T, row.names=1, sep="\t")  
floristic.crop<-read.table("floristic_crop.txt", header=T, row.names=1, sep="\t")

##### calaulation of floristic diversity within crops ####
#### overview, calculate value for method description ####
floristic.crop.c<-subset(floristic.crop,floristic.crop$type=="conventional")
floristic.crop.o<-subset(floristic.crop,floristic.crop$type=="organic ")

test.crop.c<-as.matrix(colSums(floristic.crop.c[,-1:-4]))
test.crop.o<-as.matrix(colSums(floristic.crop.o[,-1:-4]))
test.bsh<-as.matrix(colSums(floristic.Bsht))
test.bh<-as.matrix(colSums(floristic.Bht))

test.crop1.c<-subset(test.crop.c,test.crop.c>0)
test.crop1.o<-subset(test.crop.o,test.crop.o>0)
test.bsh1<-subset(test.bsh,test.bsh>0)
test.bh1<-subset(test.bh,test.bh>0)

pec.crop<-as.matrix(colSums(floristic.crop[,-1:-4]))/sum(colSums(floristic.crop[,-1:-4]))
pec.bsh<-as.matrix(colSums(floristic.Bsht))/sum(colSums(floristic.Bsht))
pec.bh<-as.matrix(colSums(floristic.Bht))/sum(colSums(floristic.Bht))

pec.crop1<-subset(pec.crop,pec.crop>0.03)
pec.bsh1<-subset(pec.bsh,pec.bsh>0.03)
pec.bh1<-subset(pec.bh,pec.bh>0.03)


mean(rowSums(floristic.crop[,-1:-4]))
min(rowSums(floristic.crop[,-1:-4]))
max(rowSums(floristic.crop[,-1:-4]))
floristic.t<-as.data.frame(floristic.crop[,-1:-4])

library(vegan)
otu.shannon.f<- diversity(floristic.t, index = "shannon") 
otu.richness.f<-rowSums(floristic.t>0)
otu.evenness.f<- otu.shannon.f/log(otu.richness.f)        ## Pielou's evenness

richness.f.c<-cbind(floristic.crop[,1:4],otu.richness.f,otu.shannon.f,otu.evenness.f)

##### focus on 8 weed species with crops ####
floristic.crop.a<-aggregate(floristic.crop, list(floristic.crop$field,floristic.crop$type), FUN=mean)    ## pool data
test<-stack(floristic.crop.a[,-1:-2])

#floristic.crop.sa <- floristic.crop[,c("Group.1","Group.2","Lamium_purpureum","Matricaria_sp.","Poa_annua","Papaver_rhoeas",
                                        # "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_sativa")]

floristic.crop.sa <- floristic.crop.a[,c("Group.1","Group.2","Lamium_purpureum","Matricaria_sp.","Poa_annua","Papaver_rhoeas",
                                "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_sativa")]

colnames(floristic.crop.sa)<-c("field","type","Lamium_purpureum","Matricaria_sp.","Poa_annua","Papaver_rhoeas",
                               "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_sativa")

test<-stack(floristic.crop.s[,-1:-2])
test1<-rep(floristic.crop.s[,1],8)
test2<-rep(floristic.crop.s[,2],8)

floristic.crop.s1<-cbind(test1,test2,test)

floristic.crop.01<-cbind(floristic.crop[,1:4],1*(floristic.crop[,-1:-4]>0))
floristic.crop.s01 <- floristic.crop.01[,c("type","field","Lamium_purpureum","Matricaria_sp.","Poa_annua","Papaver_rhoeas",
                                        "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_sativa")]

floristic.crop.a01<-aggregate(floristic.crop.s01[,-1:-2], list(floristic.crop$field,floristic.crop$type), FUN=sum)    ## pool data
test01<-stack(floristic.crop.a01[,-1:-2])

floristic.crop.sa01<-cbind(test1,test2,test01)

library(ggplot2)
ggplot(floristic.crop.sa01,aes(x=test2,y=values))+geom_boxplot(aes(fill=as.factor(test2)))

##### Floristic survey at edges ######
##### shape the data format first ####
floristic.shape<-read.table("DATA_flore_AGRIM.txt", header=T, row.names=1, sep="\t")  
floristic.shape.a<-aggregate(floristic.shape, list(floristic.shape$ID,floristic.shape$Localisation,floristic.shape$Date),FUN=mean)    ## pool data
floristic.shape.s<-aggregate(floristic.shape, list(floristic.shape$Espece),FUN=mean)    ## pool data

Bh<-subset(floristic.shape.a,floristic.shape.a$Group.2=="Bh")
Bsh<-subset(floristic.shape.a,floristic.shape.a$Group.2=="Bsh")

floristic.Bh<-subset(floristic.shape,floristic.shape$Localisation=="Bh")
floristic.Bsh<-subset(floristic.shape,floristic.shape$Localisation=="Bsh")

floristic.shape.sBh<-aggregate(floristic.Bh, list(floristic.Bh$Espece),FUN=mean)    ## pool data
floristic.shape.sBsh<-aggregate(floristic.Bsh, list(floristic.Bsh$Espece),FUN=mean)    ## pool data

floristic.shape.c<-matrix(NA,nrow=420,ncol=193)

Qua<-as.data.frame(rep(c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10"),42))
qua1<-as.data.frame(Qua[order(Qua),])
colnames(qua1)<-"qua"

floristic.shape.Bh<-cbind(rep(floristic.shape.a$Group.1,5),qua1,floristic.shape.c)
colnames(floristic.shape.Bh)<-c("code","qua",floristic.shape.sBh$Group.1)

for (i in 1:193)
{ espece<-subset(floristic.Bh,floristic.Bh$Espece==colnames(floristic.shape.Bh)[i+2])
  for (k in 1:dim(floristic.shape.Bh)[1])
  {
    for (j in 1:dim(espece)[1])
       if ( espece$ID[j]==floristic.shape.Bh[,1][k] & espece$Quadrat[j]==floristic.shape.Bh[,2][k])

            {print ("yes") 
              floristic.shape.Bh[k,i+2]<-espece[j,10]}
  }
}

floristic.shape.cs<-matrix(NA,nrow=420,ncol=176)

Qua<-as.data.frame(rep(c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10"),42))
qua1<-as.data.frame(Qua[order(Qua),])
colnames(qua1)<-"qua"

#colnames(floristic.shape.c)<-floristic.shape.s$Group.1
floristic.shape.Bsh<-cbind(rep(floristic.shape.a$Group.1,5),qua1,floristic.shape.cs)
colnames(floristic.shape.Bsh)<-c("code","qua",floristic.shape.sBsh$Group.1)

for (i in 1:176)
{ espece<-subset(floristic.Bsh,floristic.Bsh$Espece==colnames(floristic.shape.Bsh)[i+2])
for (k in 1:dim(floristic.shape.Bsh)[1])
{
  for (j in 1:dim(espece)[1])
    if ( espece$ID[j]==floristic.shape.Bsh[,1][k] & espece$Quadrat[j]==floristic.shape.Bsh[,2][k])
      
    {print ("yes") 
      floristic.shape.Bsh[k,i+2]<-espece[j,10]}
  #else {floristic.shape.Bh[k,i+2]<-0}
}
}

floristic.shape.Bsh[is.na(floristic.shape.Bsh)] <- 0 
floristic.shape.Bh[is.na(floristic.shape.Bh)] <- 0 

#floristic.shape.con<-rbind(as.data.frame(floristic.shape.Bh),as.data.frame(floristic.shape.Bsh))
#write.table(floristic.shape.con,file="floristic.shape.con.txt",sep="\t")
#### calculate floristic diversity at edges ####
#floristic.edge.new<-read.table("floristic.shape.con1.txt", header=T, row.names=1, sep="\t")

mean(rowSums(floristic.shape.Bh[,-1:-2]))
min(rowSums(floristic.shape.Bh[,-1:-2]))
max(rowSums(floristic.shape.Bh[,-1:-2]))
floristic.Bht<-as.data.frame(floristic.shape.Bh[,-1:-2])
floristic.Bsht<-as.data.frame(floristic.shape.Bsh[,-1:-2])

#library(vegan)
#rarecurve(floristic.Bht,step=1,xlab="Sequence depth", ylab="Number of weed species in field",col="blue",label=FALSE)

otu.shannon.Bh<- diversity(floristic.Bht, index = "shannon") 
otu.richness.Bh<-rowSums(floristic.Bht>0)
otu.evenness.Bh<- otu.shannon.Bh/log(otu.richness.Bh)        ## Pielou's evenness

otu.shannon.Bsh<- diversity(floristic.Bsht, index = "shannon") 
otu.richness.Bsh<-rowSums(floristic.Bsht>0)
otu.evenness.Bsh<- otu.shannon.Bsh/log(otu.richness.Bsh)        ## Pielou's evenness

richness.f.e<-cbind(floristic.shape.Bh[,1:2],otu.richness.Bh,otu.shannon.Bh,otu.evenness.Bh,otu.richness.Bsh,otu.shannon.Bsh,otu.evenness.Bsh)

##### Focus on 8 weed species at edges ####
floristic.Bh.lp<-floristic.shape.Bh$"Lamium purpureum"
floristic.Bh.ma<-floristic.shape.Bh$"Matricaria chamomilla"
floristic.Bh.pr<-floristic.shape.Bh$"Papaver rhoeas"
floristic.Bh.pa<-floristic.shape.Bh$"Poa annua"
floristic.Bh.pt<-floristic.shape.Bh$"Poa trivialis"
floristic.Bh.tr<-floristic.shape.Bh$"Trifolium repens"
floristic.Bh.vp<-floristic.shape.Bh$"Veronica persica"
floristic.Bh.vs<-floristic.shape.Bh$"Vicia sativa"
floristic.Bh.s<-cbind(floristic.shape.Bh[,1:2],floristic.Bh.lp,floristic.Bh.ma,floristic.Bh.pa,floristic.Bh.pr,floristic.Bh.pt,floristic.Bh.tr,floristic.Bh.vp,floristic.Bh.vs)

floristic.Bsh.lp<-floristic.shape.Bsh$"Lamium purpureum"
floristic.Bsh.ma<-floristic.shape.Bsh$"Matricaria chamomilla"
floristic.Bsh.pr<-floristic.shape.Bsh$"Papaver rhoeas"
floristic.Bsh.pa<-floristic.shape.Bsh$"Poa annua"
floristic.Bsh.pt<-floristic.shape.Bsh$"Poa trivialis"
floristic.Bsh.tr<-floristic.shape.Bsh$"Trifolium repens"
floristic.Bsh.vp<-floristic.shape.Bsh$"Veronica persica"
floristic.Bsh.vs<-floristic.shape.Bsh$"Vicia sativa"
floristic.Bsh.s<-cbind(floristic.shape.Bsh[,1:2],floristic.Bsh.lp,floristic.Bsh.ma,floristic.Bsh.pa,floristic.Bsh.pr,floristic.Bsh.pt,floristic.Bsh.tr,floristic.Bsh.vp,floristic.Bsh.vs)

write.table(richness.f.e,file="richness.f.e.txt",sep="\t")
write.table(floristic.Bh.s,file="floristic.Bh.s.txt",sep="\t")
write.table(floristic.Bsh.s,file="floristic.Bsh.s.txt",sep="\t")
write.table(richness.f,file="richness.f.txt",sep="\t")

richness.f.edge<-read.table("richness.f.e1.txt", header=T, row.names=1, sep="\t")  
floristic.edge.s<-read.table("floristic.edge.s.txt", header=T, row.names=1, sep="\t") 
richness.f1<-read.table("richness.f1.txt", header=T, row.names=1, sep="\t") 


##### pool data from all quadrats together at field scale ####
floristic.Bht.a<-aggregate(floristic.shape.Bh, list(floristic.shape.Bh$code),FUN=mean)    ## pool data
floristic.Bsht.a<-aggregate(floristic.shape.Bsh, list(floristic.shape.Bsh$code),FUN=mean)    ## pool data
floristic.c.a<-aggregate(floristic.crop, list(floristic.crop$ID,floristic.crop$type),FUN=mean)    ## pool data

otu.shannon.Bhf<- diversity(floristic.Bht.a[,-1:-3], index = "shannon") 
otu.richness.Bhf<-rowSums(floristic.Bht.a[,-1:-3]>0)
otu.evenness.Bhf<- otu.shannon.Bhf/log(otu.richness.Bhf)      

otu.shannon.Bshf<- diversity(floristic.Bsht.a[,-1:-3], index = "shannon") 
otu.richness.Bshf<-rowSums(floristic.Bsht.a[,-1:-3]>0)
otu.evenness.Bshf<- otu.shannon.Bshf/log(otu.richness.Bshf)      

otu.shannon.cf<- diversity(floristic.c.a[,-1:-6], index = "shannon") 
otu.richness.cf<-rowSums(floristic.c.a[,-1:-6]>0)
otu.evenness.cf<- otu.shannon.cf/log(otu.richness.cf)

richness.ff.Bh<-cbind(floristic.Bht.a[,1:3],otu.richness.Bhf,otu.shannon.Bhf,otu.evenness.Bhf)
richness.ff.Bsh<-cbind(floristic.Bsht.a[,1:3],otu.richness.Bshf,otu.shannon.Bshf,otu.evenness.Bshf)
richness.ff.c<-cbind(floristic.c.a[,1:3],otu.richness.cf,otu.shannon.cf,otu.evenness.cf)

write.table(richness.ff.Bh,file="richness.ff.Bh.txt",sep="\t")
write.table(richness.ff.Bsh,file="richness.ff.Bsh.txt",sep="\t")
write.table(richness.ff.c,file="richness.ff.c.txt",sep="\t")

richness.ff.all<-read.table("richness.ff.all.txt", header=T, row.names=1, sep="\t")

##### Selected weed species ####
floristic.edge.sa<-aggregate(floristic.edge.s[,-1:-4], list(floristic.edge.s$code,floristic.edge.s$type,floristic.edge.s$local), FUN=mean)    ## pool data

##### Some mixed linear models ####
#### loading the packages #### 
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

#mod<-lmer(otu.evenness.root~otu.richness.f*otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#mod<-glmer.nb(otu.richness.soil~otu.richness.f+(1|field)+(1|quadrats), data=richness.all1)
mod<-glmer.nb(otu.richness.root~otu.richness.f+otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.shannon.f,by="otu.shannon.soil",data=pcoa.root.r)
lsmeans(mod,pairwise~otu.shannon.soil,by="otu.shannon.f",data=richness.all1)

mod<-lmer(otu.evenness.root~otu.richness.f+otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.shannon.f,by="otu.shannon.soil",data=pcoa.root.r)
lsmeans(mod,pairwise~otu.shannon.soil,by="otu.shannon.f",data=richness.all1)

############################################
####' Check the correlation #####
library(corrplot)
coverage.all<-cbind(floristic.crop.s1,floristic.edge.bh2,floristic.edge.bsh2)
write.table(coverage.all,file="coverage.all.txt",sep="\t")
coverage.cor<-read.table("coverage.all1.txt", header=T, row.names=1, sep="\t")
cor_data<-coverage.cor[,-1:-2]

res2 <- cor(as.matrix(cor_data))
res2
#res2$r
#res2$P

corrplot(res2,type="upper",order="hclust",method='number')
library(arm)




