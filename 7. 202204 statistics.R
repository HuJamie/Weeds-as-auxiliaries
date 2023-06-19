#################################################
###Statistical test ####
library(ggplot2)
library(ggprism)
#### statistics for weed identity effect on aboveground dry weight, microbiome diversity and composition ####
# 1 -----------------------------------------------------------------------
ONE_Tukey_HSD1 <- function(data,group,compare,value){
  
  library(multcomp)#Tukey检验需要用到这个包来标显著性字母标记
  
  a <- data.frame(stringsAsFactors = F)#做一个空的数据框
  type <- unique(data[,group])#统计需要运行多重比较的次数
  for (i in type)#进行type次多重比较
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]#根据指定的i去取相应的数据集出来
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    names(sub_dat)[names(sub_dat)==compare] <- 'g1' ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==value] <- 'value' ## 重命名方便后面使用
    sub_dat$g1 <- factor(sub_dat$g1)#将列转化成因子以进行多重比较
    
    fit <- aov(value ~ g1,data = sub_dat )#方差分析
    #Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)#Tukey检验多重比较
    Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)#获取多重比较字母标注
    Tukey.labels$compare = rownames(Tukey.labels)## 提取字母分组行名为group组名
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),#获取数据标准差
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"#获取数据均值
    )
    names(mean_sd) <- c('compare','std','mean')#列名重命名
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))#合并数据
  }
  
  names(a) <- c(compare,'std','mean','Letters',group)#列名重命名
  return(a)
}

#### Weed richness in each crop field ####
library(agricolae)
pcoa.weed$weed.neighbour<-as.factor(pcoa.weed$weed.neighbour)
stat.pc1.weed<-aov(`pcoaVS.weed$vectors[, 1]`~weed.neighbour,pcoa.weed)
comparison.we1<-LSD.test(stat.pc1.weed,"weed.neighbour")

## floristic coverage of each weed species
floristic.crop.s1.o<-subset(floristic.crop.s1,floristic.crop.s1$test2=="organic "&floristic.crop.s1$values<15)
sta.o.fc <- ONE_Tukey_HSD1(floristic.crop.s1.o,'test2','ind','values')

pdf("Figure S2. floristic coverage of each weed species.pdf",height=5,width=7)
ggplot(floristic.crop.s1.o,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(ind)))+
  geom_point(position=position_dodge(width=0.75),aes(group=ind,color=as.factor(ind)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.o.fc,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~local,scales = "free_y")+ 
  labs(x='',y='Floristic coverage of each weed species in organic fields')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

## At quadrat level
sta.rich.fc <- ONE_Tukey_HSD1(richness.f1,'code','type','otu.richness.f')
sta.rich.fet <- ONE_Tukey_HSD1(richness.f.edge,'type','local','otu.richness')
sta.rich.fel <- ONE_Tukey_HSD1(richness.f.edge,'local','type','otu.richness')
sta.even.fet <- ONE_Tukey_HSD1(richness.f.edge,'type','local','otu.evenness')
sta.even.fel <- ONE_Tukey_HSD1(richness.f.edge,'local','type','otu.evenness')

rich.con<-subset(richness.f.edge,richness.f.edge$type=="conventional")
rich.orga<-subset(richness.f.edge,richness.f.edge$type=="organic")
  
pdf("Figure 1 weed richness in different locations in fields.pdf",height=4,width=6)
ggplot(richness.f.edge,aes(x=type,y=otu.richness))+geom_boxplot(aes(fill=as.factor(local)))+
  geom_point(position=position_dodge(width=0.75),aes(group=local,color=as.factor(local)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.rich.fet,aes(x=type,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~local,scales = "free_y")+ 
  labs(x='',y='Weed richness at different locations in fields')+
  ggprism::theme_prism()#+theme(axis.text.x = element_text(angle=90))
dev.off()

## At field scale
sta.rich.ff <- ONE_Tukey_HSD1(richness.ff.all,'local','type','otu.richness')
sta.rich.ff1 <- ONE_Tukey_HSD1(richness.ff.all,'type','local','otu.richness')

pdf("Figure 1 weed richness at field scale.pdf",height=4,width=6)
ggplot(richness.ff.all,aes(x=type,y=otu.richness))+geom_boxplot(aes(fill=as.factor(local)))+
  geom_point(position=position_dodge(width=0.75),aes(group=local,color=as.factor(local)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.rich.ff,aes(x=type,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~local,scales = "free_y")+ 
  labs(x='',y='Weed richness at field scale')+
  ggprism::theme_prism()#+theme(axis.text.x = element_text(angle=90))
dev.off()

#### selected eight weed species ####

####
sta.crop.w <- ONE_Tukey_HSD1(floristic.crop.s1,'test2','ind','values')
sta.crop.t <- ONE_Tukey_HSD1(floristic.crop.s1,'ind','test2','values')

pdf("Figure 1.1. floristic coverage within crops.pdf",height=5,width=10)
ggplot(floristic.crop.s1,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(ind)))+
  geom_point(position=position_dodge(width=0.75),aes(group=ind,color=as.factor(ind)))+ 
  ylim(0,22)+
  geom_text(data=sta.crop.w,aes(x=ind,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~test2,scales = "free_y")+ labs(x='',y='Weed coverage in crop field')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("Figure 1.2.pdf",height=5,width=8)
ggplot(floristic.crop.s1,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(test2)))+
  geom_point(position=position_dodge(width=0.75),aes(group=test2,color=as.factor(test2)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.crop.t,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~test2,scales = "free_y")+ 
  labs(x='Different Weeds',y='Weed abundance in crop field')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()


#######
sta.crop.w <- ONE_Tukey_HSD1(floristic.crop.sa01,'test2','ind','values')
sta.crop.t <- ONE_Tukey_HSD1(floristic.crop.sa01,'ind','test2','values')

pdf("Figure 1.5.pdf",height=5,width=10)
ggplot(floristic.crop.sa01,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(ind)))+
  geom_point(position=position_dodge(width=0.75),aes(group=ind,color=as.factor(ind)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.crop.w,aes(x=ind,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~test2,scales = "free_y")+ labs(x='Different Weeds',y='Number of present quadrats in each field')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("Figure 1.6.pdf",height=5,width=8)
ggplot(floristic.crop.sa01,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(test2)))+
  geom_point(position=position_dodge(width=0.75),aes(group=test2,color=as.factor(test2)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.crop.t,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~test2,scales = "free_y")+ 
  labs(x='Different Weeds',y='Number of present quadrats in each field')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

###### weed occurrence in edge fields ####

floristic.edge.bh<-subset(floristic.edge.sa,floristic.edge.sa$Group.3=="Bh")
floristic.edge.bsh<-subset(floristic.edge.sa,floristic.edge.sa$Group.3=="Bsh")

floristic.edge.bh1<-stack(floristic.edge.bh[,-1:-3])
floristic.edge.bh[,1:3]
test<-as.data.frame(rbind(floristic.edge.bh[,1:3],floristic.edge.bh[,1:3],floristic.edge.bh[,1:3],floristic.edge.bh[,1:3]
                          ,floristic.edge.bh[,1:3],floristic.edge.bh[,1:3],floristic.edge.bh[,1:3],floristic.edge.bh[,1:3]))
floristic.edge.bh2<-cbind(test,floristic.edge.bh1)
sta.edge.w <- ONE_Tukey_HSD1(floristic.edge.bh2,'Group.2','ind','values')

floristic.edge.bsh1<-stack(floristic.edge.bsh[,-1:-3])
test1<-as.data.frame(rbind(floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3]
                          ,floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3],floristic.edge.bsh[,1:3]))
floristic.edge.bsh2<-cbind(test1,floristic.edge.bsh1)
sta.edge.w1 <- ONE_Tukey_HSD1(floristic.edge.bsh2,'Group.2','ind','values')

pdf("Figure 1.2. floristic coverage at borders with hedges.pdf",height=5,width=10)
ggplot(floristic.edge.bh2,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(ind)))+
  geom_point(position=position_dodge(width=0.75),aes(group=ind,color=as.factor(ind)))+ 
  ylim(0,20)+
  geom_text(data=sta.edge.w,aes(x=ind,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~Group.2,scales = "free_y")+ labs(x='',y='Weed coverage in edge field with hedges(%)')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("Figure 1.3. floristic coverage at borders without hedges.pdf",height=5,width=10)
ggplot(floristic.edge.bsh2,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(ind)))+
  geom_point(position=position_dodge(width=0.75),aes(group=ind,color=as.factor(ind)))+ 
  ylim(0,20)+
  geom_text(data=sta.edge.w1,aes(x=ind,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~Group.2,scales = "free_y")+ labs(x='',y='Weed coverage in edge field without hedges(%)')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("202209 Figure 1.1.pdf",height=5,width=8)
ggplot(floristic.crop.s1,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(test2)))+
  geom_point(position=position_dodge(width=0.75),aes(group=test2,color=as.factor(test2)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.crop.t,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~test2,scales = "free_y")+ 
  labs(x='Different Weeds',y='Weed coverage in crop field centers')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("202209 Figure 1.3.pdf",height=5,width=8)
ggplot(floristic.edge.bh2,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(Group.2)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.2,color=as.factor(Group.2)))+ 
  ylim(0,20)+
  geom_text(data=sta.edge.w,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~Group.2,scales = "free_y")+ 
  labs(x='',y='Weed coverage in fields with hedgerow edges(%)')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

floristic.edge.bh2.c<-subset(floristic.edge.bh2,floristic.edge.bh2$Group.2=="conventional")
summary(aov(values~ind,floristic.edge.bh2.c))

floristic.edge.bh2.o<-subset(floristic.edge.bh2,floristic.edge.bh2$Group.2=="organic")
summary(aov(values~ind,floristic.edge.bh2.o))

summary(aov(values~ind+Group.2,floristic.edge.bh2))
chisq.test()

pdf("202209 Figure 1.2.pdf",height=5,width=8)
ggplot(floristic.edge.bsh2,aes(x=ind,y=values))+geom_boxplot(aes(fill=as.factor(Group.2)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.2,color=as.factor(Group.2)))+ 
  ylim(0,20)+
  geom_text(data=sta.edge.w1,aes(x=ind,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~Group.2,scales = "free_y")+ 
  labs(x='',y='Weed coverage in fields with grassy edges(%)')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

floristic.edge.bsh2.c<-subset(floristic.edge.bsh2,floristic.edge.bsh2$Group.2=="conventional")
summary(aov(values~ind,floristic.edge.bsh2.c))

floristic.edge.bsh2.o<-subset(floristic.edge.bsh2,floristic.edge.bsh2$Group.2=="organic")
summary(aov(values~ind,floristic.edge.bsh2.o))



###### 2.farmer's perception ####
sta.score.t <- ONE_Tukey_HSD1(perception.s2,'weed','type','values')
sta.score.w <- ONE_Tukey_HSD1(perception.s2,'type','weed','values')

pdf("Figure 2.2.pdf",height=5,width=10)
ggplot(perception.s2,aes(x=weed,y=values))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  geom_text(data=sta.score.w,aes(x=weed,y=mean+1.3*std,label=Letters))+
  ylim(0,10)+
  facet_wrap(.~type,scales = "free_y")+ labs(x='Different Weeds',y='Scores for each weed based on farmers')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf("Figure 2.3.pdf",height=5,width=8)
ggplot(perception.s2,aes(x=weed,y=values))+geom_boxplot(aes(fill=as.factor(type)))+
  geom_point(position=position_dodge(width=0.75),aes(group=type,color=as.factor(type)))+ 
  geom_text(data=sta.score.t,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='Different Weeds',y='Scores for each weed based on farmers')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()


###### 3.1 Weed root mycobiota sequence cluster richness #####
sta.weed <- ONE_Tukey_HSD1(richness.weed,'farm','weed','otu.richness.weed')

pdf("Figure 3.1.pdf",height=5,width=6)
ggplot(richness.weed,aes(x=weed,y=otu.richness.weed))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  #ylim(0,1.8)+
  geom_text(data=sta.weed,aes(x=weed,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~test2,scales = "free_y")+ 
  labs(x='Different Weeds',y='Weed root mycobiota sequence cluster richness')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))
dev.off()

V.sativa.r<-subset(richness.weed,richness.weed$weed=="Vicia sativa")
mean(V.sativa.r$otu.richness.weed)

L.pur.r<-subset(richness.weed,richness.weed$weed=="Lamium purpureum")
mean(L.pur.r$otu.richness.weed)

Mat.p.r<-subset(richness.weed,richness.weed$weed=="Matricaria sp")
mean(Mat.p.r$otu.richness.weed)

V.per.r<-subset(richness.weed,richness.weed$weed=="Veronica persica")
mean(V.per.r$otu.richness.weed)

(138.75-94.6)/94.6
(135.3333-94.6)/94.6
(141.4286-94.6)/94.6

library(ggplot2)
library(gridExtra)
###### 3.2 Weed root each phylum diversity ####
sta.weed.a <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.richness.asco.r')
sta.weed.b <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.richness.basi.r')
sta.weed.c <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.richness.chyt.r')
sta.weed.g <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.richness.glom.r')
sta.weed.z <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.richness.zygo.r')

pdf("Figure S3 Richness in weed root phylum v.s each weed.pdf",height=5,width=28)
a<-ggplot(data=richness.phylum,aes(x=weed,y=otu.richness.asco.r))+ 
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,80)+
  geom_text(data=sta.weed.a,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster richness of weed root Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

b<-ggplot(data=richness.phylum,aes(x=weed,y=otu.richness.basi.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,80)+
  geom_text(data=sta.weed.b,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster richness of weed root Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

c<-ggplot(data=richness.phylum,aes(x=weed,y=otu.richness.chyt.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,80)+
  geom_text(data=sta.weed.c,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster richness of weed root Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

g<-ggplot(data=richness.phylum,aes(x=weed,y=otu.richness.glom.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,80)+
  geom_text(data=sta.weed.g,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster richness of weed root Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

z<-ggplot(data=richness.phylum,aes(x=weed,y=otu.richness.zygo.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,80)+
  geom_text(data=sta.weed.z,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster richness of weed root Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

grid.arrange(a,b,c,g,z,ncol=5)
dev.off()

###### 3.2 Weed root each phylum abundance ####
sta.weed.aa <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.abundance.asco.r')
sta.weed.ba <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.abundance.basi.r')
sta.weed.ca <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.abundance.chyt.r')
sta.weed.ga <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.abundance.glom.r')
sta.weed.za <- ONE_Tukey_HSD1(richness.phylum,'farm','weed','otu.abundance.zygo.r')

pdf("Figure S6.2 Abundance in weed root phylum v.s each weed.pdf",height=5,width=28)
aa<-ggplot(data=richness.phylum,aes(x=weed,y=otu.abundance.asco.r))+ 
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,6000)+
  geom_text(data=sta.weed.aa,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster abundance of weed root Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

ba<-ggplot(data=richness.phylum,aes(x=weed,y=otu.abundance.basi.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,6000)+
  geom_text(data=sta.weed.ba,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster abundance of weed root Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

ca<-ggplot(data=richness.phylum,aes(x=weed,y=otu.abundance.chyt.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,6000)+
  geom_text(data=sta.weed.ca,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster abundance of weed root Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

ga<-ggplot(data=richness.phylum,aes(x=weed,y=otu.abundance.glom.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,6000)+
  geom_text(data=sta.weed.ga,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster abundance of weed root Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

za<-ggplot(data=richness.phylum,aes(x=weed,y=otu.abundance.zygo.r))+
  geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(weed)))+
  ylim(0,6000)+
  geom_text(data=sta.weed.za,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Sequence cluster abundance of weed root Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 90))

grid.arrange(aa,ba,ca,ga,za,ncol=5)
dev.off()

###### 3.3 Weed root pathogenic or symbiotic sequence cluster richness #####
sta.path.w <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.richness.path')
sta.sym.w <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.richness.sym')
sta.sap.w <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.richness.sap')

pdf("Figure 3.2 FUNGuild of weed root mycobiota.pdf",height=5,width=18)
path<-ggplot(richness.fun,aes(x=weed,y=otu.richness.path))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,35)+
  geom_text(data=sta.path.w,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster richness in Pathotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

sym<-ggplot(richness.fun,aes(x=weed,y=otu.richness.sym))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,35)+
  geom_text(data=sta.sym.w,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster richness in Symbiotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

sap<-ggplot(richness.fun,aes(x=weed,y=otu.richness.sap))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,35)+
  geom_text(data=sta.sap.w,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster richness in Saprotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

grid.arrange(path,sym,sap,ncol=3)
dev.off()

###### 3.3 Weed root pathogenic or symbiotic sequence cluster abundance #####
sta.path.wa <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.abundance.path')
sta.sym.wa <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.abundance.sym')
sta.sap.wa <- ONE_Tukey_HSD1(richness.fun,'farm','weed','otu.abundance.sap')

pdf("Figure S6.2 FUNGuild of weed root mycobiota abundance.pdf",height=5,width=18)
path.a<-ggplot(richness.fun,aes(x=weed,y=otu.abundance.path))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,3000)+
  geom_text(data=sta.path.wa,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster abundance in Pathotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

sym.a<-ggplot(richness.fun,aes(x=weed,y=otu.abundance.sym))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,3000)+
  geom_text(data=sta.sym.wa,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster abundance in Symbiotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

sap.a<-ggplot(richness.fun,aes(x=weed,y=otu.abundance.sap))+geom_boxplot(aes(fill=as.factor(weed)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed,color=as.factor(weed)))+ 
  ylim(0,6500)+
  geom_text(data=sta.sap.wa,aes(x=weed,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Weed root mycobiota sequence cluster abundance in Saprotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

grid.arrange(path.a,sym.a,sap.a,ncol=3)
dev.off()

###### 4.1 Transferred phylum by weed to wheat #####
sta.share.a <- ONE_Tukey_HSD1(richness.share.phy,'condition','weed.neighbour','richness.a')
sta.share.b <- ONE_Tukey_HSD1(richness.share.phy,'condition','weed.neighbour','richness.b')
sta.share.c <- ONE_Tukey_HSD1(richness.share.phy,'condition','weed.neighbour','richness.c')
sta.share.g <- ONE_Tukey_HSD1(richness.share.phy,'condition','weed.neighbour','richness.g')
sta.share.z <- ONE_Tukey_HSD1(richness.share.phy,'condition','weed.neighbour','richness.z')

s.a<-ggplot(richness.share.phy,aes(x=weed.neighbour,y=richness.a))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,13)+
  geom_text(data=sta.share.a,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.b<-ggplot(richness.share.phy,aes(x=weed.neighbour,y=richness.b))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,13)+
  geom_text(data=sta.share.b,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.c<-ggplot(richness.share.phy,aes(x=weed.neighbour,y=richness.c))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,13)+
  geom_text(data=sta.share.c,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.g<-ggplot(richness.share.phy,aes(x=weed.neighbour,y=richness.g))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,13)+
  #geom_text(data=sta.share.g,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.z<-ggplot(richness.share.phy,aes(x=weed.neighbour,y=richness.z))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,13)+
  #geom_text(data=sta.share.z,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

pdf("Figure S4. Shared weed root mycobiota richness in Phylum.pdf",height=5,width=28)
grid.arrange(s.a,s.b,s.c,s.g,s.z,ncol=5)
dev.off()

###### 4.2 Transferred FUNGuilds by weed to wheat #####
sta.share.num <- ONE_Tukey_HSD1(richness.share.fun,'condition','weed.neighbour','share.num')
sta.share.path <- ONE_Tukey_HSD1(richness.share.fun,'condition','weed.neighbour','richness.path')
sta.share.sym <- ONE_Tukey_HSD1(richness.share.fun,'condition','weed.neighbour','richness.sym')
sta.share.sap <- ONE_Tukey_HSD1(richness.share.fun,'condition','weed.neighbour','richness.sap')

s.num<-ggplot(richness.share.fun,aes(x=weed.neighbour,y=share.num))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,50)+
  geom_text(data=sta.share.num,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.path<-ggplot(richness.share.fun,aes(x=weed.neighbour,y=richness.path))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,50)+
  geom_text(data=sta.share.path,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Pathotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.sym<-ggplot(richness.share.fun,aes(x=weed.neighbour,y=richness.sym))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,8)+
  geom_text(data=sta.share.sym,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Symbiotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

s.sap<-ggplot(richness.share.fun,aes(x=weed.neighbour,y=richness.sap))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  ylim(0,8)+
  geom_text(data=sta.share.sap,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  labs(x='',y='Shared weed root mycobiota sequence cluster richness in Saprotroph')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

pdf("Figure 4 FUNGuild of shared weed root mycobiota.pdf",height=9,width=12)
grid.arrange(s.num,s.path,s.sym,s.sap,ncol=2)
dev.off()


###### 4.weed-wheat biomass competition, growth promotion ####
com.ind2<-subset(com.ind1,com.ind1$com2>(-1)&com.ind1$weed.neighbour!="Galium aparine")
com.ind2$cd<-(-com.ind2$com1)

sta.com.1 <- ONE_Tukey_HSD1(com.ind2,'type','weed.neighbour','cd')
sta.com.2 <- ONE_Tukey_HSD1(com.ind2,'type','weed.neighbour','com2')

com1<-ggplot(com.ind2,aes(x=weed.neighbour,y=-com1))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  geom_text(data=sta.com.1,aes(x=weed.neighbour,y=mean+2.3*std,label=Letters))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x='',y='Growth promotion effect of each weed on wheat')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

com2<-ggplot(com.ind2,aes(x=weed.neighbour,y=com2))+geom_boxplot(aes(fill=as.factor(weed.neighbour)))+
  geom_point(position=position_dodge(width=0.75),aes(group=weed.neighbour,color=as.factor(weed.neighbour)))+ 
  geom_text(data=sta.com.2,aes(x=weed.neighbour,y=mean+1.3*std,label=Letters))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x='',y='Competition index of each weed with single wheat')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90))

pdf("Figure 4.2 weed and wheat competition.pdf",height=5,width=12)
grid.arrange(com1,com2,ncol=2)
dev.off()

pdf("Figure 5 wheat growth promotion by weed.pdf",height=5,width=7)
com1
dev.off()

V.per.p<-subset(com.ind2,com.ind2$weed.neighbour=="Veronica persica")
mean(V.per.p$Wheat.weight)

Pa.r.p<-subset(com.ind2,com.ind2$weed.neighbour=="Papaver rhoeas")
mean(Pa.r.p$Wheat.weight)

Po.a.p<-subset(com.ind2,com.ind2$weed.neighbour=="Poa annua")
mean(Po.a.p$Wheat.weight)

(0.4872-0.8376778)/0.8376778
(0.51059-0.8376778)/0.8376778

###### new method to do statistics from Cendrine ####
library(agricolae)

model<-aov(Wheat.weight~com2,data=com.ind1)

comparison<- LSD.test(model,"weed.neighbour",alpha=0.05,group=TRUE)

###### identity effect in the field ####
#### loading the packages 
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

mod<-glmer.nb(otu.richness.weed~weed+(1|field), data=richness.weed)
#mod<-lmer(otu.richness.weed~weed+(1|field), data=richness.weed)

summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.asco.r~weed+(1|field), data=richness.phylum)
summary(mod)
Anova(mod)
AIC(mod)
#mod<-lmer(otu.richness.asco.r~weed+(1|field), data=richness.phylum)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.basi.r~weed+(1|field), data=richness.phylum)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.chyt.r~weed+(1|field), data=richness.phylum)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.glom.r~weed+(1|field), data=richness.phylum)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.zygo.r~weed+(1|field), data=richness.phylum)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals


mod<-glmer.nb(otu.richness.path~weed+(1|field), data=richness.fun)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.sym~weed+(1|field), data=richness.fun)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(otu.richness.sap~weed+(1|field), data=richness.fun)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

mod<-glmer.nb(values~weed+(1|field), data=perception.s2.o)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

perception.s2.c<-subset(perception.s2,perception.s2$type=="conventional")
mod<-glmer.nb(values~weed+(1|field), data=perception.s2.c)
summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

###### identity effect in the control experiment ####

lm1<-lm(com1~weed.neighbour,com.ind2)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(share.num~weed.neighbour,richness.share.fun)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.a~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.b~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.c~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.g~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.z~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.path~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.sym~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)

lm1<-lm(richness.sap~weed.neighbour,richness.share.phy)
anova(lm1)
summary(lm1)
AIC(lm1)



