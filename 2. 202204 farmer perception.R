################################################################
######NGS data analysis for weed root mycobiota in fields#######
knowledge<-read.table("farmers_knowledge.txt", header=T, row.names=1, sep="\t")
perception<-read.table("farmers_perception.txt", header=T, row.names=1, sep="\t")  

perception.f<-subset(perception,perception$Diplome!="NONE")
#####Farmers' knowledge to recognize weed species ####
ggplot(knowledge,aes(x=type,y=JNC))+geom_boxplot(aes(fill=as.factor(type)))+
  geom_point(aes(color=as.factor(type)),size=2)+
  theme(axis.text.x = element_text(angle = 60))+ 
  ylab("Number of farmeres do not know at least one weed species")

pdf("Figure 8. farmers perception.pdf",height=5,width=7)
ggplot(knowledge,aes(x=type,y=JNC.p*100))+geom_boxplot(aes(fill=as.factor(type)))+
  geom_point(aes(color=as.factor(type)),size=2)+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 60))+ 
  xlab("")+ylab("Percentage of farmeres do not know at least one weed species (%)")
dev.off()

knowledge.c<-subset(knowledge,knowledge$type=="conventional")
knowledge.o<-subset(knowledge,knowledge$type=="organic")

mean(knowledge.c$JNC.p)
mean(knowledge.o$JNC.p)

sd(knowledge.c$JNC.p)
sd(knowledge.o$JNC.p)

ggplot(knowledge,aes(x=type,y=JNC.p))+geom_boxplot(aes(fill=as.factor(type)))+
  geom_point(aes(color=as.factor(type)),size=2)+
  theme(axis.text.x = element_text(angle = 60))+ 
  ylab("Percentage of farmeres do not know at least one weed species")

summary(aov(JNC.p~type,knowledge))

chi<-as.data.frame(cbind(knowledge[1:8,]$JNC,knowledge[9:16,]$JNC))
chisq.test(chi$V1,chi$V2)
chisq.test(knowledge[1:8,]$JNC.p,knowledge[9:16,]$JNC.p,simulate.p.value = TRUE)
fisher.test(knowledge[1:8,]$JNC,knowledge[9:16,]$JNC)
fisher.test(knowledge[1:8,]$JNC,knowledge[9:16,]$JNC)

pdf("Figure S3. farmers knowledge about eight species.pdf",height=5,width=7)
ggplot(knowledge[1:8,], aes(x = weed, y=JNC.p)) +
  geom_col(fill = "red", color = "white") +
  geom_col(data = knowledge[9:16,],aes(x = weed, y = -JNC.p))+
  xlab("") + ylab("Percentage of farmeres do not know this weed species")+
  theme(axis.text.x = element_text(angle = 60))
dev.off()

summary(aov(values~ind+Group.2,floristic.edge.bh2))
chisq.test()

knowledge$all<-knowledge$JNC+knowledge$NP
summary.aov(lm(all~type,data=knowledge))
t.test(knowledge[1:8,]$all,knowledge[9:16,]$all,paired=TRUE)
t.test(knowledge[1:8,]$JNC,knowledge[9:16,]$JNC,paired=TRUE)
t.test(knowledge[1:8,]$NP,knowledge[9:16,]$NP,paired=TRUE)
c<-c(0.125,0.8125,0.3125,0.75,0.375,0.5625,0.3125,0.0625)
o<-c(1,1,1,1,1,0.7857,0.8571,0.7143)
t.test(c,o,paired=TRUE)
chisq.test(c,o)


ggplot(knowledge[1:8,], aes(x = weed, y=NP.p)) +
  geom_col(fill = "red", color = "white") +
  geom_col(data = knowledge[9:16,],aes(x = weed, y = -JNC.p))+
  xlab("") + ylab("Percentage of farmeres do not know this weed species")+
  theme(axis.text.x = element_text(angle = 60))

##### Troublesome of weed species based on farmers' view #####
perception.s<-perception[,c("Lamium_purpureum","Matricaria_sp","Poa_annua","Papaver_rhoeas",
                            "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_Sativa")]
perception.s1<-cbind(perception[,1:7],perception.s)

test.p<-stack(perception.s)
test.p1<-rep(perception[,1],8)
test.p2<-rep(perception[,3],8)
test.p3<-rep(perception[,4],8)

perception.s1<-cbind(test.p1,test.p2,test.p3,test.p)
write.table(perception.s1,file="perception.s1.txt",sep="\t")

perception.s2<-read.table("perception.s2.txt", header=T, row.names=1, sep="\t")  

perception.s2.c<-subset(perception.s2,perception.s2$type=="conventional")
perception.s2.o<-subset(perception.s2,perception.s2$type=="organic")
summary(aov(values~weed,perception.s2.c))
summary(aov(values~weed,perception.s2.o))


