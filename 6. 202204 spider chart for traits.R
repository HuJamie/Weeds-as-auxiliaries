################################################
#### Create Spider chart for each functions ####

library(fmsb)

## https://stackoverflow.com/questions/47644156/how-to-measure-the-area-of-a-polygon-in-ggplot2
## calculate the spider area

# Create data
set.seed(99)
data<-unstack(multi_std.b,value~variable)

#data1<-cbind(lab.mean[,1],data)
data1<-t(df_short)
#data2<-as.data.frame(data2[,-(8:9)])
colnames(data1)=c("Lamium purpureum","Matricaria chamomilla","Papaver rhoeas","Poa annua",
                  "Poa trivialis","Trifolium repens","Veronica persica","Viscia sativa")
data2<-t(data1)


# To use the fmsb package, I have to add 2 lines to the data frame: the max and min of each topic to show on the plot!
data.r=as.data.frame(rbind(rep(1.2,9) , rep(0,9) , df_short))

#==================
# Plot 1: Default radar chart proposed by the library:
radarchart(data.r)

radarchart(data.r[c(1,2,3),])
radarchart(data.r[c(1,2,4),])
radarchart(data.r[c(1,2,5),])
radarchart(data.r[c(1,2,6),])
radarchart(data.r[c(1,2,7),])
radarchart(data.r[c(1,2,8),])
radarchart(data.r[c(1,2,9),])
radarchart(data.r[c(1,2,10),])


pdf("Radar chart.pdf", height=10, width=14)

colors_border=c( rgb(0.3,0.6,0.1,0.9),rgb(0.8,0.2,0.5,0.9),rgb(0.7,0.5,0.1,0.9),rgb(0.2,0.5,0.5,0.9),rgb(0.1,0.8,0.5,0.9),
                 rgb(0.2,0.1,0.5,0.8),rgb(0.4,0.1,0.2,0.9),rgb(0.2,0.5,0.1,0.9),rgb(0.1,0.1,0.7,0.9),rgb(0.5,0.1,0.5,0.9))
colors_in=c( rgb(0.3,0.6,0.1,0.1),rgb(0.8,0.2,0.5,0.1),rgb(0.7,0.5,0.1,0.1),rgb(0.2,0.5,0.5,0.1),rgb(0.1,0.8,0.5,0.1),
             rgb(0.2,0.1,0.5,0.1),rgb(0.4,0.1,0.2,0.1),rgb(0.2,0.5,0.1,0.1),rgb(0.1,0.1,0.7,0.1),rgb(0.5,0.1,0.5,0.1))
radarchart(data.r, axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=1,
           #custom labels
           vlcex=1)
legend(x=1.2, y=1, legend = c("Lamium purpureum","Matricaria chamomilla","Papaver rhoeas","Poa annua",
                              "Poa trivialis","Trifolium repens","Veronica persica","Viscia sativa"), 
       bty = "n", pch=20 , col=colors_border, text.col = "grey", cex=1.2, pt.cex=3)

dev.off()



pdf("2023 Spider graph.pdf", height=8, width=20)

par(mfrow=c(2,4))
##1
colors_border=c( rgb(0.3,0.6,0.1,0.9))
colors_in=c( rgb(0.3,0.6,0.1,0.4))
radarchart(data.r[c(1,2,3),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##2
colors_border=c( rgb(0.8,0.2,0.5,0.9) )
colors_in=c( rgb(0.8,0.2,0.5,0.4) )
radarchart(data.r[c(1,2,4),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##3
colors_border=c(rgb(0.7,0.5,0.1,0.9))
colors_in=c(rgb(0.7,0.5,0.1,0.4))
radarchart(data.r[c(1,2,5),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##4
colors_border=c(rgb(0.2,0.5,0.5,0.9))
colors_in=c(rgb(0.2,0.5,0.5,0.4))
radarchart(data.r[c(1,2,6),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##5
colors_border=c(rgb(0.1,0.8,0.5,0.9))
colors_in=c(rgb(0.1,0.8,0.5,0.4))
radarchart(data.r[c(1,2,7),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##6
colors_border=c(rgb(0.2,0.1,0.5,0.8))
colors_in=c(rgb(0.2,0.1,0.5,0.4))
radarchart(data.r[c(1,2,8),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##7
colors_border=c(rgb(0.4,0.1,0.2,0.9))
colors_in=c(rgb(0.4,0.1,0.2,0.4))
radarchart(data.r[c(1,2,9),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
##8
colors_border=c(rgb(0.2,0.5,0.1,0.9))
colors_in=c(rgb(0.2,0.5,0.1,0.4))
radarchart(data.r[c(1,2,10),], axistype=1 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.2,0.3), cglwd=0.8,
           #custom labels
           vlcex=0.8)
dev.off()

############################################################################################################
#### Calculate the area and evenness ####
library(vegan)
shannon.m<- diversity(df_short,index="shannon") 
richness.m<-rowSums(df_short>0)
evenness.m<- as.data.frame(shannon.m/log(richness.m))
colnames(evenness.m)<-"area"
#area<-as.data.frame(c(0.2716312,0.6431483,0.1498839,0.7998977,0.6787655,0.9635457,1.100601,0.01346223))
area<as.data.frame(c(0.3559474,0.54689,0.1575,0.8283652,0.7913697,0.8884152,1.144626,0.249483))
colnames(area)<-"area"

summa<-rbind(evenness.m,area)
summa$weed<-rownames(evenness.m)

tes<-c("even","even","even","even","even","even","even","even",
       "area","area","area","area","area","area","area","area")

summa<-cbind(summa,tes)

library(dplyr)
library(ggplot2)

pdf("2023 Spider Figure summary.pdf",height=5,width=9)

  ggplot(summa,aes(x=weed, y=area, fill=weed))+
  geom_col(position = "dodge")+
  #scale_y_continuous(sec.axis = sec_axis(~ . * 1, name = "Evenness of selected traits"))+
  labs(y="Relative area size of spider graph")+
    facet_wrap(.~tes,scales = "free_y")+
    ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle=90))

dev.off()

#### For Lamium purpureum
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,1]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Matricaria sp.
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,2]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))


#### For Papaver rhoeas
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,3]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Poa annua
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,4]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Poa trivialis
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,5]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Trifolium repens
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,6]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Veronica persica
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,7]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### For Vicia sativa
df <- structure(
  list(
    V1 = 1:6,
    V2 = t(df_short)[,8]
  ),
  .Names = c("V1", "V2"),
  class = "data.frame",
  row.names = c(NA, -6L)) 

areas <- df %>% 
  setNames(c("variable", "value")) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6) * (2*pi),
         # change 1/n to number of variables
         area = value*nextval*sin(angle)/2)

areas %>% summarise(total = sum(area))

#### A bit advanced script ####
df_i<-stack(as.data.frame(t(df_short)))

df <- expand.grid(var = 1:6, grp = c("Lamium_purpureum","Matricaria_sp.","Poa_annua","Papaver_rhoeas",
                                     "Poa_trivialis","Trifolium_repens","Veronica_persica","Vicia_sativa")) %>% 
  mutate(value = df_i) %>% 
  
  group_by(grp) %>% 
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/6)*(2*pi),
         area = value*nextval*sin(angle)/2) %>% 
  mutate(total = sum(area)) 






