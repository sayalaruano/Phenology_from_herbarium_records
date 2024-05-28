rm(list=ls(all=TRUE))

library(dplyr)
library(tidyverse)
library(plotrix)
library(cowplot)
library(circular)
library(bpnreg)
library(bpDir)
library(CircMLE)

####read data herbarium ##########
setwd("C:/Users/Jenny/Documents/UDLA/propuestas DGIV/2021-10 Fenologia/paper fenologia/Phenology v5_15-5-2024/scripts and data")

all.gbif.pred_trop<- read.csv("S2-data herbarium records.csv",header=T)
dim(all.gbif.pred_trop)
colnames(all.gbif.pred_trop)

spp_tbl <- all.gbif.pred_trop %>% group_by(acceptedScientificName) %>% 
  summarise(total_count=n())
#write.csv(spp_tbl,"summ_spp.csv")

fam_tbl <- all.gbif.pred_trop %>% group_by(family, acceptedScientificName) %>% 
  summarise(total_count=n())
#write.csv(fam_tbl,"summ_fam.csv")

##in the table records by families and species there are 429 species (in place of 427 uniques)
##becuase 2 species have duplicates since they have been assigned to 2 different families, possibly
##due to changes in taxonomy of the species over time
# these cases are: 
## Cordia hebeclada	n=4, assigned to the Cordiaceae family (current)
## Cordia hebeclada	n=24, assigned to the Boraginaceae familiy 
## Tournefortia fuliginosa	n=37 assigned to Heliotropiaceae family (current)
## Tournefortia fuliginosa	n=20 assigned to Boraginaceae family

#### what about flowering data#################################

all.gbif.pred_trop<-merge(all.gbif.pred_trop,spp_tbl,by="acceptedScientificName" )

for (i in 1:nrow(all.gbif.pred_trop[1])){
      if (all.gbif.pred_trop$total_count[i] >=500) {all.gbif.pred_trop$nrecords_interval[i]<-">500"} 
        else{
            if (all.gbif.pred_trop$total_count[i] >=100 & 
                all.gbif.pred_trop$total_count[i] < 500) {all.gbif.pred_trop$nrecords_interval[i]<-"100-499"}
            else {
                if (all.gbif.pred_trop$total_count[i] >=50 & 
                    all.gbif.pred_trop$total_count[i] < 100) {all.gbif.pred_trop$nrecords_interval[i]<-"50-99"}
                else {
                    if (all.gbif.pred_trop$total_count[i] >=20 & 
                        all.gbif.pred_trop$total_count[i] < 50) {all.gbif.pred_trop$nrecords_interval[i]<-"20-49"}
                  else{ all.gbif.pred_trop$nrecords_interval[i]<-"<20"} 
                     }
                 }
            }
        }

##check for elevation and latitude##########################

a<-gsub(" ", "", all.gbif.pred_trop$elevation)
a<-sub("m.*", "",a)
all.gbif.pred_trop$elevation1 <-as.numeric(sub("-.*", "",a))
all.gbif.pred_trop$elevation2 <-as.numeric(sub(".*-", "",a)) 
all.gbif.pred_trop$mean_elevation<-(all.gbif.pred_trop$elevation1+all.gbif.pred_trop$elevation2)/2

summary(all.gbif.pred_trop$mean_elevation)

alt.data<-all.gbif.pred_trop %>%
    filter(mean_elevation>0,mean_elevation<=4000)

alt.data$nrecords_interval <- 
  factor(alt.data$nrecords_interval,                                    # Change ordering manually
  levels = c("<20", "20-49", "50-99", "100-499", ">500"))

ggplot(alt.data, aes(x=nrecords_interval,y=mean_elevation, fill=nrecords_interval))+
  geom_violin()+
  geom_boxplot(width=0.05, fill="white") + 
  labs(x="n records", y="Elevation (m asl)", colour="Elevation")+
  theme_half_open(12) +
  theme(legend.position="none")

spp_tbl_flw_alt <- alt.data %>%
  #filter(nrecords_interval==">500")%>%
  group_by(acceptedScientificName) %>% 
  summarise(total_count=n(),min_alt=min(mean_elevation, na.rm=T),
            max_alt=max(mean_elevation, na.rm=T),mean_alt=mean(mean_elevation, na.rm=T))

#write.csv(spp_tbl_flw_alt, "species_byelevation.csv")

## check distribution of n records by latitude
  
summary(all.gbif.pred_trop$decimalLatitude)

all.gbif.pred_trop$nrecords_interval <- 
  factor(all.gbif.pred_trop$nrecords_interval,                                    # Change ordering manually
         levels = c("<20", "20-49", "50-99", "100-499", ">500"))

ggplot(all.gbif.pred_trop, aes(x=nrecords_interval,y=decimalLatitude, fill=nrecords_interval))+
  geom_violin()+
  geom_boxplot(width=0.02, fill="white") + 
  labs(x="n records", y="Latitude (decimal degree)", colour="Elevation")+
  theme_half_open(12) +
  theme(legend.position="none")


####dataset with phenology###################

summary(all.gbif.pred_trop$year)
summary(all.gbif.pred_trop$month)
summary(all.gbif.pred_trop$day)

var.include <-c("data_source"  ,"gbifID" , "hcnqID" , "eventDate" , "year",                    
                "month" , "day" , "decimalLongitude",     
                "decimalLatitude" ,  "mean_elevation" , "country",
                "family" , "genus" , "species" , "acceptedScientificName",
                "total_count"   ,  "nrecords_interval"  , 
                "Flowering_pred_RF", "Flowering_pred_prob_RF" ,  
                "clusters")

col.index<-which(colnames(all.gbif.pred_trop) %in% var.include)

##create a dataset for further analysis with selected columns

linked.data<-all.gbif.pred_trop[,col.index]

##transform dates with day 0 to day 1 to avoid lossing data #####

summary(linked.data$day)

###create a day 1 for records with no day in date, 
linked.data$new.day<-linked.data$day
linked.data$new.day[linked.data$day == 0]<-1
linked.data$new.day[is.na(linked.data$day) == TRUE]<-1

summary(linked.data$new.day)

#format date and create start date to calculate number of days in the year and transform to circular scale

linked.data$startdate <- as.Date(paste("01/01/",linked.data$year,sep=""),"%d/%m/%Y")

linked.data$FlwDate <- as.Date (paste(linked.data$new.day,"/",
                                      linked.data$month,"/",
                                      linked.data$year,sep=""),"%d/%m/%Y" )

linked.data$NumDays  <- difftime(linked.data$FlwDate,linked.data$startdate, units="days")

linked.data$fld.date.grad<- as.numeric(linked.data$NumDays)*360/365
linked.data$fld.date.rad<- linked.data$fld.date.grad * pi / 180

linked.data$fld.date.grad<- circular(linked.data$fld.date.grad, units = "degrees",zero=circular(pi/2), rotation="clock")
linked.data$fld.date.rad<- circular(linked.data$fld.date.rad, units = "radians",zero=circular(pi/2), rotation="clock")

###check coverage of data by clusters and select only flowering records ##########

cluster_spp_flw <- linked.data %>% 
  group_by(clusters, acceptedScientificName,Flowering_pred_RF) %>% 
  summarise(cluster_spp_flw=n())
#write.csv(spp_tbl,"cluster_spp_flw.csv")

cluster_spp_flw_20<-cluster_spp_flw%>%
  filter(Flowering_pred_RF=="Yes" , cluster_spp_flw>=20 )

###run best models distribution for all species with more than 20 flowering records#######

linked.data.flw<-merge(linked.data, cluster_spp_flw_20, by=c("acceptedScientificName","clusters", "Flowering_pred_RF"))
linked.data.flw$spp.x.region<-paste(linked.data.flw$acceptedScientificName,"_",linked.data.flw$clusters,sep="" )

spp.reg.flw<-unique(linked.data.flw$spp.x.region)

#####100 runs per species to determine best model for each species x region combination 
best.model<- data.frame(NULL)
for (j in 1:100){
  for (i in 1:162){
    spp.i<- linked.data.flw[linked.data.flw$spp.x.region==spp.reg.flw[i],]
    dist.mod<- circ_mle(
      spp.i$fld.date.rad,
      criterion = "AIC",
      nchains = 5,
      BadStart = 10^9,
      niter = 5000,
      method = "BFGS",
      lambda.min = 0.25,
      exclude = NULL)

    best.model.i <- dist.mod$results[which(rownames(dist.mod$results) == dist.mod$bestmodel), ]
    best.model.i<-cbind(spp.reg.flw[i],dist.mod$bestmodel,paste("run_",i,sep=""),best.model.i )
    best.model <-rbind(best.model,best.model.i)
  }
  print(j)
}
write.csv(best.model,"bestmodel100.csv")

spp.cluster<-linked.data.flw[,c(1,2,27,28)]
spp.cluster<-unique(spp.cluster)

best_model_100<-merge(spp.cluster,best_model_100,by="spp.x.region")
best_model_100$distribution<-NA

best_model_100$distribution[which(best_model_100$bestmodel=="M1")]<-"Uniform"
best_model_100$distribution[which(best_model_100$bestmodel=="M2A")]<-"Unimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M2B")]<-"Unimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M2C")]<-"Unimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M3A")]<-"Bimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M3B")]<-"Bimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M4A")]<-"Bimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M4B")]<-"Bimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M5A")]<-"Bimodal"
best_model_100$distribution[which(best_model_100$bestmodel=="M5B")]<-"Bimodal"

as.factor(best_model_100$distribution)
summary(best_model_100$distribution)
unique(best_model_100$distribution)

model.consistency<-best_model_100%>% 
  group_by(clusters, acceptedScientificName,distribution,bestmodel) %>% 
  summarise(total_count=n())

n.models<-model.consistency %>%
  group_by(clusters, acceptedScientificName) %>% 
  summarise(n_model=n())
  
model.consistency<-merge(model.consistency,n.models, by=c("clusters" , "acceptedScientificName"))

model.consistency$dist_mod_consist<-NA

model.consistency$dist_mod_consist[which(model.consistency$total_count>=75 &
                                         model.consistency$n_model<=2)]<-"consistent"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=25 &
                                           model.consistency$n_model<=2)]<-"delete"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=75 & model.consistency$total_count>25 &
                                           model.consistency$n_model<=2)]<-"inconsistent"

model.consistency$dist_mod_consist[which(model.consistency$total_count>=75 &
                                           model.consistency$n_model==3)]<-"consistent"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=20 &
                                           model.consistency$n_model==3)]<-"delete"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=75 & model.consistency$total_count>20 &
                                           model.consistency$n_model==3)]<-"inconsistent"

model.consistency$dist_mod_consist[which(model.consistency$total_count>=75 &
                                           model.consistency$n_model==4)]<-"consistent"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=18 &
                                           model.consistency$n_model==4)]<-"delete"

model.consistency$dist_mod_consist[which(model.consistency$total_count<=75 & model.consistency$total_count>18 &
                                           model.consistency$n_model==4)]<-"inconsistent"

#write.csv(model.consistency,"model.consistency.csv",row.names = F)

main.dist<-model.consistency%>% 
  filter(dist_mod_consist=="consistent")

variable.dist<-model.consistency%>% 
  filter(dist_mod_consist=="inconsistent")

best_model_consist_spp<-merge(best_model_100,main.dist,
                             by=c("acceptedScientificName","clusters","bestmodel"))
best_model_variable_spp <- merge(best_model_100,variable.dist,
                                 by=c("acceptedScientificName","clusters","bestmodel"))

##create a dataset of 162 species x regions combinations and their distribution
dist.cons<-main.dist %>%
  group_by(clusters, acceptedScientificName, distribution) %>% 
  summarise(n_model=n())

dist.var<- variable.dist %>%
  group_by(clusters, acceptedScientificName, distribution) %>% 
  summarise(n_model=n())

n.dist.var<-dist.var %>%
  group_by(clusters, acceptedScientificName) %>% 
  summarise(n_dist=n())

dist.var<-merge(dist.var,n.dist.var,by=c("clusters", "acceptedScientificName"))

dist.cons<-rbind(dist.cons,
                 dist.var[which(dist.var$n_dist==1),c(1:2,4,3)])

n.dist.var$dist.cons<-NA
n.dist.var$dist.cons[which(n.dist.var$n_dist==1 )]<-"one.dist"
n.dist.var$dist.cons[which(n.dist.var$n_dist>1 )]<-"two.dist"

n.dist.var<-n.dist.var[which(n.dist.var$n_dist==2),c(1:2,4,3)]
colnames(n.dist.var)<-c("clusters", "acceptedScientificName" ,"distribution" , "n_model")

dist.cons<-rbind(dist.cons,n.dist.var)

#join data to linked data to work figures

linked.data.flw<-merge(linked.data.flw,dist.cons[,c(1:3)], by=c("clusters", "acceptedScientificName"))

## generate parameters data for figures 

param_consist_spp<-best_model_consist_spp %>% 
  distinct(acceptedScientificName, clusters, bestmodel, .keep_all = TRUE)

param_variable_spp<-best_model_variable_spp %>% 
  distinct(acceptedScientificName, clusters, bestmodel, .keep_all = TRUE)

#table summary models by species x region combinations #######
best_model_100$best.model.spp.reg<-paste(best_model_100$spp.x.region, "_" ,
                                         best_model_100$bestmodel ,sep="")

model.consistency$best.model.spp.reg<-paste(model.consistency$acceptedScientificName, "_" ,
                                            model.consistency$clusters, "_" ,
                                            model.consistency$bestmodel,sep="")

head(best_model_100$best.model.spp.reg)
head(model.consistency$best.model.spp.reg)

regions<-paste0("region",unique(model.consistency$clusters))

table.param.mod.spp.reg<-data.frame(NULL)
for (i in 1:length(regions)){ 
  
  param.data <- best_model_100%>%
    filter(clusters==i) 
  
  summ.data <- model.consistency %>%
    filter(clusters==i) 
  
  for (j in 1:nrow(summ.data)){ 
    #    for (j in 4:5){ 
    
    param.data.j <- param.data%>%
      filter(best.model.spp.reg==summ.data$best.model.spp.reg[j]) 
    
    mod.consis<-summ.data$total_count[j]
    
    param.data.jj<-filter(param.data.j,AIC==min(param.data$AIC))
    if (nrow(param.data.jj) != 1) {
      param.data.jj<-param.data.j[1,]
    }else{
      param.data.jj=param.data.jj
    }
    param.data.jj<-cbind(param.data.jj,mod.consis)
    table.param.mod.spp.reg<-rbind(table.param.mod.spp.reg,param.data.jj)
  }
}  

table.param.mod.spp.reg$day.flw.1<-(table.param.mod.spp.reg$q1*(180/pi))*365/360
table.param.mod.spp.reg$day.flw.2<-(table.param.mod.spp.reg$q2*(180/pi))*365/360

table.param.mod.spp.reg$flw.date.1<-as.Date(table.param.mod.spp.reg$day.flw.1,origin = "2023-01-01")
table.param.mod.spp.reg$flw.date.2<-as.Date(table.param.mod.spp.reg$day.flw.2,origin = "2023-01-01")

#write.csv(table.param.mod.spp.reg,"table_modcoeff_spp_reg.csv",row.names = F)


##Modify function for plot circMLE JO tu run with table generated by Jenny parameters (selected from 100 runs)####

plot_circMLE_JO<-function (data, table, model, bins, shrink, col, lwd, lty, col.points, sep, cex, main.t) 
{
  if (is.null(data) | missing(data)) 
    stop("Please provide input data vector")
  if (is.null(table) | missing(table)) 
    stop("Please provide the output list from the \"circ_mle\" function")
  if (missing(model)) 
    model = rownames(table)[1]
  if (length(model) != 1) 
    stop("Only 1 model can be specified")
  if (!all(model %in% rownames(table))) 
    stop("Model not specified correctly.")
  if (missing(bins)) 
    bins = 18
  else bins = bins
  if (length(bins) != 1 | !is.numeric(bins)) 
    stop("Must specify the number of bins as a numeric vector of length 1")
  if (missing(shrink)) 
    shrink = 1.5
  else shrink = shrink
  if (length(shrink) != 1 | !is.numeric(shrink)) 
    stop("Must specify \"shrink\" as a numeric vector of length 1")
  if (missing(col)) 
    col = c("grey", "red", "black", "black")
  else col = col
  if (missing(lwd)) 
    lwd = c(2, 2, 2)
  else lwd = lwd
  if (missing(lty)) 
    lty = c("solid", "dashed", "dashed")
  else lty = lty
  if (missing(col.points)) 
    col.point= "grey30"
  else col.point= col.points
  if (missing(sep)) 
    sep = 0.05
  else sep = sep
  if (missing(cex)) 
    cex = 1
  else cex = cex
  if (missing(main.t)) 
    main = ""
  else main = main.t
  data = check_data(data)
  params = circularp(data)
  rose.diag(data, bins = bins, col = col[1], shrink = shrink, 
            main=main, axes=F)
  axis.circular(at=circular(seq(0, 11*pi/6,pi/6)), 
                labels=c("Apr","Mar","Feb","Jan","Dec","Nov",
                         "Oct", "Sep","Aug","Jul","Jun","May"),cex=1)
  points(x=data, cex = cex, stack = T, start.sep=0, sep = sep, 
         shrink = shrink, bins = NULL, col = col.point, zero = pi/2)
  arrows.circular(mean.circular(data, na.rm = T), 
                  y = rho.circular(data, na.rm = T), 
                  col = col[2], lwd = lwd[1], lty = lty[1])
 
  cols<-c("q1","k1","lamda","q2","k2")
  
  col.index<-which(colnames(table) %in% cols)
  
  model.vector <- table[,col.index]
  
  #model.vector = table[which(rownames(table) == model), 3:7]
  q1 = suppressWarnings(as.circular(model.vector[1], 
                                    control.circular = params))
  q2 = suppressWarnings(as.circular(model.vector[4], 
                                    control.circular = params))
  k1 = as.numeric(model.vector[2])
  k2 = as.numeric(model.vector[5])
  l = as.numeric(model.vector[3])
  plot.function.circular(function(x) dmixedvonmises(x, prop=l,
                                                    as.circular(q1, control.circular = params), 
                                                    as.circular(q2, control.circular = params), 
                                                    k1, k2), 
                                                    add = T, 
                                                    lwd = lwd[2], 
                                                    shrink = shrink, 
                                                    lty = lty[2], 
                                                    col = col[3])
  
  if (any(c("M2A", "M2B", "M2C") == model)) {
    arrows.circular(q1, col = col[4], lwd = lwd[3], lty = lty[3])
  }
  if (any(c("M3A", "M3B", "M4A", "M4B", "M5A", "M5B") == model)) {
    arrows.circular(q1, col = col[4], lwd = lwd[3], lty = lty[3])
    arrows.circular(q2, col = col[4], lwd = lwd[3], lty = lty[3])
  }
}

###############by spp in region 1 #########
##Unimodal
spp.plot<- linked.data.flw %>%
  filter(clusters==1 & Flowering_pred_RF=="Yes", total_count>=20)

spp<-unique(spp.plot$acceptedScientificName )

spp.plot.unim<- spp.plot %>%
  filter(distribution == "Unimodal" )

spp.unim.r1<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:15){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r1[i] & clusters==1)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r1[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#fac228",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.unim.r1[i],"-R1"))
  
  savePlot(filename = paste("R1spp_plot_unimodal_",2,sep=""), # Name of the file to be saved
            type = c("tiff"), # Type of file to be saved
            device = dev.cur(), # Number of the device to be saved
            restoreConsole = TRUE)
  
  }

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(distribution == "Bimodal" )

spp.bim.r1<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 15:17){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r1[i] & clusters==1)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r1[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#fac228",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r1[i],"-R1"))
  
  savePlot(filename = paste("R1spp_plot_bimodal_",2,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
  
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r1<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,4))

for (i in 13:15){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r1[i] & clusters==1)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r1[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#fac228",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.unif.r1[i],"-R1"))
  
  savePlot(filename = paste("R1spp_plot_uniform_",2,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r1<-unique(spp.plot.incons$acceptedScientificName)

param_variable_spp$acceptedScientificName

windows(10,4)
par(mar = c(0,1,1,1),mfrow=c(2,3))

for (i in 1:2){
  length.spp<-length(param_variable_spp[param_variable_spp$acceptedScientificName==spp.incos.r1[1],1 ])
  # windows(10,2)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r1[i] & clusters==1)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r1[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#fac228",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r1[i],"-R1",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R1spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}

#########by spp in region 2 ##############
spp.plot<- linked.data.flw %>%
  filter(clusters==2 & Flowering_pred_RF=="Yes")

spp.plot<-merge(spp.plot,main.dist,by=c("acceptedScientificName","clusters"))

spp<-unique(spp.plot$acceptedScientificName )

##Unimodal

spp.plot.unim<- spp.plot %>%
  filter(Dist.mod.cons == "Unimodal" )

spp.unim.r2<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:6){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r2[i] & clusters==2)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r2[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#f57d15",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.unim.r2[i],"-R2"))
  
  savePlot(filename = paste("R2spp_plot_unimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(Dist.mod.cons == "Bimodal" )

spp.bim.r2<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:4){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r2[i] & clusters==2)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r2[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#f57d15",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r2[i],"-R2"))
  
  savePlot(filename = paste("R2spp_plot_bimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r2<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(2,5))

for (i in 1:1){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r2[i] & clusters==2)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r2[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#f57d15",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.unif.r2[i],"-R2"))
  
  savePlot(filename = paste("R2spp_plot_uniform_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r2<-unique(spp.plot.incons$acceptedScientificName)

for (i in 1:2){
  length.spp<-length(param_variable_spp[param_variable_spp$acceptedScientificName==spp.incos.r2[1],1 ])
  # windows(10,3)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r2[i] & clusters==2)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r2[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#f57d15",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r2[i],"-R2",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R2spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}

##############by spp in region 3 ###################
spp.plot<- linked.data.flw %>%
  filter(clusters==3 & Flowering_pred_RF=="Yes")

spp.plot<-merge(spp.plot,main.dist,by=c("acceptedScientificName","clusters"))

spp<-unique(spp.plot$acceptedScientificName )

spp.plot.unim<- spp.plot %>%
  filter(Dist.mod.cons == "Unimodal" )

spp.unim.r3<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:13){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r3[i] & clusters==3)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r3[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#d44842",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unim.r3[i],"-R3"))
  
  savePlot(filename = paste("R3spp_plot_unimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(Dist.mod.cons == "Bimodal" )

spp.bim.r3<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 3:11){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r3[i] & clusters==3)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r3[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#d44842",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r3[i],"-R3"))
  
  savePlot(filename = paste("R3spp_plot_bimodal_",2,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r3<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:1){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r3[i] & clusters==3)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r3[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#d44842",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.unif.r3[i],"-R3"))
  
  savePlot(filename = paste("R3spp_plot_uniform_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r3<-unique(spp.plot.incons$acceptedScientificName)

param_variable_spp_3<-param_variable_spp[param_variable_spp$clusters==3,]

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,4))

for (i in 1:3){
  length.spp<-length(param_variable_spp_3[param_variable_spp_3$acceptedScientificName==spp.incos.r3[i],1 ])
  # windows(10,3)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r3[i] & clusters==3)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r3[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#d44842",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r3[i],"-R3",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R3spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}

#### by spp in region 4 ##########
spp.plot<- linked.data.flw %>%
  filter(clusters==4 & Flowering_pred_RF=="Yes")

spp.plot<-merge(spp.plot,main.dist,by=c("acceptedScientificName","clusters"))

spp<-unique(spp.plot$acceptedScientificName )

spp.plot.unim<- spp.plot %>%
  filter(Dist.mod.cons == "Unimodal" )

spp.unim.r4<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:15){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r4[i] & clusters==4)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r4[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#9f2a63",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unim.r4[i],"-R4"))
  
  savePlot(filename = paste("R4spp_plot_unimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(Dist.mod.cons == "Bimodal" )

spp.bim.r4<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:15){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r4[i] & clusters==4)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r4[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#9f2a63",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r4[i],"-R4"))
  
  savePlot(filename = paste("R4spp_plot_bimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r4<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:14){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r4[i] & clusters==4)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r4[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#9f2a63",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unif.r4[i],"-R4"))
  
  savePlot(filename = paste("R4spp_plot_uniform_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r4<-unique(spp.plot.incons$acceptedScientificName)

param_variable_spp_4<-param_variable_spp[param_variable_spp$clusters==4,]

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,4))

for (i in 1:length(spp.incos.r4)){
  length.spp<-length(param_variable_spp_4[param_variable_spp_4$acceptedScientificName==spp.incos.r4[i],1 ])
  # windows(10,3)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r4[i] & clusters==4)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r4[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#9f2a63",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r4[i],"-R4",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R4spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}

#### by spp in region 5 ##########
spp.plot<- linked.data.flw %>%
  filter(clusters==5 & Flowering_pred_RF=="Yes")

spp.plot<-merge(spp.plot,main.dist,by=c("acceptedScientificName","clusters"))

spp<-unique(spp.plot$acceptedScientificName )

spp.plot.unim<- spp.plot %>%
  filter(Dist.mod.cons == "Unimodal" )

spp.unim.r5<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:9){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r5[i] & clusters==5)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r5[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#65156e",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unim.r5[i],"-R5"))
  
  savePlot(filename = paste("R5spp_plot_unimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(Dist.mod.cons == "Bimodal" )

spp.bim.r5<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:3){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r5[i] & clusters==5)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r5[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#65156e",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r5[i],"-R5"))
  
  savePlot(filename = paste("R5spp_plot_bimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r5<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:1){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r5[i] & clusters==5)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r5[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#65156e",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unif.r5[i],"-R5"))
  
  savePlot(filename = paste("R5spp_plot_uniform_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r5<-unique(spp.plot.incons$acceptedScientificName)

param_variable_spp_5<-param_variable_spp[param_variable_spp$clusters==5,]

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,4))

for (i in 1:length(spp.incos.r5)){
  length.spp<-length(param_variable_spp_5[param_variable_spp_5$acceptedScientificName==spp.incos.r5[i],1 ])
  # windows(10,3)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r5[i] & clusters==5)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r5[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#65156e",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r5[i],"-R5",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R5spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}


#### by spp in region 6  ###########
spp.plot<- linked.data.flw %>%
  filter(clusters==6 & Flowering_pred_RF=="Yes")

spp.plot<-merge(spp.plot,main.dist,by=c("acceptedScientificName","clusters"))

spp<-unique(spp.plot$acceptedScientificName )

spp.plot.unim<- spp.plot %>%
  filter(Dist.mod.cons == "Unimodal" )

spp.unim.r6<-unique(spp.plot.unim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:1){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unim.r6[i] & clusters==6)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unim[spp.plot.unim$acceptedScientificName==spp.unim.r6[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#280b53",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unim.r6[i],"-R6"))
  
  savePlot(filename = paste("R6spp_plot_unimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Bimodal
spp.plot.bim<- spp.plot %>%
  filter(Dist.mod.cons == "Bimodal" )

spp.bim.r6<-unique(spp.plot.bim$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:1){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.bim.r6[i] & clusters==6)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.bim[spp.plot.bim$acceptedScientificName==spp.bim.r6[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#280b53",sep=0.08, cex=1, shrink=1.4,
                  main.t=paste(spp.bim.r6[i],"-R6"))
  
  savePlot(filename = paste("R6spp_plot_bimodal_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

##Uniform
spp.plot.unif<- spp.plot %>%
  filter(Dist.mod.cons == "Uniform" )

spp.unif.r6<-unique(spp.plot.unif$acceptedScientificName)

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,5))

for (i in 1:6){
  param.data <- param_consist_spp%>%
    filter(acceptedScientificName==spp.unif.r6[i] & clusters==6)
  
  row.names(param.data)<-param.data$bestmodel
  param.data<-param.data[,-c(1:5)]
  
  spp.i<- spp.plot.unif[spp.plot.unif$acceptedScientificName==spp.unif.r6[i],]
  
  plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                  col.points="#280b53",sep=0.08, cex=0.8, shrink=1.4,
                  main.t=paste(spp.unif.r6[i],"-R6"))
  
  savePlot(filename = paste("R6spp_plot_uniform_",1,sep=""),    # Name of the file to be saved
           type = c("tiff"), # Type of file to be saved
           device = dev.cur(), # Number of the device to be saved
           restoreConsole = TRUE)
}

###Inconsistent

spp.plot.incons<- spp.plot %>%
  filter(Dist.mod.cons == "Inconsistent" )

spp.incos.r6<-unique(spp.plot.incons$acceptedScientificName)

param_variable_spp_6<-param_variable_spp[param_variable_spp$clusters==6,]

windows(60,30)
par(mar = c(0,1,1,1),mfrow=c(3,4))
for (i in 1:length(spp.incos.r6)){
  length.spp<-length(param_variable_spp_6[param_variable_spp_6$acceptedScientificName==spp.incos.r6[i],1 ])
  # windows(10,3)
  # par(mar = c(0,1,2,1),mfrow=c(1,length.spp))
  for(j in 1:length.spp){
    
    param.data<- param_variable_spp%>%
      filter(acceptedScientificName==spp.incos.r6[i] & clusters==6)
    BM<-param.data$bestmodel
    porc<-param.data$total_count
    param.data<-param.data[j,]
    row.names(param.data)<-param.data$bestmodel
    param.data<-param.data[,-c(1:5)]
    
    spp.i<- spp.plot.incons[spp.plot.incons$acceptedScientificName==spp.incos.r6[i],]
    
    
    plot_circMLE_JO(spp.i$fld.date.rad, param.data, bins=12,
                    col.points="#280b53",sep=0.08, cex=1, shrink=1.6,
                    main.t=paste(spp.incos.r6[i],"-R6",BM[j],"-",porc[j],"%"))
    
    savePlot(filename = paste("R6spp_plot_variable_",1,sep=""),    # Name of the file to be saved
             type = c("tiff"), # Type of file to be saved
             device = dev.cur(), # Number of the device to be saved
             restoreConsole = TRUE)
  }
  
}

####table flowering dates#####
param_consist_spp$day.flw.1<-(param_consist_spp$q1*(180/pi))*365/360
param_consist_spp$day.flw.2<-(param_consist_spp$q2*(180/pi))*365/360

param_consist_spp$flw.date.1<-as.Date(param_consist_spp$day.flw.1,origin = "2023-01-01")
param_consist_spp$flw.date.2<-as.Date(param_consist_spp$day.flw.2,origin = "2023-01-01")

param_variable_spp$day.flw.1<-(param_variable_spp$q1*(180/pi))*365/360
param_variable_spp$day.flw.2<-(param_variable_spp$q2*(180/pi))*365/360

param_variable_spp$flw.date.1<-as.Date(param_variable_spp$day.flw.1,origin = "2023-01-01")
param_variable_spp$flw.date.2<-as.Date(param_variable_spp$day.flw.2,origin = "2023-01-01")

##write tables parameters for selected models
#write.csv(param_consist_spp,"param_consist_spp.csv",row.names = F)
#write.csv(param_variable_spp,"param_variable_spp.csv",row.names = F)


###################### models with climate

### read data flowering with climate#$####

all.spp.clim<- read.csv("S2_2-data herbairum and climate 1900.csv",header=T)

colnames(all.spp.clim)

all.spp.clim<- all.spp.clim %>%
  filter(filter== "PASS")

var.include <-c("gbifID","hcnqID" , "eventDate" , "year",                    
                "month" , "day" , "decimalLongitude", "decimalLatitude" ,
                "family" , "genus" , "species" , "acceptedScientificName",
                "clusters","filter" , "n_records_flw_by_cluster" ,   
                "prec" , "prec_1" , "prec_2" , "prec_3" ,
                "tmin" , "tmin_1" , "tmin_2" , "tmin_3" ,
                "tmax" , "tmax_1",  "tmax_2" , "tmax_3" ,  
                "trange" , "trange_1" , "trange_2" , "trange_3" ,
                "new_day" , "new_month" , "new_date",                
                "photoperiod", "photoperiod_1", 
                "photoperiod_2" , "photoperiod_3")


col.index<-which(colnames(all.spp.clim) %in% var.include)

all.spp.climatology<-all.spp.clim[,col.index]

linked.data.flw[linked.data.flw$gbifID==2270689987, ]
all.spp.climatology[all.spp.climatology$gbifID==2270689987, ]  

summary(linked.data.flw$gbifID)
summary(all.spp.climatology$gbifID)

linked.data.clim<-merge(linked.data.flw,all.spp.climatology, by=c("clusters", "acceptedScientificName",                
                                                                  "gbifID","hcnqID", "eventDate" ,"year" ,                 
                                                                  "month" , "day","decimalLongitude"  ,    
                                                                  "decimalLatitude","family", "genus" ,                
                                                                  "species" ))
####Best models by region #####
datasets<-ls(pattern = "spp.plot")

best.model_byregion <- data.frame(NULL)

for (i in 1:6){
  for (j in 1:10){
  print(paste(i,j))
      spp.plot.i<- linked.data.flw %>%
      filter(clusters==i & Flowering_pred_RF=="Yes", total_count>=20)
      dist.mod<- circ_mle(
      spp.plot.i$fld.date.rad,
      criterion = "AIC",
      nchains = 5,
      BadStart = 10^9,
      niter = 5000,
      method = "BFGS",
      lambda.min = 0.25,
      exclude = NULL)

    best.model.i <- dist.mod$results[which(rownames(dist.mod$results) == dist.mod$bestmodel), ]
    stat.i<- dist.mod$rt[[1]]
    prob.i<- dist.mod$rt[[2]]
    best.model.i<-cbind(paste("region",i,sep=""),dist.mod$bestmodel,paste("run_",j,sep=""),best.model.i,stat.i,prob.i )
    colnames(best.model.i)[1]<-"region"
    colnames(best.model.i)[2]<-"best.model"
    colnames(best.model.i)[3]<-"run"
    best.model_byregion <-rbind(best.model_byregion,best.model.i)

  }

}

write.csv(best.model_byregion,"modelos_reg.csv",row.names = F)


summ.mod.region<-best.model_byregion%>%
  group_by(region, best.model) %>% 
  summarise(total_count=n(),minAIC=min(AIC),maxAIC=max(AIC),meanAIC=mean(AIC))

#write.csv(summ.mod.region,"modelos_reg_consist.csv",row.names = F)

summ.mod.region<-summ.mod.region[summ.mod.region$total_count>20,] ####esto debe ser 20 cuando se dan 100 corridas

best.model_byregion$best.model.reg<-paste(best.model_byregion$region,best.model_byregion$best.model,sep="")
summ.mod.region$best.model.reg<-paste(summ.mod.region$region,summ.mod.region$best.model,sep="")

regions<-paste0("region",unique(model.consistency$clusters))

table.param.mod.reg<-data.frame(NULL)
for (i in 1:length(regions)){ 
  param.data <- best.model_byregion%>%
    filter(region==regions[i]) 
  
  summ.data <- summ.mod.region%>%
    filter(region==regions[i]) 
  
    for (j in 1:nrow(summ.data)){ 
    
      param.data.j <- param.data%>%
        filter(best.model.reg==summ.data$best.model.reg[j]) 
      
      mod.consis<-summ.data$total_count[j]
      
      param.data.jj<-filter(param.data.j,AIC==min(param.data$AIC))
      if (nrow(param.data.jj) != 1) {
        param.data.jj<-param.data.j[1,]
      }else{
        param.data.jj=param.data.jj
      }
      param.data.jj<-cbind(param.data.jj,mod.consis)
      table.param.mod.reg<-rbind(table.param.mod.reg,param.data.jj)
  }
}

#write.csv(table.param.mod.reg, "table_modcoeff_reg.csv",row.names = F)
###### regional figures with climate ###################### 

dry.month<-read.csv("dry_months_region.csv", header=T)
month_day<-read.csv("month_day.csv",header=T)

regions<-unique(summ.mod.region$region)
color.reg<- c("#fac228","#f57d15","#d44842","#9f2a63","#65156e","#280b53")

for (i in 1:length(regions)){ 
  
  spp.plot<- linked.data.flw %>%
    filter(clusters==i & Flowering_pred_RF=="Yes", total_count>=20)
  
  param.data <- best.model_byregion%>%
  filter(region==regions[i]) 
  
  summ.data <- summ.mod.region%>%
    filter(region==regions[i]) 
  
  month_prec<- dry.month %>%
    filter(region==i)
  
  color.reg.i<-color.reg[i]
  
 for (j in 1:nrow(summ.data)){ 
    
    param.data.j <-  param.data%>%
    filter(best.model.reg==summ.data$best.model.reg[j]) 

    param.data.j<-filter(param.data,AIC==min(param.data$AIC))
      if (nrow(param.data.j) != 1) {
      param.data.j<-param.data.j[1,]
      }else{
      param.data.j=param.data.j
      }
    row.names(param.data.j)<-param.data.j$best.model
    param.data.j<-param.data.j[,-c(1:3)]

      if ( nrow(month_prec)==1){
      windows(5,5)
      plot_circMLE_JO(spp.plot$fld.date.rad, param.data.j, bins=12,
                    col.points=color.reg.i,sep=0.08, cex=0.8, shrink=2,
                    main.t=summ.data$best.model.reg[j])
    
      }else{
      days_x<-NULL
        for (i in 1:length(month_prec$month)){
        days.i<-seq(month_day[month_day$monyh==month_prec$month[i],3],
                  month_day[month_day$monyh==month_prec$month[i],4])
        days_x<-c(days_x,days.i)
        }
      windows(5,5)
      plot_circMLE_JO(spp.plot$fld.date.rad, param.data.j, bins=12,
                    col.points=color.reg.i,sep=0.08, cex=0.8, shrink=2,
                    main.t=summ.data$best.model.reg[j])
      lines.circular(x=circular(days_x*360/365*pi/180, zero = pi/2, rotation = "clock" ), 
                   y=rep(1,length(days_x)),join=F, 
                   col="red",lwd=3)
      }
  }
}  


######Analysis best model by latitudinal band within region#######


###intialize dataset to store restults
select.cluster<-c(1:6)

best.model_reg_lat <- data.frame(NULL)

spp.plot<- linked.data.flw %>%
    filter(clusters==6 & Flowering_pred_RF=="Yes", total_count>=20)
  
spp.plot.clim<- linked.data.clim %>%
    filter(clusters==select.cluster[6])
  
  spp.plot$lat.bin<-NA
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=15)]<-"from n15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=15)]<-"from s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<5&spp.plot$decimalLatitude>-5)]<-"from s5 to n5"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=-5&spp.plot$decimalLatitude>-15)]<-"from s5 to s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=5&spp.plot$decimalLatitude<15)]<-"from n5 to n15"
  
  spp.plot.clim$lat.bin<-NA
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude>=15)]<-"from n15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<=15)]<-"from s15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<5&spp.plot.clim$decimalLatitude>-5)]<-"from s5 to n5"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<=-5&spp.plot.clim$decimalLatitude>-15)]<-"from s5 to s15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude>=5&spp.plot.clim$decimalLatitude<15)]<-"from n5 to n15"
  
  spp.plot.clim<-merge(spp.plot.clim,main.dist,by=c("acceptedScientificName","clusters"))
  
lat.group<-unique(spp.plot$lat.bin)

spp_reg_lat_6<- spp.plot %>%
  group_by(lat.bin,acceptedScientificName) %>% 
  summarise(total_count=n())

#write.csv(spp_reg_lat_6,"species_per_reg_lat6.csv",row.names = F)

for (i in 1:5){
    spp.plot_lat<- spp.plot %>%
    filter(lat.bin == lat.group[i])

    for (j in 1:1){
    print(paste(i,j))
    dist.mod<- circ_mle(
      spp.plot_lat$fld.date.rad,
      criterion = "AIC",
      nchains = 5,
      BadStart = 10^9,
      niter = 5000,
      method = "BFGS",
      lambda.min = 0.25,
      exclude = NULL)

      best.model.i <- dist.mod$results[which(rownames(dist.mod$results) == dist.mod$bestmodel), ]
      stat.i<- dist.mod$rt[[1]]
      prob.i<- dist.mod$rt[[2]]
      best.model.i<-cbind(paste("region",6,sep=""),lat.group[i],dist.mod$bestmodel,paste("run_",j,sep=""),best.model.i,stat.i,prob.i  )
      colnames(best.model.i)[1]<-"region"
      colnames(best.model.i)[2]<-"subgroup"
      colnames(best.model.i)[3]<-"best.model"
      colnames(best.model.i)[4]<-"run"
      best.model_reg_lat <-rbind(best.model_reg_lat,best.model.i)
    }
}

####write restults
write.csv(best.model_reg_lat,"modelos_reg_lat.csv",row.names = F)

summ.mod.region.lat<-best.model_reg_lat%>%
  group_by(region,subgroup,best.model) %>% 
  summarise(total_count=n(),minAIC=min(AIC),maxAIC=max(AIC),meanAIC=mean(AIC))

#write.csv(summ.mod.region.lat,"modelos_reg_lat_consist.csv",row.names = F)

summ.mod.region.lat<-summ.mod.region.lat[summ.mod.region.lat$total_count>30,]

best.model_reg_lat$best.model.reg<-paste(best.model_reg_lat$region,
                                         best.model_reg_lat$subgroup,
                                         best.model_reg_lat$best.model,sep="")
summ.mod.region.lat$best.model.reg<-paste(summ.mod.region.lat$region,
                                          summ.mod.region.lat$subgroup,
                                          summ.mod.region.lat$best.model,sep="")


#table summary regions latitude #######

table.param.mod.reg.lat<-data.frame(NULL)
for (i in 1:length(regions)){ 

#for (i in 2:2){ 
  
  param.data <- best.model_reg_lat%>%
    filter(region==regions[i]) 
  
  summ.data <- summ.mod.region.lat%>%
    filter(region==regions[i]) 
  
for (j in 1:nrow(summ.data)){ 
#    for (j in 4:5){ 
    
    param.data.j <- param.data%>%
      filter(best.model.reg==summ.data$best.model.reg[j]) 
    
    mod.consis<-summ.data$total_count[j]
    
    param.data.jj<-filter(param.data.j,AIC==min(param.data$AIC))
    if (nrow(param.data.jj) != 1) {
      param.data.jj<-param.data.j[1,]
    }else{
      param.data.jj=param.data.jj
    }
    param.data.jj<-cbind(param.data.jj,mod.consis)
    table.param.mod.reg.lat<-rbind(table.param.mod.reg.lat,param.data.jj)
  }
}  

#write.csv(table.param.mod.reg.lat,"table_modcoeff_reg_lat.csv",row.names = F)

####estimate n per latitude subgroup ########
n_lat_reg_group<-spp.plot<-data.frame(NULL)
for (i in 1:length(regions)){ 
  spp.plot<- linked.data.flw %>%
    filter(clusters==i & Flowering_pred_RF=="Yes", total_count>=20)
  
  spp.plot$lat.bin<-NA
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=15)]<-"from n15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=15)]<-"from s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<5&spp.plot$decimalLatitude>-5)]<-"from s5 to n5"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=-5&spp.plot$decimalLatitude>-15)]<-"from s5 to s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=5&spp.plot$decimalLatitude<15)]<-"from n5 to n15"
  
  n_subgroup<-spp.plot%>%
    group_by(clusters,lat.bin) %>% 
    summarise(total_count=n())
  n_lat_reg_group<-rbind(n_lat_reg_group,n_subgroup)
}
#### figuras reg lat #############

regions<-unique(summ.mod.region.lat$region)
color.reg<- c("#fac228","#f57d15","#d44842","#9f2a63","#65156e","#280b53")

month_day<-read.csv("month_day.csv",header=T)


#for (i in 1:length(regions)){ 

for (i in 1:1){ 
  spp.plot<- linked.data.flw %>%
    filter(clusters==i & Flowering_pred_RF=="Yes", total_count>=20)
  
  spp.plot.clim<- linked.data.clim %>%
    filter(clusters==select.cluster[i])
  
  spp.plot$lat.bin<-NA
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=15)]<-"from n15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=15)]<-"from s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<5&spp.plot$decimalLatitude>-5)]<-"from s5 to n5"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude<=-5&spp.plot$decimalLatitude>-15)]<-"from s5 to s15"
  spp.plot$lat.bin[which(spp.plot$decimalLatitude>=5&spp.plot$decimalLatitude<15)]<-"from n5 to n15"
  
  spp.plot.clim$lat.bin<-NA
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude>=15)]<-"from n15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<=15)]<-"from s15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<5&spp.plot.clim$decimalLatitude>-5)]<-"from s5 to n5"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude<=-5&spp.plot.clim$decimalLatitude>-15)]<-"from s5 to s15"
  spp.plot.clim$lat.bin[which(spp.plot.clim$decimalLatitude>=5&spp.plot.clim$decimalLatitude<15)]<-"from n5 to n15"
  
  #spp.plot.clim<-merge(spp.plot.clim,main.dist,by=c("acceptedScientificName","clusters"))
  
  param.data <- best.model_reg_lat%>%
    filter(region==regions[i]) 
  
  summ.data <- summ.mod.region.lat%>%
    filter(region==regions[i]) 
  
  color.reg.i<-color.reg[i]
  
    for (j in 1:nrow(summ.data)){ 
    #for (j in 1:1){ 
    
    spp.plot.lat<- spp.plot %>%
      filter(lat.bin==summ.data$subgroup[j])
    
    spp.plot.clim.lat<- spp.plot.clim %>%
      filter(lat.bin==summ.data$subgroup[j])
    
    param.data.j <- param.data%>%
      filter(best.model.reg==summ.data$best.model.reg[j]) 
    
    param.data.jj<-filter(param.data.j,AIC==min(param.data$AIC))
    if (nrow(param.data.jj) != 1) {
      param.data.jj<-param.data.j[1,]
    }else{
      param.data.jj=param.data.jj
    }
    row.names(param.data.jj)<-param.data.jj$best.model
    param.data.jj<-param.data.jj[,-c(1:4)]
    
    windows(5,5)
    plot_circMLE_JO(spp.plot.lat$fld.date.rad, param.data.jj, bins=12,
                    col.points=color.reg.i,sep=0.08, cex=0.8, shrink=1.8,lwd=c(2,2,2),
                    main.t=summ.data$best.model.reg[j])
    
    month.photo<-spp.plot.clim.lat%>%
      group_by(month) %>% 
      summarise(month_photo=mean(photoperiod, na.rm=T))
    month.photo$month_photo<-month.photo$month_photo-12
    month.photo<- month.photo%>%
      filter(month_photo>0)
    
    days_x<-NULL
      for (i in 1:length(month.photo$month)){
        days.i<-seq(month_day[month_day$monyh==month.photo$month[i],3],
                  month_day[month_day$monyh==month.photo$month[i],4])
        days_x<-c(days_x,days.i)
      }
      sun_y<-NULL
      for (i in 1:length(month.photo$month)){
      sun.i<-rep(month.photo$month_photo[i] ,
                 month_day[month_day$monyh==month.photo$month[i],2])
            
      sun_y<-c(sun_y,sun.i)
      }
    
    lines.circular(x=circular(days_x*360/365*pi/180, zero = pi/2, rotation = "clock" ), 
                   y=sun_y+0.5,join=F,shrink=1.8,
                   col="grey50",lwd=2)
    
    month_prec <- spp.plot.clim.lat%>%
      group_by(month) %>% 
      summarise(month_prec=mean(prec))
    month_prec<- month_prec%>%
      filter(month_prec<100)
    
      days_x<-NULL
      for (i in 1:length(month_prec$month)){
      days.i<-seq(month_day[month_day$monyh==month_prec$month[i],3],
                  month_day[month_day$monyh==month_prec$month[i],4])
      days_x<-c(days_x,days.i)
      }
    
     lines.circular(x=circular(days_x*360/365*pi/180, zero = pi/2, rotation = "clock" ), 
                   y=rep(0.5,length(days_x)),join=F, 
                   col="red",lwd=2)
    }
}







