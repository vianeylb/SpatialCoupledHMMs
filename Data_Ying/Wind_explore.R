#############################################
#This code provides:
# 1) Basic Exploration: Map of the 20 sites, Monthly/Annual Windroses, and Time series plots
#############################################


#############################################
##Preliminaries:  Loading packages, Setting the working directory, 
## Reading in datasets
#############################################

rm(list=ls())

setwd("/Users/Karen/Documents/2016-Summer/KAUST/KAUST_CRG_DATASETS_v1/3. Northwest Wind")

library("fields")
library("chron")
library("plyr")
library("xtable")
library("leaps")
library("xts")
library("ggplot2")
library("ggmap")


#Hourly averaged data at 20 sites:
#-------------
load(file="Wind_data_2012to2014",verbose=T)
str(D);length(D);names(D)
ldply(D,function(x)x[1,]) #first rows of each list entry
ldply(D,function(x)x[dim(x)[1],]) #last rows of each list entry
#-------------


#Site Info:
#-------------
info=read.csv(file="BPA_Anemometer_Network_v2013.csv")
info
info$Datasets=names(D.60.12) #Corresponding names of datasets
names(info)[c(4,7,8,9)]=c("elev.ft","Lv1.ft","Lv2.ft","Lv3.ft")
#All datasets relate to Level 1, so just remove Levels 2 and 3
info=info[,-which(colnames(info)%in%c("Lv2.ft","Lv3.ft"))]
#Convert ft to meters
info$elev.m=info$elev.ft*.3048 
info$Lv1.m=info$Lv1.ft*.3048 
info=info[,c(1,8,2:7,9,10)]
info
#-------------


#############################################
##Basic Exploration: Terrain Map, Windroses
## Time series Plots
#############################################



#Terrain Map
#-------------
#Jitter Acronym Labels for ggplot2
#--------
relable_above<-c(1,2,6,8,13)
relable_below<-c(3,4,7,10,14,19,20,11)
relable<-c(relable_above,relable_below) #,relable_left
lable<-setdiff(seq(1,20,1),relable)

acr.loc<-cbind(info$Longitude,info$Latitude)
colnames(acr.loc)<-c("acr.x","acr.y")
acr.loc[lable,]=acr.loc[lable,]+cbind(rep(.2,length(lable)),rep(.1,length(lable)))
acr.loc[relable_above,]=acr.loc[relable_above,]+cbind(rep(-.01,length(relable_above)),rep(.14,length(relable_above)))
acr.loc[relable_below,]=acr.loc[relable_below,]+cbind(rep(-.01,length(relable_below)),rep(-.14,length(relable_below)))
loc=cbind(info,acr.loc)
#--------

lat.long <- c(lon=info$Longitude[7], lat=info$Latitude[7])
terrain.map <-get_googlemap(center = lat.long, zoom=7,maptype=c("terrain"))
dev.new(width=6,height=4)
ggmap(terrain.map,extent='device') + xlim(-125,-117) + ylim(44,47) + geom_point(data=loc, aes(x = Longitude, y = Latitude)) + geom_text(data=loc, aes(x = acr.x, y = acr.y, label=loc$Acronym), size=4, colour="black")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())
#-------------


# 10 locations w/ fully imputed datasets 
#-------------
NA.sets=ldply(D,function(x)any(is.na(x)))
comp.sites=NA.sets[which(!NA.sets$V1),]$.id #sites with complete datasets
as.matrix(comp.sites)
 # [1,] "BiddleButte" 
 # [2,] "ForestGrove" 
 # [3,] "HoodRiver"   
 # [4,] "HorseHeaven" 
 # [5,] "Megler"      
 # [6,] "NaselleRidge"
 # [7,] "Roosevelt"   
 # [8,] "Shaniko"     
 # [9,] "Sunnyside"   
# [10,] "Tillamook"  
#-------------

#Windroses:
#-------------
#Windrose Function:
#------
WR<-function(dataset,plot.title,new.win=F,wid=4,hgt=3){
#dataset - dataset to plot the wind rose for
	
	#So colors correspond to the same range of speeds across all plots:
	myColors=rainbow(6)
	# hist(rnorm(700),col=myColors)
	names(myColors)=c("0-5","5-10","10-15","15-20","20-25",">25")
	speed.scale=function(x){
		if(is.na(x)){
			scl="NA"
		}else{	
			if(0<=x && x<5) scl="0-5"
			if(5<=x&& x<10) scl="5-10"
			if(10<=x&& x<15) scl="10-15"
			if(15<=x&& x<20) scl="15-20"
			if(20<=x&& x<25) scl="20-25"
			if(25<= x) scl=">25"
		}
		return(scl)
	}
	fillScale<-scale_fill_manual(name="SPD",values = myColors)
	
	WR.data=dataset
	SPD=adply(WR.data$Speed,1,speed.scale)[,2]
	SPD=factor(SPD,levels=c("NA","0-5","5-10","10-15","15-20","20-25",">25"),ordered=T)
	WR.data=as.data.frame(cbind(WR.data,SPD))
	colnames(WR.data)[dim(WR.data)[2]]="SPD"
	head(WR.data)
	if(new.win)dev.new(width=wid,height=hgt)
	print(ggplot(WR.data, aes(x=Dir,fill=SPD))+fillScale+geom_histogram(aes(Speed=..density..,fill=SPD),binwidth=10)+coord_polar(theta="x",start=0,direction=1)+theme(panel.grid.major = element_line(colour = "black", size = 0.2, linetype=2))+labs(x="",y="") + scale_x_continuous(limits=c(0,360),breaks=seq(0,330,30))+theme(legend.position = "left")+theme(plot.title=element_text(size=15))+ggtitle(paste(plot.title,sep="")))	
}
#------

#Plot Windroses:
#------
# For each month of each year and for each full year
# as.matrix(names(D))

i=2;names(D)[i]  #Select a site
# for(i in 1:length(D)){ #Loop might take a while	
	for(yr in c(2012,2013,2014)){ #yr=2012
		for(mnth in 1:12){ #mnth=1
			mnth.yr.set=which(D[[i]]$month==mnth & D[[i]]$year==yr)					
			WR(D[[i]][mnth.yr.set,],paste(names(D)[i]," ",levels(months(1))[mnth]," ",yr,sep=""),new.win=T)	
		}
		yr.set=which(D[[i]]$year==yr)
		WR(D[[i]][yr.set,],paste(names(D)[i]," ",yr,sep=""),new.win=T)	
	}		
# }
#------
#-------------


#Time series Plots
#-------------
i=1
for(i in 1:length(D)){
	dev.new(width=9,height=6)
	par(mfrow=c(4,1),mai=c(.4,.7,.4,.4))
	
	#Dates on axes:
	mnth.change=c(1,which(!diff(D[[i]]$month)==0)+1)
	labels=apply(D[[i]][mnth.change,],1,function(x)paste(x["month"],"/",substring(x["year"],3,4),sep=""))
	
	plot(D[[i]]$Speed,type="l",ylab="speed",xaxt="n")	
	axis(1,at=mnth.change,labels=labels)
	plot(D[[i]]$Dir,type="l",ylab="direction",xaxt="n")
	axis(1,at=mnth.change,labels=labels)
	plot(D[[i]]$U,type="l",ylab="U",xaxt="n")
	axis(1,at=mnth.change,labels=labels)
	plot(D[[i]]$V,type="l",ylab="V",xaxt="n")
	axis(1,at=mnth.change,labels=labels)
	title(names(D)[i],line=-2,cex.main=2.5,outer=T)
}
#-------------






















