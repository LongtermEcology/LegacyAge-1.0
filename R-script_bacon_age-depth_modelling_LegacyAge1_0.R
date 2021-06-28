
#########################################################################################
###
### Supplement to the "Harmonized chronologies of a global late Quaternary pollen dataset (LegacyAge 1.0)"
###
### R-Script for establishing bacon age-depth-models
###
### by Alexander Karl Postl
###
### Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Germany 2021
###
#########################################################################################


### -----General Beginning-----

{
  
  # clean memory and set the working directory:
  rm(list=ls())
  folder<-"set your working directory here"
  setwd(folder)
  
  # libraries:
  library(tidyverse)
  library(neotoma)
  library(stringr)
  library(ggpubr)
  library(rlist)
  library(rbacon)
  library(clam)
  library(IntCal)
  
  # additional functions:
  `%notin%` <- Negate(`%in%`)
  numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
  }
  
  
  # loading files:
  AgeDepthPollen<-read.csv2(paste0(folder,"/AgeDepthPollen_dataset.csv"),stringsAsFactors = FALSE,sep = ";",dec=".")
  composed.cores<-read.csv2(paste0(folder,"/composed_cores.csv"),stringsAsFactors = FALSE,sep = ";",dec=".")
  cao<-read.csv2(paste0(folder,"/Sites_by_Cao.csv"),stringsAsFactors = FALSE,sep = ";",dec=".")
  parameter<-read.csv2(paste0(folder,"/Bacon_parameters.csv"),header = TRUE,stringsAsFactors = FALSE, dec = ".",sep = ";")
  
  # define values:
  delete.chrono<-c("Radiocarbon years BP","No chronology",NA)
  c14<-c("Carbon-14","14C","C14","14c","c14","","Radiocarbon years BP")
  datings<-NULL
  
  # error reporting:
  bacon.error<-read.csv2(paste0(folder,"/bacon_error.csv"),header = TRUE, stringsAsFactors = FALSE,sep = ";",dec=".")[,1]
  r.abort<-read.csv2(paste0(folder,"/r_abort.csv"),stringsAsFactors = FALSE,sep = ";",dec=".")[,1]
  
}

# ---------------------------------------------------------------------------------------------

### -----Loop for every ID-----

# defining ID list:
parameter.ID<-parameter[which(parameter$rerun.==TRUE),]$ID

# loop run for every ID:
for (ID in parameter.ID){
  
  # error protocol if R collapses:
  r.abort<-unique(c(r.abort,ID))
  write.csv2(r.abort,"r_abort.csv",row.names = FALSE)
  r.abort<-r.abort[-which(r.abort==ID)]
  
  # parameter settings for the specific ID:
  parameter.subset<-parameter[which(parameter$ID==ID),]
  if (is.na(parameter.subset$reservoir)){parameter.subset$reservoir<-0}
  if (is.na(parameter.subset$waterline)){parameter.subset$waterline<-0}  
  
  # reset of values and defining the download ID's for collecting data in the case of composed cores:
  calibration=list_ages=AgeDepthPollen.subset=cal.age=cal.age.se=cal.sp<-NULL
  if (ID<70000|ID>=100000){neotoma.download<-ID}else{neotoma.download<-composed.cores[composed.cores$ID.complete==ID,2]}
  
  # ---------------------------------------------------------------------------------------------
  
  ### -----Data download & collection-----  
  
  {
    
    # loop for downloading and collecting data:
    for(dataset.id in neotoma.download){
      
      # reset values:
      AgeDepthPollen.subset.new=core=cali=calibration.new<-NULL
      
      for (i in 1:6){
        assign(paste0("chronology.name_",i),NA)
        assign(paste0("chronology.age_",i),NA)
        assign(paste0("chronology.error_",i),NA)
        assign(paste0("chronology.type_",i),NA)
      }    
      
      # ---------------------------------------------------------------------------------------------
      
      # collecting pollen-data from Neotoma or from Cao et al.:
      
      {
        
        if(ID<100000){  
          
          # repeat loop to download the pollen data from Neotoma:
          repeat{
            try(core<-get_download(dataset.id),silent = TRUE)
            if(is.list(core)){break}
          }
          
          # loop for adding Neotoma chronologies:
          for (i in 1:min(6,length(core[[1]]$chronologies))){
            assign(paste0("chronology.name_",i),core[[1]]$chronologies[[i]]$chronology.name)
            assign(paste0("chronology.age_",i),core[[1]]$chronologies[[i]]$age)
            assign(paste0("chronology.error_",i),(core[[1]]$chronologies[[i]]$age.older-core[[1]]$chronologies[[i]]$age.younger)/2)
            assign(paste0("chronology.type_",i),core[[1]]$chronologies[[i]]$age.type)
            
            # detecting the chronologies with calibrated ages:
            if(any(core[[1]]$chronologies[[i]]$age.type %notin% delete.chrono)){list_ages<-c(list_ages,core[[1]]$chronologies[[i]]$age)}      
          }
          
          # mean age of Neotoma calibrated chronologies as reference:
          if (length(list_ages)==0){mean.age<-NA}else{mean.age<-mean(list_ages)}
          
          # correction if depth is given in meters:
          if (ID%in%c()){core[[1]]$sample.meta$depth<-core[[1]]$sample.meta$depth*100}
          
          # pollen depths of composed cores:
          if (ID>=70000){core[[1]]$sample.meta$depth<-core[[1]]$sample.meta$depth+composed.cores[which(composed.cores$dataset.id==dataset.id),3]}
          
          # use data from the AgeDepthPollen table if there are entries of the same length (just in case there have been manual modifications based on publications):
          if(length(core[[1]]$sample.meta$depth)==length(AgeDepthPollen[which(AgeDepthPollen$dataset.id==ID),16])){core[[1]]$sample.meta$depth<-AgeDepthPollen[which(AgeDepthPollen$dataset.id==ID),16]}
          
          # creating the AgeDepthPollen table for the specific ID:
          AgeDepthPollen.subset.new<-data.frame(core[[1]]$dataset$site.data$site.id,core[[1]]$dataset$site.data$site.name,core[[1]]$dataset$site.data$long,core[[1]]$dataset$site.data$lat,core[[1]]$dataset$site.data$elev,parameter.subset$Region,NA,NA,NA,NA,core[[1]]$dataset$dataset.meta$dataset.id,core[[1]]$dataset$dataset.meta$dataset.name,core[[1]]$dataset$dataset.meta$collection.type,core[[1]]$dataset$dataset.meta$collection.handle,core[[1]]$dataset$dataset.meta$dataset.type,core[[1]]$sample.meta$depth,mean.age,NA,NA,NA,NA,NA,NA,NA,NA,chronology.name_1,chronology.type_1,chronology.age_1,chronology.error_1,chronology.name_2,chronology.type_2,chronology.age_2,chronology.error_2,chronology.name_3,chronology.type_3,chronology.age_3,chronology.error_3,chronology.name_4,chronology.type_4,chronology.age_4,chronology.error_4,chronology.name_5,chronology.type_5,chronology.age_5,chronology.error_5,chronology.name_6,chronology.type_6,chronology.age_6,chronology.error_6,core[[1]]$dataset$site.data$description)
          colnames(AgeDepthPollen.subset.new)<-names(AgeDepthPollen)
          AgeDepthPollen.subset<-rbind(AgeDepthPollen.subset,AgeDepthPollen.subset.new) 
          
        }else{
          
          # read pollen data from Cao et al.:
          cao.id<-cao[which(cao$Asia.Site.ID==(ID-100000)),][1,]
          file<-read.csv2(paste0("N:/bioing/data/NorthernHemisphericPollendataset/AgeDepth_Modelling/Asia/sites/",ID-100000,".csv"),header = FALSE, stringsAsFactors = FALSE)
          pollen<-data.frame(t(file[1,-1]),t(file[2,-1]),stringsAsFactors = FALSE)      
          
          # creating the AgeDepthPollen table for specific ID:
          AgeDepthPollen.subset<-data.frame(paste0("Asia_",ID-100000),cao.id$SiteName,as.numeric(cao.id$Longitude),as.numeric(cao.id$Latitude),as.numeric(cao.id$Elevation),"Asia",NA,NA,NA,NA,ID,NA,NA,NA,NA,as.numeric(pollen[,1]),as.numeric(pollen[,2]),stringsAsFactors = FALSE)
          AgeDepthPollen.subset[,(18:50)]<-NA
          colnames(AgeDepthPollen.subset)<-names(AgeDepthPollen)
          
        }
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
      # collecting dating-data from Neotoma or from Cao et al.:
      
      {
        
        if (ID<100000){
          
          # repeat loop for downloading the datasets of the site in Neotoma:
          repeat{
            try(cali<-get_dataset(get_site(as.integer(core[[1]]$dataset$site.data$site.id))),silent = TRUE)
            if(is.list(cali)){break}
          }
          
          # repeat loop to detect the corresponding geochronology ID from Neotoma:
          element_c14<-if (ID==4311){2}else{0}
          repeat {
            element_c14<-element_c14+1
            if(cali[[element_c14]]$dataset.meta$dataset.type=="geochronologic" & core[[1]]$dataset$dataset.meta$collection.handle==cali[[element_c14]]$dataset.meta$collection.handle){break}
          }
          
          # repeat loop for downloading the corresponding geochronology data from Neotoma:
          repeat{
            try(cali<-get_geochron(cali[[element_c14]]$dataset.meta$dataset.id)[[1]]$geochron,silent = TRUE)
            if(is.list(cali)){break}
          }
          
          # correction if depth is given in meters:
          if (ID%in%c(47608)){cali$depth<-cali$depth*100}
          
          # dating depths of composed cores:
          if (ID>=70000){cali$depth<-cali$depth+composed.cores[which(composed.cores$dataset.id==dataset.id),3]}
          
          # AD to BP correction:
          if (any(cali[which(cali$age.type=="Calendar years AD/BC"),5]>=1950)){cali[which(cali$age.type=="Calendar years AD/BC"),5]<-1950-cali[which(cali$age.type=="Calendar years AD/BC"),5]}
          
          # creating the dating table for calibration (Neotoma):
          calibration.new<-data.frame(depth=cali$depth,age=cali$age,e.older=cali$e.older,age.type=cali$age.type,cal.age=cali$depth,cal.age.se=cali$depth,stringsAsFactors = FALSE)
          calibration<-rbind(calibration,calibration.new)
          cali<-cbind(datasetID=ID,cali)
          datings.new<-rbind(datings,cali)
          
        }else{
          
          # creating the dating table for calibration (Cao et al.):
          calibration<-cao[which(cao$Asia.Site.ID+100000==ID),]
          calibration<-data.frame(depth=calibration$Depth..cm.,age=calibration$Age.BP,e.older=calibration$Error,age.type=calibration$age.type,cal.age=calibration$Depth..cm.,cal.age.se=calibration$Depth..cm.,stringsAsFactors = FALSE)
          cali<-cao[which(cao$Asia.Site.ID+100000==ID),]
          datings.new<-data.frame(datasetID=ID,sample.id=cali$Asia.Site.ID,depth=cali$Depth..cm.,thickness=cali$Thickness,age.type=NA,age=cali$Age.BP,e.older=cali$Error,e.young=cali$Error,delta13C=NA, lab.no=NA, material.dated=cali$Description, geo.chron.type=cali$age.type,notes=NA, infinite=FALSE,stringsAsFactors = FALSE)
        }
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
    } 
    
    # ---------------------------------------------------------------------------------------------
    
    # eliminate datings defined by the parameters.csv:
    delete.data<-as.numeric(parameter.subset[,(which(names(parameter)=="delete.1"):dim(parameter)[2])])
    keep.data<-as.numeric(parameter.subset[,(which(names(parameter)=="keep.age.1"):(which(names(parameter)=="delete.1")-1))])
    calibration<-calibration[which(calibration$depth %notin% delete.data | calibration$age %in% keep.data),]
    datings.new<-datings.new[which(datings.new$depth %notin% delete.data | datings.new$age %in% keep.data),]
    
    # add datings defined by the parameters.csv:
    add.data<-data.frame(NA,NA,NA,NA,NA,NA)
    names(add.data)<-names(calibration)
    
    for (a in (1:((which(names(parameter)=="keep.age.1")-which(names(parameter)=="add.depth.1"))/4))) {
      add.new<-data.frame(parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4],parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4+1],parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4+2],parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4+3],parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4],parameter.subset[,which(names(parameter)=="add.depth.1")+(a-1)*4])
      names(add.new)<-names(add.data)
      add.new[,c(1:3,5:6)]<-as.numeric(add.new[,c(1:3,5:6)])
      add.data<-rbind(add.data,add.new)
    }
    
    add.data<-add.data[-which(is.na(add.data[,1])),]
    add.data[which(add.data[,4]==TRUE),4]<-"Radiocarbon years BP"
    add.data[which(add.data[,4]==FALSE),4]<-"Calendar years BP"
    calibration<-rbind(calibration,add.data)
    
    # prepare datings csv:
    add.datings<-data.frame(ID,NA,add.data[,1],NA,add.data[,4],add.data[,2],add.data[,3],add.data[,3],NA,NA,NA,"added dating",NA,NA)
    names(add.datings)<-names(datings.new)
    data<-rbind(datings,datings.new) 
    write.table(datings,"datings.csv", row.names = FALSE,col.names = TRUE, sep=";",dec = ".")
    
  }
  
  # ---------------------------------------------------------------------------------------------
  
  ### -----Calibration----- 
  
  {
    
    # defining the calibration curve and their borders:
    if(AgeDepthPollen.subset$lat[1]>=0){
      cc<-1
      cc1<-min(IntCal::copyCalibrationCurve(1)[,2]-IntCal::copyCalibrationCurve(1)[,3])
      cc2<-max(IntCal::copyCalibrationCurve(1)[,2]+IntCal::copyCalibrationCurve(1)[,3])}
    
    if(AgeDepthPollen.subset$lat[1]<0){
      cc<-3
      cc1<-min(IntCal::copyCalibrationCurve(3)[,2]-IntCal::copyCalibrationCurve(3)[,3])
      cc2<-max(IntCal::copyCalibrationCurve(3)[,2]+IntCal::copyCalibrationCurve(3)[,3])}
    
    if(parameter.subset$marine){
      cc<-2
      cc1<-min(IntCal::copyCalibrationCurve(2)[,2]-IntCal::copyCalibrationCurve(2)[,3])
      cc2<-max(IntCal::copyCalibrationCurve(2)[,2]+IntCal::copyCalibrationCurve(2)[,3])}
    
    
    # eliminating NA and separating uncalibrated C14 from the other datings:
    if (is.na(parameter.subset$waterline)){parameter.subset$waterline<-0}
    
    for (j in (1:nrow(calibration))){
      if (is.na(calibration[j,3]))   {if (calibration[j,2]<11460){calibration[j,3]<-30}else{calibration[j,3]<-60}}
      if (calibration[j,4]%in%c14 & (calibration[j,2]+calibration[j,3]-parameter.subset$reservoir)>cc1 & (calibration[j,2]-calibration[j,3]-parameter.subset$reservoir)<cc2) {calibration[j,4]<-cc}else{calibration[j,4]<-0}
    }
    
    calibration.0<-calibration[which(calibration$age.type==0),]
    calibration.0[,(5:6)]<-calibration.0[,(2:3)]
    calibration.1<-calibration[which(calibration$age.type==cc),]
    
    # defining the postbomb:
    if (AgeDepthPollen.subset$lat[1]>40) {postbomb<-1}
    if (AgeDepthPollen.subset$lat[1]<40) {postbomb<-2}
    if (AgeDepthPollen.subset$lat[1]<0) {postbomb<-4} 
    
    # are there non calibrated C14 datings? & calibration process:
    if (dim(calibration.1)[1]>=1){
      
      # calibration process for every non-calibrated C14 dating:
      for (clam in 1:dim(calibration.1)[1]){
        
        #calibration by clam:  
        calibration.calc<-calibrate(
          cage = calibration.1$age[clam],
          error = calibration.1$e.older[clam],
          cc = as.numeric(calibration.1$age.type[clam]),
          postbomb = postbomb,
          reservoir= parameter.subset$reservoir,
          graph = TRUE,
          title=paste0(ID,"_",calibration.1$depth[clam]),
          mar=c(3.5,3,2,1),
          mgp = c(1.7,0.8,0),
          bty = "n"
        )
        
        dev.print(pdf,paste0(folder,"/calibration/",ID,"_",calibration.1$depth[clam],".pdf") )
        
        # standardization of the deviation & weighted mean age:
        standart<-(1/sum(calibration.calc$calib[,2]))*calibration.calc$calib[,2]
        new.memory<-sum(standart*calibration.calc$calib[,1])
        new.se<-sum(abs(standart*(calibration.calc$calib[,1]-new.memory)))
        spread<-data.frame(calibration.1$depth[clam],calibration.calc$hpd[,1],calibration.calc$hpd[,2])
        
        # setting values (weighted mean age, weighted standard deviation):
        cal.age<-c(cal.age,new.memory)
        cal.age.se<-c(cal.age.se,new.se)
        cal.sp<-rbind(cal.sp,spread)
        
      }
      
      #finalize calibrated data for the calibration table:
      calibration.1$cal.age<-cal.age
      calibration.1$cal.age.se<-cal.age.se
      
    }
    
    # write calibration table:
    calibration<-rbind(calibration.0,calibration.1)[order(calibration[,1]),]
    write.table(cbind(calibration,reservoir=parameter.subset$reservoir),paste0(folder,"/calibration/",ID,".csv"),row.names = FALSE,col.names = TRUE,sep = ";",dec = ".")
    
  }
  
  # ---------------------------------------------------------------------------------------------
  
  ### -----Bacon modelling-----
  
  {
    
    # create csv for bacon:
    data<-data.frame(sample_id=AgeDepthPollen.subset$dataset.id[1],age=calibration$age ,error=calibration$e.older ,depth=calibration$depth,cc=calibration$age.type,stringsAsFactors = FALSE)
    zero<-data.frame(AgeDepthPollen.subset$dataset.id[1],max(-70,parameter.subset$starting.age,na.rm=TRUE),30,parameter.subset$waterline,0)
    colnames(zero)<-colnames(data)
    
    # is in the parameter csv a surface added or not:
    if (parameter.subset$add.surface.) {data<-rbind(zero,data)}
    
    # write csv for bacon:
    dir.create(paste0(folder,"/sites/",ID)) 
    file.remove(paste0(paste0(folder,"/sites/",ID,"/",list.files(path = paste0(folder,"/sites/",ID)))))
    write.table(data[order(data$depth),],file=paste0(folder,"/sites/",ID,"/",ID,".csv"),row.names=FALSE,sep = ",",dec = ".")  
    
    # calculate accumulation rate:
    if (is.na(parameter.subset$acc.rate)){
      acc.rate<-(max(data[,2])-min(data[,2]))/(max(data$depth)-min(data$depth))
    }else{acc.rate<-parameter.subset$acc.rate}
    
    # AgeDethModelling by Bacon:
    try(Bacon(as.character(ID)
              ,ask=FALSE
              ,suggest=FALSE
              ,thick = if(is.na(parameter.subset$resolution.cm)){(max(AgeDepthPollen.subset$depth,calibration$depth)-parameter.subset$waterline)/parameter.subset$resolution.sections}else{parameter.subset$resolution.cm}
              ,acc.shape = 1.5
              ,acc.mean = acc.rate
              ,mem.strength = 20
              ,mem.mean = parameter.subset$memory
              ,d.min = min(parameter.subset$waterline,AgeDepthPollen.subset$depth)
              ,d.max = max(AgeDepthPollen.subset$depth,calibration$depth)
              ,MinYr = parameter.subset$starting.age
              ,postbomb = postbomb
              ,hiatus.depths = unique(c(as.numeric(parameter.subset[1,which(names(parameter.subset)=="hiatus.1"):(which(names(parameter)=="add.depth.1")-1)])))
              ,coredir = paste0(folder,"/sites/")
              ,close.connections = TRUE
              ,delta.R = parameter.subset$reservoir
    ),silent = TRUE)
    
  }
  
  # ---------------------------------------------------------------------------------------------
  
  ### -----Bacon review-----  
  
  {
    
    # if bacon was not able to create a model:
    if (file.exists(paste0(folder,"/sites/",ID,"/",list.files(path = paste0(folder,"/sites/",ID),".pdf")))==FALSE) {
      
      # protocol the ID:
      bacon.error<-read.csv2(paste0(folder,"/bacon_error.csv"),stringsAsFactors = FALSE,sep = ";",dec=".")[,1]
      bacon.error<-unique(c(bacon.error,ID))
      write.csv2(bacon.error,"bacon_error.csv",row.names = FALSE)
      
      # delete the ID from the list for the outer loop and go to the next ID:
      parameter.ID<-parameter.ID[-1]
      write.csv2(parameter.ID,"bacon_ID.csv",row.names = FALSE)
      
    }else{
      
      # copy model files to new folders:
      file.copy(from = paste0(folder,"/sites/",ID,"/",list.files(path = paste0(folder,"/sites/",ID),".pdf")), to = paste0(folder,"/pdf/",ID,".pdf"),overwrite = TRUE)
      file.copy(from = paste0(folder,"/sites/",ID,"/",list.files(path = paste0(folder,"/sites/",ID),"ages.txt")), to = paste0(folder,"/txt/",ID,".txt"),overwrite = TRUE)    
      
      # ---------------------------------------------------------------------------------------------
      
      ### -----Age Allocation-----
      
      {
        
        # define basic values for the age allocation:
        model.AWI<-read.table(paste0(folder,"/txt/",ID,".txt"),header = TRUE)  
        model.rest <- round(model.AWI$depth[length(model.AWI$depth)]-floor(model.AWI$depth[length(model.AWI$depth)]),4)
        new.memory <- data.frame(min=NA, max=NA, median=NA, mean=NA)
        
        # loop for every pollen entry:
        for(g in 1:length(AgeDepthPollen.subset$depth)){
          
          # define basic values for the specific pollen entry: 
          pollen.num <- as.numeric(as.character(floor(AgeDepthPollen.subset$depth[g]-model.rest)+model.rest))  
          pollen.rest <- AgeDepthPollen.subset$depth[g] - pollen.num
          
          # Case I:
          if(AgeDepthPollen.subset$depth[g] < model.AWI$depth[1]){
            lower.age <- model.AWI[1, c(2:5)]
            upper.age <- model.AWI[2, c(2:5)]
            m <- (upper.age-lower.age)/(model.AWI$depth[2]-model.AWI$depth[1])
            b <- upper.age - m*model.AWI$depth[2]
            new.memory[g, ] <- m*AgeDepthPollen.subset$depth[g]+b
          }
          
          # Case II:
          if(AgeDepthPollen.subset$depth[g] >= model.AWI$depth[1] & AgeDepthPollen.subset$depth[g] < model.AWI$depth[length(model.AWI$depth)]){
            lower.age <- model.AWI[model.AWI$depth == pollen.num, c(2:5)]
            upper.age <- model.AWI[model.AWI$depth == pollen.num+1, c(2:5)]
            new.memory[g, ] <- lower.age*(1-pollen.rest) + upper.age*pollen.rest
          }
          
          # Case III:
          if(AgeDepthPollen.subset$depth[g] == model.AWI$depth[length(model.AWI$depth)]){new.memory[g, ] <- model.AWI[length(model.AWI$depth), c(2:5)]}
          
          #Case IV:
          if(AgeDepthPollen.subset$depth[g] > model.AWI$depth[length(model.AWI$depth)]){
            lower.age <- model.AWI[length(model.AWI$depth)-1, c(2:5)]
            upper.age <- model.AWI[length(model.AWI$depth), c(2:5)]
            m <- (upper.age-lower.age)/(model.AWI$depth[length(model.AWI$depth)]-model.AWI$depth[length(model.AWI$depth)-1])
            b <- upper.age - m*model.AWI$depth[length(model.AWI$depth)]
            new.memory[g, ] <- m*AgeDepthPollen.subset$depth[g]+b
          }
        }
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
      ### -----Preparing AgeDepthPollen table-----
      
      {
        
        # statistics:
        difference_age.perc<-100*((new.memory$mean+as.numeric(format(Sys.time(),"%Y"))-1950)/(AgeDepthPollen.subset$mean.age+as.numeric(format(Sys.time(),"%Y"))-1950))
        min.max<-data.frame(AgeDepthPollen.subset$mean.age>=new.memory$min & AgeDepthPollen.subset$mean.age<=new.memory$max,AgeDepthPollen.subset$mean.age>=(2*new.memory$min-new.memory$mean) & AgeDepthPollen.subset$mean.age<=(2*new.memory$max-new.memory$mean))
        
        # create & write table:
        memory.sample<-data.frame(AgeDepthPollen.subset$mean.age-new.memory$mean,difference_age.perc,min.max)
        memory.core<-data.frame(all(memory.sample[,3]),all(memory.sample[,4]),mean(memory.sample[,1]),sd(memory.sample[,1]))
        AgeDepthPollen.subset[,c((7:10),(18:25))]<-data.frame(memory.core,new.memory[,4],memory.sample,new.memory[,c(3,(1:2))])
        write.table(AgeDepthPollen.subset,paste0(folder,"/Subsets/",ID,".csv"), row.names = FALSE,col.names = TRUE, sep=";",dec = ".")
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
      #-----Plotting-----
      
      {
        
        # define the Neotoma-chronologies to plot (only calibrated chronologies):
        if (any(AgeDepthPollen.subset$chronology.type_1 %notin% delete.chrono)){chronology1<-data.frame(AgeDepthPollen.subset$chronology.name_1,AgeDepthPollen.subset$chronology.age_1,AgeDepthPollen.subset$chronology.error_1)}else{chronology1<-NULL}
        if (any(AgeDepthPollen.subset$chronology.type_2 %notin% delete.chrono)){chronology2<-data.frame(AgeDepthPollen.subset$chronology.name_2,AgeDepthPollen.subset$chronology.age_2,AgeDepthPollen.subset$chronology.error_2)}else{chronology2<-NULL}
        if (any(AgeDepthPollen.subset$chronology.type_3 %notin% delete.chrono)){chronology3<-data.frame(AgeDepthPollen.subset$chronology.name_3,AgeDepthPollen.subset$chronology.age_3,AgeDepthPollen.subset$chronology.error_3)}else{chronology3<-NULL}
        if (any(AgeDepthPollen.subset$chronology.type_4 %notin% delete.chrono)){chronology4<-data.frame(AgeDepthPollen.subset$chronology.name_4,AgeDepthPollen.subset$chronology.age_4,AgeDepthPollen.subset$chronology.error_4)}else{chronology4<-NULL}
        if (any(AgeDepthPollen.subset$chronology.type_5 %notin% delete.chrono)){chronology5<-data.frame(AgeDepthPollen.subset$chronology.name_5,AgeDepthPollen.subset$chronology.age_5,AgeDepthPollen.subset$chronology.error_5)}else{chronology5<-NULL}
        if (any(AgeDepthPollen.subset$chronology.type_6 %notin% delete.chrono)){chronology6<-data.frame(AgeDepthPollen.subset$chronology.name_6,AgeDepthPollen.subset$chronology.age_6,AgeDepthPollen.subset$chronology.error_6)}else{chronology6<-NULL}
        chronology<-list(chronology1,chronology2,chronology3,chronology4,chronology5,chronology6)
        chronology<-list.clean(chronology, function(x) length(x) == 0L, TRUE)
        
        # use the Cao et al.-chronologies for Cao-sites:
        if (ID>100000){chronology<-data.frame(AgeDepthPollen.subset$site.name,AgeDepthPollen.subset$mean.age,AgeDepthPollen.subset$mean.age)}
        if (length(chronology)!=0){chronology<-data.frame(AgeDepthPollen.subset$depth,chronology)}else{chronology<-data.frame(AgeDepthPollen.subset$depth)}
        
        # set strings for ggplot:
        plot.c14.correct<-"age.dating"
        plot.c14.added<-"age.dating.added"
        plot.AWI<-"AWI.age"
        
        #collect and separate added and original datings:
        zero<-data.frame(zero$depth,zero$age,zero$error,"",zero$age,zero$error)
        colnames(zero)<-names(calibration)[1:6]
        if (parameter.subset$add.surface.) {calibration<-rbind(calibration,zero)}
        calibration.added<-calibration[which(calibration$depth%in% c(add.data$depth,zero$depth)),]
        calibration<-calibration[which(calibration$depth%notin% c(add.data$depth,zero$depth)),]
        
        #values for ggplot without a chronology:
        if (dim(chronology)[2]>=1){
          y.plot.1<-chronology[,1]
          y.plot.2<-chronology[,1]
          y.plot.3<-chronology[,1]
          y.plot.4<-chronology[,1]
          y.plot.5<-chronology[,1]
          color.plot.1<-"meanAgeBP"
          color.plot.2<-"meanAgeBP"
          color.plot.3<-"meanAgeBP"
          color.plot.4<-"meanAgeBP"
          color.plot.5<-"meanAgeBP"
          alpha.plot.1<-0
          alpha.plot.2<-0
          alpha.plot.3<-0
          alpha.plot.4<-0
          alpha.plot.5<-0 
          alpha.added<-1
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI)
        }
        
        # with one chronology:
        if (dim(chronology)[2]>=4){
          y.plot.1<-chronology[,3]
          color.plot.1<-"Neotoma.1"
          alpha.plot.1<-0.8
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI,chronology[which(!is.na(chronology[,2])),2][1])
        }
        
        # with two chronologies:
        if (dim(chronology)[2]>=7){
          y.plot.2<-chronology[,6]
          color.plot.2<-"Neotoma.2"
          alpha.plot.2<-0.8
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI,chronology[which(!is.na(chronology[,2])),2][1],chronology[which(!is.na(chronology[,5])),5][1])    
        }
        
        # with three chronologies:
        if (dim(chronology)[2]>=10){
          y.plot.3<-chronology[,9]
          color.plot.3<-"Neotoma.3"
          alpha.plot.3<-0.8
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI,chronology[which(!is.na(chronology[,2])),2][1],chronology[which(!is.na(chronology[,5])),5][1],chronology[which(!is.na(chronology[,8])),8][1])     
        }
        
        # with four chronologies:
        if (dim(chronology)[2]>=13){
          y.plot.4<-chronology[,12]
          color.plot.4<-"Neotoma.4"
          alpha.plot.4<-0.8
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI,chronology[which(!is.na(chronology[,2])),2][1],chronology[which(!is.na(chronology[,5])),5][1],chronology[which(!is.na(chronology[,8])),8][1],chronology[which(!is.na(chronology[,11])),11][1]) 
        } 
        
        # with five chronologies:
        if (dim(chronology)[2]>=16){
          y.plot.5<-chronology[,12]
          color.plot.5<-"Neotoma.5"
          alpha.plot.5<-0.8
          naming<-list(plot.c14.correct,plot.c14.added,plot.AWI,chronology[which(!is.na(chronology[,2])),2][1],chronology[which(!is.na(chronology[,5])),5][1],chronology[which(!is.na(chronology[,8])),8][1],chronology[which(!is.na(chronology[,11])),11][1],chronology[which(!is.na(chronology[,14])),14][1]) 
        }
        
        # masking datings without data:
        if (dim(calibration.added)[1]==0) {
          calibration.added<-calibration
          alpha.added<-0}
        if (dim(calibration)[1]==0) {
          calibration<-calibration.added
          alpha.added<-0}
        
        # ---------------------------------------------------------------------------------------------
        
        # ggplot:
        {
          
          output<-ggplot() +
            
            # curves:
            geom_ribbon(data=model.AWI,aes(depth,ymin=min,ymax=max),alpha=0.4) +  
            geom_ribbon(data=model.AWI,aes(depth,ymin=2*min-mean,ymax=2*max-mean),alpha=0.2) +  
            
            geom_line(aes(x=chronology[,1], y=y.plot.1, color=color.plot.1),alpha=alpha.plot.1) +
            geom_point(aes(x=chronology[,1], y=y.plot.1, color=color.plot.1),alpha=min(alpha.plot.1,0.6),shape="|",size=2) +
            
            geom_line(aes(x=chronology[,1], y=y.plot.2, color=color.plot.2),alpha=alpha.plot.2) +
            geom_point(aes(x=chronology[,1], y=y.plot.2, color=color.plot.2),alpha=min(alpha.plot.2,0.6),shape="|",size=2) +
            
            geom_line(aes(x=chronology[,1], y=y.plot.3, color=color.plot.3),alpha=alpha.plot.3) +
            geom_point(aes(x=chronology[,1], y=y.plot.3, color=color.plot.3),alpha=min(alpha.plot.3,0.6),shape="|",size=2) +
            
            geom_line(aes(x=chronology[,1], y=y.plot.4, color=color.plot.4),alpha=alpha.plot.4) +
            geom_point(aes(x=chronology[,1], y=y.plot.4, color=color.plot.4),alpha=min(alpha.plot.4,0.6),shape="|",size=2) +
            
            geom_line(aes(x=chronology[,1], y=y.plot.5, color=color.plot.5),alpha=alpha.plot.5) +
            geom_point(aes(x=chronology[,1], y=y.plot.5, color=color.plot.5),alpha=min(alpha.plot.5,0.6),shape="|",size=2) +
            
            geom_line(data=model.AWI, aes(x=depth, y=mean, color="meanAgeBP")) +
            geom_point(data=AgeDepthPollen.subset, aes(x=depth, y=AWI.meanAgeBP, color="meanAgeBP"),alpha=0.9,shape="|",size=2) +
            
            # datings:
            
            geom_errorbar(data=calibration, aes(x=depth,ymin=cal.age-cal.age.se,ymax=cal.age+cal.age.se, color="Calibrated.Age"),width=(max(calibration$depth)-min(calibration$depth))*0.01,alpha=1) +
            geom_point(data=calibration, aes(x=depth, y=cal.age, color="Calibrated.Age") , shape=1,alpha=1) + 
            geom_linerange(aes(x=cal.sp[,1],ymin=cal.sp[,2],ymax=cal.sp[,3],color="Calibrated.Age"))+
            
            geom_errorbar(data=calibration.added, aes(x=depth,ymin=cal.age-cal.age.se,ymax=cal.age+cal.age.se, color="Calibrated.Age.green"),width=(max(calibration$depth)-min(calibration$depth))*0.015,alpha=alpha.added) +
            geom_point(data=calibration.added, aes(x=depth, y=cal.age, color="Calibrated.Age.green"),shape=(max(calibration$depth)-min(calibration$depth))*0.015) +    
            
            # legend, theme and color:
            
            labs(title=paste0(AgeDepthPollen.subset$site.name[1]," | ID ",ID," | long ",AgeDepthPollen.subset$long[1], "° lat ",AgeDepthPollen.subset$lat[1],"° | elev ",AgeDepthPollen.subset$elev,"\n","acc.rate  ",round(acc.rate,digits = 2)," a/cm  ",round(1/acc.rate,digits = 4)," cm/a | resolution ",round(if(is.na(parameter.subset$resolution.cm)){(max(AgeDepthPollen.subset$depth,calibration$depth)-parameter.subset$waterline)/parameter.subset$resolution.sections}else{parameter.subset$resolution.cm},digits = 2)," cm | mem.mean ",parameter.subset$memory," | reservoir ",parameter.subset$reservoir," a  ")) +
            labs(x="depth (cm)", y="years BP") +
            scale_color_manual(name="", values=c("Neotoma.5"="#FF00FF","Neotoma.4"="orange","Neotoma.3"="#FFCCFF","Neotoma.2"="yellow","Neotoma.1"="black","meanAgeBP"="blue","Calibrated.Age.green"="#FF00FF","Calibrated.Age"="red"),labels=naming) +
            theme_bw() +
            theme(plot.title=element_text(size=10),
                  plot.margin=margin(1,1,1,1,"cm"),
                  axis.title=element_text(size=10),
                  legend.position="bottom") +
            theme(axis.title.x=element_text(margin=margin(t=15, r=0, b=0, l=0)),
                  axis.text.x=element_text(size=10)) +
            theme(axis.title.y=element_text(angle=90, margin=margin(t=0, r=15, b=00, l=0)),
                  axis.text.y=element_text(size=10)) 
          
          # write png:  
          output
          ggsave(paste0(folder,"/sites/",ID,"/",ID,".png"), width = 297, height = 210, units = "mm", dpi = 300)
          
          # flipped version:
          output.flipped<-output +
            coord_flip (xlim=NULL, ylim=NULL, expand = TRUE, clip="on")
          output.flipped
          ggsave(paste0(folder,"/sites/",ID,"/",ID,".flipped.png"), width = 297, height = 210, units = "mm", dpi = 300) 
          
        }
        
        # ---------------------------------------------------------------------------------------------
        
        # copy plots to new subfolders:
        file.copy(from = paste0(folder,"/sites/",ID,"/",ID,".png"), to = paste0(folder,"/png/",ID,".png"),overwrite = TRUE)
        file.copy(from = paste0(folder,"/sites/",ID,"/",ID,".flipped.png"), to = paste0(folder,"/png.flipped/",ID,".flipped.png"),overwrite = TRUE)
        if(is.na(AgeDepthPollen.subset$core.all.true[1])){file.copy(from = paste0(folder,"/png/",ID,".png"), to = paste0(folder,"/sites.without.neotoma.model/",ID,".png"),overwrite = TRUE)}else{
          if(AgeDepthPollen.subset$core.all.true[1]==FALSE){file.copy(from = paste0(folder,"/png/",ID,".png"), to = paste0(folder,"/sites.outside.min.max/",ID,".png"),overwrite = TRUE)}
          if(AgeDepthPollen.subset$core.all.true.2[1]==FALSE){file.copy(from = paste0(folder,"/png/",ID,".png"), to = paste0(folder,"/sites.outside.min.max_2/",ID,".png"),overwrite = TRUE)}
        }
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
      ### -----complete ID run-----    
      
      {
        
        # delete the ID from the error list:
        bacon.error<-bacon.error[-which(bacon.error==ID)] 
        write.csv2(bacon.error,"bacon_error.csv",row.names = FALSE)
        
        # delete the ID from the list for the outer loop and go to the next ID:
        parameter.ID<-parameter.ID[-1]
        write.csv2(parameter.ID,"bacon_ID.csv",row.names = FALSE)
        
      }
      
      # ---------------------------------------------------------------------------------------------
      
    }
    
  } 
  
  # ---------------------------------------------------------------------------------------------
  
} # end loop run for every ID:


