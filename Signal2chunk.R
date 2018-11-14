####### Parameters #######

args=commandArgs(trailingOnly = TRUE)

spanning_reads_file=args[1]
theo_squiggle_file=args[2]

if(file.exists(spanning_reads_file)==F){stop("Spanning_reads file not found")}
if(file.exists(theo_squiggle_file)==F){stop("Reference_squiggle file not found")}


mvnStepPattern=25
digit=1

####### Libraries #######
library(rhdf5)
library(dtw,quietly = T)
library(ggplot2)
library(dplyr,quietly = T)
library(tidyr,quietly = T)

####### Functions #######
fast5_to_current=function(f5){
a=h5ls(f5)
b=as.integer(h5read(f5,a$group[grep("Read_",a$group)])[[1]])
d=h5readAttributes(f5,"/UniqueGlobalKey/channel_id")

df=data.frame(number=1:length(b),signal=b)
df$sec=df$number / d$sampling_rate
#df$signalz=as.numeric(scale(df$signal))

H5close()

df
}

####### Load sequence read data #######
sr=read.table(spanning_reads_file,header=T,sep="\t",stringsAsFactors = F)
sr=unique(sr[order(sr$name),])

sr_squiggle=setNames(lapply(sr$path,fast5_to_current),sr$name)


####### Process theoretical squiggles #######

index=scan(paste(theo_squiggle_file,".index",sep=""),"character",quiet=T)
index_split=strsplit(index,":#|;")
index_split[["extra"]]=-1

theo_data=list()

for(i in 1:length(index)){
x=read.table(theo_squiggle_file,
             skip=as.integer(index_split[[i]][1]),
             nrows=as.integer(index_split[[i+1]][1])-as.integer(index_split[[i]][1])-2,
             header=T,
             stringsAsFactors = F,
             sep="\t")

x$currentz=as.numeric(scale(x$current))
theo_data[[index_split[[i]][2]]]=x
}

theo_meta=data.frame(start=as.integer(sapply(index_split,function(x) x[1])),
                    theo_name=sapply(index_split,function(x) x[2]),
                    strand=gsub("left_|right_|flank_|_plus|TR_|sequence_","",sapply(index_split,function(x) x[2])),
                    flank_size=as.integer(sapply(index_split,function(x) gsub("flank_size=","",x[grep("flank_size=",x)]))),
                    theo_squiggle_length=c(sapply(theo_data,nrow),NA),
                    TR_number=as.integer(sapply(index_split,function(x) gsub("TR_number=","",x[grep("TR_number=",x)]))),
                    TR_number_plus=as.integer(sapply(index_split,function(x) gsub("TR_number_plus=","",x[grep("TR_number_plus=",x)]))),
                    stringsAsFactors = F)

theo_meta$TR_unit_length=(theo_meta$theo_squiggle_length-theo_meta$flank_size)/theo_meta$TR_number
theo_meta[is.na(theo_meta$flank_size),"TR_unit_length"]=theo_meta[is.na(theo_meta$flank_size),"theo_squiggle_length"]/theo_meta[is.na(theo_meta$flank_size),"TR_number_plus"]

for(i in c("TR_sequence_positive","TR_sequence_negative")){
  x=theo_data[[i]]
  a=theo_meta[theo_meta$theo_name==i & is.na(theo_meta$theo_name)==F,"TR_unit_length"]
  b=(theo_meta[theo_meta$theo_name==i & is.na(theo_meta$theo_name)==F,"TR_number_plus"]-theo_meta[theo_meta$theo_name==i & is.na(theo_meta$theo_name)==F,"TR_number"])/2
  x2=x[-c(1:(b*a),(nrow(x)-(b*a)+2):nrow(x)),]
  theo_data[[i]]=x2
}


##### Digits per TR unit to select #####

a=theo_meta[is.na(theo_meta$TR_number)==F,c("theo_name","flank_size","TR_number","TR_unit_length","theo_squiggle_length")]

sep_digits=mapply(function(v,w,x,y,z){
  if(is.na(w)){
    seq(digit,digit+(x)*y,y)
  } else if(grepl("left",v)){
    seq(w+digit,w+digit+(x)*y,y)
  } else if(grepl("right",v)){
    z-seq(w+digit,w+digit+(x)*y,y)
  }
},a[,1],a[,2],a[,3],a[,4],a[,5],SIMPLIFY = F)



######Reduce the amount of repeat units in TR_sequence theo's and seps######

rep_units=5

theo_data_less=setNames(lapply(c("TR_sequence_positive","TR_sequence_negative"),function(x) {
  y=theo_data[[x]]
  z=theo_meta[theo_meta$theo_name==x & is.na(theo_meta$theo_name)==F,"TR_unit_length"]
  y[y$pos <= min(y$pos)+rep_units * z,]
}),c("TR_sequence_positive","TR_sequence_negative"))

for(theo_name in c("TR_sequence_positive","TR_sequence_negative")){
  sep_digits[[theo_name]]=sep_digits[[theo_name]][1:(rep_units+1)]
}

####### Theoretical squiggle VS read DTW #######

##### Round 1: Delineating the TR within the read #####

### Divide each squiggle in overlapping chunks with chunk-based scaling ###

max_flank_size=max(theo_meta$flank_size,na.rm = T)
window_size=max_flank_size*mvnStepPattern*2

sr_squiggle_win1=lapply(sr_squiggle,function(x){
  split(x,cut(x$number,unique(c(seq(0,max(x$number),by=window_size),max(x$number)))))
})

sr_squiggle_win2=lapply(sr_squiggle,function(x){
  split(x,cut(x$number,unique(c(seq(-window_size/2,max(x$number)-window_size/2,by=window_size),max(x$number)))))
})

keys <- unique(c(names(sr_squiggle_win1), names(sr_squiggle_win2)))
sr_squiggle_win=setNames(mapply(c, sr_squiggle_win1[keys], sr_squiggle_win2[keys]), keys)

sr_squiggle_win=lapply(sr_squiggle_win,function(x) lapply(x,function(y){
  signalz=as.numeric(scale(y$signal))
  cbind(y,signalz)}))

sr_squiggle_win=mapply(function(x,y) x[y==FALSE],sr_squiggle_win,sapply(sr_squiggle_win,function(x) sapply(x,nrow) < max_flank_size+1))


### Start DTW alignment and record distances ###

name_list=sapply(sr_squiggle_win,names)

dtw_master=data.frame(name=rep(names(name_list),times=sapply(name_list,length)),window=unlist(name_list),stringsAsFactors = F)

dtw_master2=left_join(dtw_master,left_join(sr[,c("name","strand")],theo_meta[is.na(theo_meta$flank_size)==F & is.na(theo_meta$TR_number),c("strand","theo_name")]))

dtw_master2$distance=apply(dtw_master2,1,function(x){
  print(paste(x[1],x[2]))
  dtw(theo_data[[x[4]]]$currentz,
      sr_squiggle_win[[x[1]]][[x[2]]]$signalz,
      open.end = T,
      open.begin = T,
      step.pattern = mvmStepPattern(mvnStepPattern),
      distance.only = T)$distance
})

best_matches=as.data.frame(dtw_master2 %>% group_by(name,theo_name) %>% summarise (min_dist=min(distance),window=window[which(distance==min(distance))[1]]),stringsAsFactors=F)

best_matches[grep("left_flank",best_matches$theo_name),"number"]=sapply(strsplit(gsub("\\(|\\]","",best_matches[grep("left_flank",best_matches$theo_name),"window"]),","),function(x) as.numeric(x[1]))

best_matches[grep("right_flank",best_matches$theo_name),"number"]=sapply(strsplit(gsub("\\(|\\]","",best_matches[grep("right_flank",best_matches$theo_name),"window"]),","),function(x) as.numeric(x[2]))

best_matches$flank=gsub("_negative|_positive","",best_matches$theo_name)


### Select region of interest in each squiggle ###
best_matches=best_matches[duplicated(best_matches[,c("name","flank")])==F,]
 
best_matches_selector=spread(best_matches[,c("name","flank","number")],flank,number)

best_matches_selector=best_matches_selector[best_matches_selector$left_flank < best_matches_selector$right_flank,]

sr_squiggle_sel=mapply(
  function(x,y,z) {
    xxx=x[x$number > y & x$number < z,]
    signalz=as.numeric(scale(xxx$signal))
    cbind(xxx,signalz)
    }
  ,
  sr_squiggle[best_matches_selector$name],best_matches_selector$left_flank,best_matches_selector$right_flank,SIMPLIFY = F
  )

### Clean-up ###

rm(sr_squiggle_win1,sr_squiggle_win2,sr_squiggle_win,keys,name_list,dtw_master,best_matches_selector)


##### Round 2: Divide squiggle into chunks #####

### Start and stops ###

dtw_master_round2=left_join(data.frame(name=names(sr_squiggle_sel),stringsAsFactors = F),left_join(sr[,c("name","strand")],theo_meta[is.na(theo_meta$flank_size)==F & is.na(theo_meta$TR_number)==F,c("strand","theo_name")]))


chunk_sep_list=apply(dtw_master_round2[,c("name","theo_name")],1,function(x){
  print(paste(x[1],x[2]));
  alignment=dtw(theo_data[[x[2]]]$currentz,
      sr_squiggle_sel[[x[1]]]$signalz,
      open.end = T,
      open.begin = T,
      step.pattern = mvmStepPattern(mvnStepPattern));
  data.frame(name=x[1],theo_name=x[2],index1=alignment$index1[alignment$index1 %in% sep_digits[[x[2]]][1]],index2=alignment$index2[alignment$index1 %in% sep_digits[[x[2]]][1]],number=sr_squiggle_sel[[x[1]]]$number[alignment$index2[alignment$index1 %in% sep_digits[[x[2]]][1]]],   distance=alignment$distance, normdist=alignment$normalizedDistance, stringsAsFactors = F,row.names = NULL);
})

chunk_sep=do.call("rbind",chunk_sep_list)

chunk_sep_meta=data.frame(name=unique(chunk_sep$name),stringsAsFactors = F)
chunk_sep_meta=cbind(chunk_sep_meta,do.call("rbind",lapply(chunk_sep_meta$name,function(x) data.frame(start=min(chunk_sep[chunk_sep$name==x,"number"]),end=max(chunk_sep[chunk_sep$name==x,"number"])))))

chunk_sep_meta2=chunk_sep_meta[chunk_sep_meta$end-chunk_sep_meta$start > nrow(theo_data$TR_sequence_positive),]
sr_squiggle_sel2=sr_squiggle_sel[chunk_sep_meta2$name]

### Center piece ###

# Create new selection for unbiased scaling #

sr_squiggle_sel_center=mapply(
  function(x,y,z) {
    xxx=x[x$number > y & x$number < z,]
    xxx$signalz=as.numeric(scale(xxx$signal))
    xxx
    }
  ,
  sr_squiggle_sel2,chunk_sep_meta2$start,chunk_sep_meta2$end,SIMPLIFY = F
  )

window_size_round2=round(nrow(theo_data$TR_sequence_positive)/450*4000*1.3*2)

dtw_master_center=data.frame(name=names(sr_squiggle_sel_center),stringsAsFactors = F)

dtw_master_center2=left_join(dtw_master_center,left_join(sr[,c("name","strand")],theo_meta[grep("TR_sequence",theo_meta$theo_name),c("strand","theo_name")]))

# Run DTW #

chunk_sep_center=data.frame()

for(i in 1:nrow(dtw_master_center2)){
  a=data.frame(sr_squiggle_sel_center[[i]],row.names = NULL,stringsAsFactors = F)
  theo_name=dtw_master_center2[i,"theo_name"]
  
  left2right=1
  right2left=nrow(a)
  seps=sep_digits[[theo_name]][-c(1,length(sep_digits[[theo_name]]))]

  print(paste(dtw_master_center2[i,"name"],dtw_master_center2[i,"theo_name"]))

  while(left2right < right2left){
    b=data.frame(a[left2right:(left2right+window_size_round2),],row.names = NULL,stringsAsFactors = F)
  
    if(right2left-window_size_round2 > 0){
      b2=data.frame(a[(right2left-window_size_round2):right2left,],row.names = NULL,stringsAsFactors = F)
    } else {
      b2=data.frame(a[0:right2left,],row.names = NULL,stringsAsFactors = F)
    }
    alignment=dtw(
      theo_data_less[[theo_name]]$currentz,
      b$signalz,
      open.end = T,
      open.begin = F,
      step.pattern = mvmStepPattern(mvnStepPattern))

    alignment2=dtw(
      theo_data_less[[theo_name]]$currentz,
      b2$signalz,
      open.end = F,
      open.begin =T,
      step.pattern = mvmStepPattern(mvnStepPattern))

    chunk_sep_tmp=data.frame(name=dtw_master_center2[i,"name"],
                             theo_name=dtw_master_center2[i,"theo_name"],
                             index1=alignment$index1[alignment$index1 %in% seps[1:ceiling(length(seps)/2)]],
                             index2=alignment$index2[alignment$index1 %in% seps[1:ceiling(length(seps)/2)]],
                             number=b$number[alignment$index2[alignment$index1 %in% seps[1:ceiling(length(seps)/2)]]],
                             distance=alignment$distance,
                             normdist=alignment$normalizedDistance,
                             direction="left2right",
                             stringsAsFactors = F,row.names = NULL)
   
    chunk_sep_tmp2=data.frame(name=dtw_master_center2[i,"name"],
                             theo_name=dtw_master_center2[i,"theo_name"],
                             index1=alignment2$index1[alignment2$index1 %in% seps[ceiling(length(seps)/2):length(seps)]],
                             index2=alignment2$index2[alignment2$index1 %in% seps[ceiling(length(seps)/2):length(seps)]],
                             number=b2$number[alignment2$index2[alignment2$index1 %in% seps[ceiling(length(seps)/2):length(seps)]]],
                             distance=alignment2$distance,
                             normdist=alignment2$normalizedDistance,
                             direction="right2left",
                             stringsAsFactors = F,row.names = NULL)
  
    left2right=left2right+alignment$index2[alignment$index1 == seps[ceiling(length(seps)/2)]]
    right2left=right2left-(nrow(b2)-alignment2$index2[alignment2$index1 == seps[ceiling(length(seps)/2)]])
    

    chunk_sep_center=rbind(chunk_sep_center,chunk_sep_tmp,chunk_sep_tmp2)
  }
}

# Remove overlapping data points #

chunk_sep_center2_list=lapply(dtw_master_center2$name,function(x){
  y=chunk_sep_center[chunk_sep_center$name==x,]
  left_max=max(y[y$direction=="left2right","number"])
  right_min=min(y[y$direction=="right2left","number"])
  y_overlap=y[(y$direction=="left2right" & y$number >= right_min) | (y$direction=="right2left" & y$number <= left_max),]
  left_dist=mean(y_overlap[y_overlap$direction=="left2right","distance"])
  right_dist=mean(y_overlap[y_overlap$direction=="right2left","distance"])
  if(right_dist<left_dist){
    exclude="left2right"
  } else {
    exclude="right2left"
  }
  y2=y[!(y$direction==exclude & y$number %in% y_overlap[y_overlap$direction==exclude,"number"]),]
})

chunk_sep_center2=do.call("rbind",chunk_sep_center2_list)

### Merge chunk data and generate QC data ###

chunk_sep_total=rbind(chunk_sep,chunk_sep_center2[,-ncol(chunk_sep_center2)])
chunk_sep_total=chunk_sep_total[order(chunk_sep_total$name,chunk_sep_total$number),]
chunk_sep_total2=left_join(chunk_sep_total,sr[,c("name","strand")])

diff_list=sapply(dtw_master_center2$name,function(x){
  diff(chunk_sep_total[chunk_sep_total$name==x,"number"])
})

chunk_sep_diff=data.frame(name=rep(names(diff_list),times=sapply(diff_list,length)),difference=unlist(diff_list),stringsAsFactors = F)
chunk_sep_diff2=left_join(chunk_sep_diff,sr[,c("name","strand")])


chunk_sep_rep_units=setNames(as.data.frame(table(chunk_sep_total$name),stringsAsFactors = F),c("name","repeat_units"))
chunk_sep_rep_units2=left_join(chunk_sep_rep_units,sr[,c("name","strand")])

qc_center=as.data.frame(chunk_sep_center2 %>% group_by(name) %>% summarise(avg_center_normdist=mean(normdist),avg_center_dist=mean(distance)),stringsAsFactors=F)

qc_flank=as.data.frame(chunk_sep %>% group_by(name) %>% summarise(avg_flank_normdist=mean(normdist),avg_flank_dist=mean(distance)),stringsAsFactors=F)

chunk_sep_rep_units3=left_join(chunk_sep_rep_units2,qc_center)
chunk_sep_rep_units3=left_join(chunk_sep_rep_units3,qc_flank)

### Write chunks and metadata ###

folder1=gsub("\\..*","",strsplit(theo_squiggle_file,"/")[[1]][length(strsplit(theo_squiggle_file,"/")[[1]])])
folder2=gsub("\\..*","",strsplit(spanning_reads_file,"/")[[1]][length(strsplit(spanning_reads_file,"/")[[1]])])
write_dir=paste(getwd(),"/chunks",folder1,folder2,sep="/")

dir.create(write_dir,recursive = T)

for(i in 1:(nrow(chunk_sep_total2)-1)){
  if(chunk_sep_total2[i,"name"]==chunk_sep_total2[i+1,"name"]){
    a=sr_squiggle_sel_center[[chunk_sep_total2[i,"name"]]]
    b=a[a$number >= chunk_sep_total2[i,"number"] & a$number <= chunk_sep_total2[i+1,"number"],]
    write.table(b[,c("number","signal")],paste(write_dir,"/",chunk_sep_total2[i,"name"],"_",chunk_sep_total2[i,"strand"],"_",i,".chunk",sep=""),sep="\t",row.names = F,eol="\n",quote=F)
  }
}

write.table(chunk_sep_rep_units3,paste(write_dir,"/read_metadata.table",sep=""),sep="\t",eol="\n",row.names = F,quote=F)


### Generate report ###
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

rmarkdown::render(paste0(script.basename,"/scripts/Generated_Chunk_report.Rmd"),params=list(sr_squiggle_sel=sr_squiggle_sel,chunk_sep_total2=chunk_sep_total2,title=paste(folder1,folder2),chunk_sep_diff2=chunk_sep_diff2,chunk_sep_rep_units2=chunk_sep_rep_units2),output_file = paste(write_dir,"/chunk_report.html",sep=""),quiet=TRUE)
