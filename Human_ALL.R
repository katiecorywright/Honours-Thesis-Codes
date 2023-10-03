

setwd("/HANGER/katie/raw_viromescan_bbmap_human/final")


  .libPaths("/home/aaron/R/x86_64-pc-linux-gnu-library/3.6/")
library(knitr)
library(dplyr)
library(ggplot2)
library(DESeq2)


save.image("Final_filtered_files.RData")


#####filtering viruses that contained 5 or less read counts

filenames<-list.files(path="/HANGER/katie/raw_viromescan_bbmap_human/coverage_human_all/",pattern="*_coverage.text$", full.names=TRUE)


path<-"/HANGER/katie/envs/viromescan/viromescan/var/HumanALLcomplete.txt"
    IDs2=read.delim(path, header = F, sep = ';')

  IDs2$ViralGenomes<-unlist(sapply(strsplit(IDs2$V3,"- *"),tail,1))
    colnames(IDs2)[1:3]=c('Family', 'Genera','Species')

    for (i in 1:length(filenames))
    {
      oneFile <- filenames[i]
      oneFileData <- as.data.frame(read.delim(oneFile, sep="\t", header=T))
    oneFileData$Read_depth<-oneFileData$Plus_reads + oneFileData$Minus_reads

      #extract virus name and aligned columns
     oneFileDataSmall <- oneFileData[,c("X.ID","Read_depth")]

      x2 <- unlist(strsplit(oneFile, "/"))
      #get only file name - that last item
      oneFileName <- x2[length(x2)]
    oneFileName = gsub(x= oneFileName, "_coverage.text", "", fixed = T)
      colnames(oneFileDataSmall) <- c("ViralGenomes",oneFileName ) 


         names<-unlist(sapply(strsplit(oneFileDataSmall[,1],"NC_"),tail,1))
    names<-gsub("|","",names,fixed=T)
    names<-paste0("NC",sep="_",names)
    oneFileDataSmall$ViralGenomes<-names


  IDs2 <- cbind(IDs2, oneFileDataSmall[match(IDs2$ViralGenomes, oneFileDataSmall$ViralGenomes), 2,drop=FALSE])

oneFileDataSmall=NULL
names=NULL
oneFileName=NULL
    }


head(IDs2)

rownames(IDs2)<-IDs2$ViralGenomes
IDs2[is.na(IDs2)] <- as.numeric("0")

IDs2<-IDs2[order(rowSums(IDs2[,5:length(IDs2)]),decreasing = TRUE),]

IDs<-IDs2[,1:4]

df<-IDs2[,5:length(IDs2)]

df2<-mutate_all(df,function(x) (as.numeric(as.character(x))))

IDs3<-apply(df2, 2, function(x) ifelse(x<5, 0, x))

IDs3<-IDs3[rowSums(IDs3)>0,]

passed_RD<-IDs3


######passed_RD is our final list of viruses that passed our read depth filter


#### we are now going to obtain the FPKMs for these for plotting abundance






filenames<-list.files(path="/HANGER/katie/raw_viromescan_bbmap_human/rpkm_human_all",pattern="*_rpkm.text$", full.names=TRUE)


path<-"/HANGER/katie/envs/viromescan/viromescan/var/HumanALLcomplete.txt"
    IDs2=read.delim(path, header = F, sep = ';')


    colnames(IDs2)[1:3]=c('Family', 'Genera','Species')

    for (i in 1:length(filenames))
    {
      oneFile <- filenames[i]
      oneFileData <- as.data.frame(read.delim(oneFile, sep="\t", skip=4, header=T))
    
      #extract virus name and aligned columns
      oneFileDataSmall <- oneFileData[,c(1,8)]
      x2 <- unlist(strsplit(oneFile, "/"))
      #get only file name - that last item
      oneFileName <- x2[length(x2)]
    oneFileName = gsub(x= oneFileName, "_rpkm.text", "", fixed = T)
      colnames(oneFileDataSmall) <- c("ViralGenomes",oneFileName ) 


         names<-unlist(sapply(strsplit(oneFileDataSmall[,1],"NC_"),tail,1))
    names<-gsub("|","",names,fixed=T)
    names<-paste0("NC",sep="_",names)
    oneFileDataSmall$ViralGenomes<-names


  IDs2 <- cbind(IDs2, oneFileDataSmall[match(IDs2$ViralGenomes, oneFileDataSmall$ViralGenomes), 2,drop=FALSE])

oneFileDataSmall=NULL
names=NULL
oneFileName=NULL
    }


head(IDs2)

rownames(IDs2)<-IDs2$ViralGenomes
IDs2[is.na(IDs2)] <- as.numeric("0")

IDs2<-IDs2[order(rowSums(IDs2[,5:length(IDs2)]),decreasing = TRUE),]

IDs<-IDs2[,1:4]

df<-IDs2[,5:length(IDs2)]



## At this point df=dataframe containing all samples, IDs is Family, Genus, species names and NCBI codes #####################

###filtering according to DESeq2 RNA expression ######
  metadata<- read.csv("/HANGER/katie/metadata_infant_mothers_meds.csv",header=TRUE)
metadata2<-metadata[match(colnames(df),metadata$sample_id),]
df2<-mutate_all(df,function(x) as.integer(as.numeric(as.character(x))))

#####removing all viral rows that are not detected in any sample##

df2<-df2[rowSums(df2)>0,]

    sapply(df2,class)
###filtering viruses to only include those that pass read depth filterting 

df_final_filtered<-df2[match(rownames(passed_RD),rownames(df2)),]
    write.csv(df_final_filtered, "Human_all_abundance.csv", row.names = T)

######frequency Calculation 
final_freq<-df_final_filtered

isnum <- sapply(final_freq, is.numeric)
final_freq[isnum] <- lapply(final_freq[isnum], function(x) ifelse(x == 0, as.numeric(0), as.numeric(1)))
final_freq
    final_freq$sums<-rowSums(final_freq[,2:ncol(final_freq)])
    final_freq<final_freq[order(final_freq$sums,decreasing = TRUE),]
head(final_freq)
    write.csv(final_freq, "Human_all_freq.csv", row.names = T)



final_tax<-IDs[match(rownames(df_final_filtered), rownames(IDs)),]

mothers_metadata<-metadata2[metadata2$PHASE=="MOTHER",]
mother_codes<-mothers_metadata$sample_id
infants_metadata<- metadata2[metadata2$PHASE=="INFANT",]
infant_codes<- infants_metadata$sample_id

###overwriting the abundane dataframe with the frequency one to save having to write codes twice####

setwd("/HANGER/katie/raw_viromescan_bbmap_human/final/frequency")
final_filtered<-df_final_filtered 
df_final_filtered<-final_freq
df_final_filtered_infants<-df_final_filtered[,match(infant_codes,colnames(df_final_filtered))]
df_final_filtered_infants$sums<-rowSums(df_final_filtered_infants)

df_final_filtered_mothers<-df_final_filtered[,match(mother_codes,colnames(df_final_filtered))]
df_final_filtered_mothers$sums<-rowSums(df_final_filtered_mothers)

######

####
##########################calculating sample abundance and ###graphing################infant
perc2<-NULL
perc<-NULL

df_final_filtered_infants<-as.data.frame(df_final_filtered_infants)
 perc=c()
  for (i in 1 : (length(rownames(df_final_filtered_infants))))
  {
    perc[i]=round(df_final_filtered_infants$sums[i]/sum(df_final_filtered_infants$sums),4)
  }


 perc2=as.matrix(perc)
  rownames(perc2)=rownames(df_final_filtered_infants)
  colnames(perc2)=c('')
  df_final_filtered_infants=as.matrix(df_final_filtered_infants)
  colnames(perc2)="Detected_frequency"
  write.table(perc2, file='infants_KT_filtered_percentages%.txt', sep='\t')
  write.table(df_final_filtered_infants, file='-infants_KT_filtered_sums.txt', sep='\t')

  filt=(perc2)<0.005
  perc2=as.data.frame(perc2[!filt,])
  table(filt)
  ##33 low abundance were filtered out, 52 were kept###
  
  infant_frequency<-as.data.frame(perc2)
  colnames(infant_frequency)<-"Detected_frequency"
  low_abundance=as.data.frame(round(1-(sum(infant_frequency)),2))
  colnames(low_abundance)<-"Detected_frequency"
  low_abundance$Species<-"Low abundance"
  infant_frequency$Species<-IDs[match(rownames(infant_frequency),rownames(IDs)),"Species"]

 infant_frequency$Species<-gsub("- .*","", infant_frequency$Species)

   infant_frequency=rbind(infant_frequency, low_abundance)
infant_frequency$Group<-"Infants"


#########################mothers##############


perc2<-NULL
perc<-NULL

df_final_filtered_mothers<-as.data.frame(df_final_filtered_mothers)
 perc=c()
  for (i in 1 : (length(rownames(df_final_filtered_mothers))))
  {
    perc[i]=round(df_final_filtered_mothers$sums[i]/sum(df_final_filtered_mothers$sums),4)
  }




 perc2=as.matrix(perc)
  rownames(perc2)=rownames(df_final_filtered_mothers)
  colnames(perc2)=c('')
  df_final_filtered_mothers=as.matrix(df_final_filtered_mothers)
  colnames(perc2)="Detected_frequency"
  write.table(perc2, file='mothers_KT_filtered_percentages%.txt', sep='\t')
  write.table(df_final_filtered_mothers, file='mothers_KT_filtered_sums.txt', sep='\t')

  filt=(perc2)<0.005
  perc2=as.data.frame(perc2[!filt,])
  table(filt)
  ##33 low abundance were filtered out, 52 were kept###
  
mothers_frequency<-as.data.frame(perc2)
  colnames(mothers_frequency)<-"Detected_frequency"
  low_abundance=as.data.frame(round(1-(sum(mothers_frequency)),2))
  colnames(low_abundance)<-"Detected_frequency"
  low_abundance$Species<-"Low abundance"
  mothers_frequency$Species<-IDs[match(rownames(mothers_frequency),rownames(IDs)),"Species"]

 mothers_frequency$Species<-gsub("- .*","", mothers_frequency$Species)

   mothers_frequency=rbind(mothers_frequency, low_abundance)
mothers_frequency$Group<-"Mothers"


freq_table<-rbind(mothers_frequency, infant_frequency)

freq_table$Species[c(2,14)]<-"Hepatitis C virus 2"

write.csv(freq_table, file="Human_abundance_filtered_final_table.csv")


ggplot(freq_table,aes( fill=Species,x=Group,y=Detected_frequency)) + geom_bar(position="stack", stat="identity", colour="white") + xlab(paste("Group")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                        panel.background = element_blank() ,
                        panel.grid.major = element_blank() , 
                        panel.grid.minor = element_blank(), 
                        panel.border = element_blank(), 
                        axis.text.x = element_text(color = 'black', size=15), 
                        axis.text.y = element_text(color = 'black', size=15),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())

ggsave(filename ="KT_filt_species_frequency.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()




    #3###these were filtering that wasn't suitable###
dds <- DESeqDataSetFromMatrix(countData = df3, colData = metadata2, design = NULL)
dds<-DESeq(dds)
dds<- estimateSizeFactors(dds)

#######could not filter using DESEq2 as the data is too sparce (too many genes cantain zero counts)

###trying EdgeR for filtering

library("edgeR")
keep.exprs<-filterByExpr(df3,group=c(metadata2$PHASE))



######################################################################################################################################################################################## WHAT WE DID TODAY (28/9/2023)


setwd("/HANGER/katie/raw_viromescan_bbmap_human/final/")
load("/HANGER/katie/raw_viromescan_bbmap_human/final/Final_filtered_files.RData")

###our taxonomy we need is final_tax## 
##3 our datatframe is df_final_filtered##
##metadata is metadata2####


library(phyloseq)
require("ape")


OTU=otu_table(df_final_filtered,taxa_are_rows=TRUE)

TAX = tax_table(as.matrix(final_tax))

rownames(metadata2)<-metadata2$sample_id
physeq=phyloseq(OTU,TAX)
sd<-sample_data(metadata2)
#plot_bar(physeq1,fill="Family")

random_tree<-rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

physeq1<-merge_phyloseq(physeq,sd, random_tree)
head(otu_table(physeq1))


xdf<-df_final_filtered[order(rowSums(df_final_filtered),decreasing=TRUE),]
head(xdf)

top_four<-rownames(xdf)[1:4]
species_names<-final_tax[match(top_four,final_tax[,"ViralGenomes"]),"Species"]

species_names<-gsub("- .*","",species_names)

species_names[2]<-"Hepatitis C Virus 2"

xdf_infants<-xdf[,metadata2$PHASE=="INFANT",]
sum(xdf_infants)
#2920209

xdf_infants<-xdf_infants[top_four,]

dim(xdf_infants)

rownames(xdf_infants)<-species_names

xdf_infants<-as.data.frame(t(xdf_infants))
xdf_infants$PHASE<-as.character("Infants")


xdf_mother<-xdf[,metadata2$PHASE=="MOTHER",]
sum(xdf_mother)
#2327649

xdf_mother<-xdf_mother[top_four,]

dim(xdf_mother)
head(xdf_mother)

rownames(xdf_mother)<-species_names

xdf_mother<-as.data.frame(t(xdf_mother))
xdf_mother$PHASE<-as.character("Mothers")

#######calculating the  contribution of each individuals virus towards the total relative abudance (overall). 

#sum(xdf_mother[,1:4])
file=xdf_mother
perc=c()
for (a in 1:4) {
  j=a
  
  for (i in 1 : (length(rownames(file))))
  {
    perc[i]=round(file[i,j]/2327649,6)*100
  }
  file[,j]<-as.data.frame(perc)
  
}
xdf_mother<-file

file=xdf_infants
perc=c()
for (a in 1:4) {
  j=a
  
  for (i in 1 : (length(rownames(file))))
  {
    perc[i]=round(file[i,j]/2920209,4)*100
  }
  file[,j]<-as.data.frame(perc)
  
}


xdf_infants<-file

total<-rbind(xdf_mother,xdf_infants)

A<- ggplot(total,aes(x=PHASE,y=total[,1],color=PHASE)) + geom_boxplot() + geom_jitter() + ggtitle(colnames(total[1]))+
  scale_color_manual("PHASE", values=c("salmon", "skyblue")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                                   panel.background = element_blank() ,
                                   panel.grid.major = element_blank() , 
                                   panel.grid.minor = element_blank(), 
                                   panel.border = element_blank(), 
                                   axis.text.x = element_text(color = 'black', size=15), 
                                   axis.text.y = element_text(color = 'black', size=15),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   plot.title = element_text(hjust = 0.5))
R<- ggplot(total,aes(x=PHASE,y=total[,2],color=PHASE)) + geom_boxplot() + geom_jitter() + ggtitle(colnames(total[2]))+
  scale_color_manual("PHASE", values=c("salmon", "skyblue")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                                   panel.background = element_blank() ,
                                   panel.grid.major = element_blank() , 
                                   panel.grid.minor = element_blank(), 
                                   panel.border = element_blank(), 
                                   axis.text.x = element_text(color = 'black', size=15), 
                                   axis.text.y = element_text(color = 'black', size=15),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   plot.title = element_text(hjust = 0.5))
O<- ggplot(total,aes(x=PHASE,y=total[,3],color=PHASE)) + geom_boxplot() + geom_jitter() + ggtitle(colnames(total[3]))+
  scale_color_manual("PHASE", values=c("salmon", "skyblue")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                                   panel.background = element_blank() ,
                                   panel.grid.major = element_blank() , 
                                   panel.grid.minor = element_blank(), 
                                   panel.border = element_blank(), 
                                   axis.text.x = element_text(color = 'black', size=15), 
                                   axis.text.y = element_text(color = 'black', size=15),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   plot.title = element_text(hjust = 0.5))
N<- ggplot(total,aes(x=PHASE,y=total[,4],color=PHASE)) + geom_boxplot() + geom_jitter() + ggtitle(colnames(total[4]))+
  scale_color_manual("PHASE", values=c("salmon", "skyblue")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                                   panel.background = element_blank() ,
                                   panel.grid.major = element_blank() , 
                                   panel.grid.minor = element_blank(), 
                                   panel.border = element_blank(), 
                                   axis.text.x = element_text(color = 'black', size=15), 
                                   axis.text.y = element_text(color = 'black', size=15),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   plot.title = element_text(hjust = 0.5))



library(ggpubr)

ggarrange(A, R, O, N , 
          labels = c("A)", "B)", "C)","D)"),
          ncol = 2, nrow = 2)


ggsave(filename ="individuals_contribution_towards_total_abundance.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()





