setwd("/HANGER/katie/raw_viromescan_bbmap_virus/final/deff_viral/")


  .libPaths("/home/aaron/R/x86_64-pc-linux-gnu-library/3.6/")
library(knitr)
library(dplyr)
library(ggplot2)
library(DESeq2)





#####filtering viruses that contained 5 or more read counts

filenames<-list.files(path="/HANGER/katie/raw_viromescan_bbmap_virus/coverage/",pattern="*_coverage.text$", full.names=TRUE)


path<-"/HANGER/katie/envs/viromescan/viromescan/var/VirusALLcomplete.txt"
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





filenames<-list.files(path="/HANGER/katie/raw_viromescan_bbmap_virus/rpkm",pattern="*_rpkm.text$", full.names=TRUE)


path<-"/HANGER/katie/envs/viromescan/viromescan/var/VirusALLcomplete.txt"
    IDs2=read.delim(path, header = F, sep = ';')

  IDs2$ViralGenomes<-unlist(sapply(strsplit(IDs2$V3,"- *"),tail,1))
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
df2<-df2[rowSums(df2)>0,]

    sapply(df2,class)

df_final_filtered<-df2[match(rownames(passed_RD),rownames(df2)),]
  write.csv(df_final_filtered, "Virus_all_abundance.csv", row.names = T)

######frequency Calculation 
final_freq<-df_final_filtered

isnum <- sapply(final_freq, is.numeric)
final_freq[isnum] <- lapply(final_freq[isnum], function(x) ifelse(x == 0, as.numeric(0), as.numeric(1)))
final_freq
    final_freq$sums<-rowSums(final_freq[,2:ncol(final_freq)])
    final_freq<final_freq[order(final_freq$sums,decreasing = TRUE),]
head(final_freq)
    write.csv(final_freq, "Virus_all_freq.csv", row.names = T)


final_tax<-IDs[match(rownames(df_final_filtered), rownames(IDs)),]

mothers_metadata<-metadata2[metadata2$PHASE=="MOTHER",]
mother_codes<-mothers_metadata$sample_id
infants_metadata<- metadata2[metadata2$PHASE=="INFANT",]
infant_codes<- infants_metadata$sample_id

setwd("/HANGER/katie/raw_viromescan_bbmap_virus/final/frequency")
final_filtered<-df_final_filtered 
df_final_filtered<-final_freq

df_final_filtered_infants<-df_final_filtered[,match(infant_codes,colnames(df_final_filtered))]
df_final_filtered_infants$sums<-rowSums(df_final_filtered_infants)


df_final_filtered_mothers<-df_final_filtered[,match(mother_codes,colnames(df_final_filtered))]
df_final_filtered_mothers$sums<-rowSums(df_final_filtered_mothers)

######

####
##########################calculating sample abundance and ###graphing################infanlllllll
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
  write.table(perc2, file='V-infants_KT_filtered_percentages%.txt', sep='\t')
  write.table(df_final_filtered_infants, file='V-infants_KT_filtered_sums.txt', sep='\t')

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
  write.table(perc2, file='V-mothers_KT_filtered_percentages%.txt', sep='\t')
  write.table(df_final_filtered_mothers, file='V-mothers_KT_filtered_sums.txt', sep='\t')

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
write.csv(freq_table, file="Virus_frequency_filtered_final_table.csv")



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

ggsave(filename ="V-KT_filt_species_frequency.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()






















################################################################################################
###############################PCA PLOT#################################
####################################################
###############################PCA PLOT################################################################
setwd("/HANGER/katie/raw_viromescan_bbmap_virus/final/deff_viral/PCA_plots/")

  .libPaths("/home/aaron/R/x86_64-pc-linux-gnu-library/3.6/")
library(ggcorrplot)
library(FactoMineR)
library(corrr)
library(factoextra)
head(final_filtered)
  metadata<- read.csv("/HANGER/katie/metadata_infant_mothers_meds.csv",header=TRUE)
metadata2<-metadata[match(colnames(final_filtered),metadata$sample_id),]
#metadata3<-metadata2[order(metadata2$PHASE),]
#metadata3<-metadata2[order(metadata2$PHASE),]
#col<-ifelse(metadata2$PHASE=="INFANT","red","blue")
#reordered_final<-final_filtered[,match(colnames(final_filtered),metadata3$sample_id)]

colnames(final_filtered)<-paste0(metadata2$PHASE,sep="_", metadata2$Siobhan_code)

va_sq<-scale(final_filtered)
corr_matrix<-cor(va_sq)




ggcorrplot(corr_matrix,hc.order=TRUE) +
ggplot2::theme(axis.text.x=element_text(size=5, angle=90, hjust =1),axis.text.y=element_text(size=5),panel.background = element_blank()

###colour=col if needed

hc.order=TRUE

ggsave(filename ="V-ordered_heat_map.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()


######repeat for human all virsues 





data.PCA<-princomp(corr_matrix)
summary(data.PCA)
fviz_eig(data.PCA, addlabels =TRUE)

ggsave(filename ="V-PCA_components.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()


fviz_pca_var(data.PCA, col.var ="black")
ggsave(filename ="V-PCA_plot.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()

fviz_cos2(data.PCA, choice ="var", axes =1:2)+
ggplot2::theme(axis.text.x=element_text(size=7, angle=90),
plot.background=element_blank()) 
ggsave(filename ="V-bar_plot.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()


fviz_pca_var(data.PCA, col.var = "cos2",
gradient.cols = c("black", "orange", "green"),
repel=TRUE)
ggsave(filename ="V-PCA_coloured_plot.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
       

save.image("final_filtere_virus.RData")



######alpha and beta diversity of bacteriophages #####
require("phyloseq")
require("ape")
OTU=otu_table(final_filtered,taxa_are_rows=TRUE)
TAX = tax_table(IDs)
colnames(TAX)<-c("Family","Genus","Species","ID Code")
 rownames(TAX)<-TAX[,4]

physeq=phyloseq(OTU,TAX)
library(ape) 
plot_bar(physeq,fill="Genus")

random_tree= rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

physeq1<-merge_phyloseq(physeq,metadata2, random_tree)
physeq2<-phyloseq(OTU,TAX,metadata2, random_tree)
identical(physeq2,physeq1)
plot_tree(physeq1, color="PHASE", label.tips="Species", ladderize="left", plot.margin=0.3)
plot_heatmap(physeq1)
