
tar -xzvf /HANGER/katie/kraken/kraken_database/k2_viral_20230605.tar.gz

kraken2 --db /HANGER/katie/kraken/kraken_database/ --threads 8 --output/HANGER/katie/kraken/AS007_kraken2.out --report /HANGER/katie/kraken2/AS_007_kraken2.report --paired  /HANGER/katie/sequences/fastq/AS007_1.fq.gz /HANGER/katie/sequences/fastq/AS007_2.fq.gz

##full for loop including all 166 samples  
for i in AS206 AS217 AS118 AS219 AS198 AS175 AS154 AS167 AS123 AS145 AS065 AS134 AS141 AS120 AS069 AS121 AS073 AS143 AS007 AS013 AS033 AS026 AS038 AS040 AS060 AS053 AS064 AS047 AS025 AS066 AS057 AS194 AS225 AS227 AS179 AS041 AS048 AS133 AS229 AS196 AS202 AS207 AS210 AS212 AS222 AS200 AS184 AS170 AS204 AS218 AS191 AS160 AS148 AS186 AS181 AS149 AS176 AS187 AS166 AS161 AS122 AS137 AS156 AS157 AS162 AS117 AS135 AS151 AS124 AS074 AS114 AS129 AS076 AS021 AS034 AS127 AS035 AS030 AS027 AS029 AS014 AS062 AS059 AS205 AS224 AS125 AS220 AS197 AS159 AS119 AS153 AS140 AS136 AS115 AS113 AS146 AS150 AS138 AS142 AS071 AS139 AS011 AS016 AS020 AS042 AS032 AS068 AS063 AS055 AS046 AS043 AS052 AS049 AS045 AS193 AS226 AS221 AS182 AS058 AS112 AS132 AS228 AS195 AS201 AS208 AS211 AS209 AS223 AS199 AS185 AS171 AS203 AS214 AS192 AS152 AS158 AS183 AS180 AS168 AS164 AS188 AS174 AS163 AS177 AS126 AS131 AS169 AS165 AS116 AS061 AS144 AS173 AS056 AS051 AS130 AS050 AS028 AS015 AS147 AS036 AS012 AS037 AS024 AS031 AS022 AS044; 
do

kraken2 --db /HANGER/katie/kraken/kraken_database/ kraken2 --paired --classified-out AS_007_#.fq AS007_1.fq.gz AS007_2.fq.gz



kraken2- --db <path_to_kraken2_database> threads --8 --classified-out <output_classified_sequences_file> --unclassified-out <output_unclassified_sequences_file> <input_kraken2_output_file>

do kraken2 --db /HANGER/katie/kraken/kraken_database/ --threads 8 --classified-out /HANGER/katie/kraken/classified/"$i"#_kraken2.out --report /HANGER/katie/kraken/"$i"_kraken2.report --use-names --gzip-compressed --paired  /HANGER/katie/sequences/fastq/"$i"_1.fq.gz /HANGER/katie/sequences/fastq/"$i"_2.fq.gz
--report /HANGER/katie/kraken2/AS_007_kraken2.report


kraken2 --db /HANGER/katie/kraken/kraken_database/ --threads 8 --classified-out /HANGER/katie/kraken/classified/AS217#_kraken2.out --report /HANGER/katie/kraken/AS217_kraken2.report --use-names --gzip-compressed --paired  /HANGER/katie/sequences/fastq/AS217_1.fq.gz /HANGER/katie/sequences/fastq/AS217_2.fq.gz


## for loop not including the 12 sequences that have very low read depth and cannpt be used. 


for i in AS206 AS217 AS118 AS219 AS198 AS175 AS154 AS167 AS123 AS145 AS065 AS134 AS141 AS120 AS069 AS121 AS073 AS143 AS007 AS013 AS038 AS040 AS060 AS053 AS064 AS047 AS066 AS057 AS194 AS225 AS227 AS179 AS041 AS048 AS133 AS229 AS196 AS202 AS207 AS210 AS212 AS222 AS200 AS184 AS170 AS204 AS218 AS191 AS160 AS186 AS181 AS176 AS187 AS166 AS161 AS122 AS137 AS156 AS157 AS162 AS117 AS135 AS151 AS124 AS074 AS114 AS129 AS076 AS021 AS127 AS027 AS014 AS062 AS059 AS205 AS224 AS125 AS220 AS197 AS159 AS119 AS153 AS140 AS136 AS115 AS113 AS146 AS150 AS138 AS142 AS071 AS139 AS011 AS016 AS020 AS042 AS032 AS068 AS063 AS055 AS046 AS043 AS052 AS049 AS045 AS193 AS226 AS221 AS182 AS058 AS112 AS132 AS228 AS195 AS201 AS208 AS211 AS209 AS223 AS199 AS185 AS171 AS203 AS214 AS192 AS152 AS158 AS183 AS180 AS168 AS164 AS188 AS174 AS163 AS177 AS126 AS131 AS169 AS165 AS116 AS061 AS144 AS173 AS056 AS051 AS130 AS050 AS015 AS147 AS012 AS037 AS024 AS022 AS044 ;
do kraken2 --db /HANGER/katie/kraken/kraken_database/ --threads 8 --classified-out /HANGER/katie/kraken/classified/"$i"#_kraken2.out --report /HANGER/katie/kraken/"$i"_kraken2.report --use-names --gzip-compressed --paired  /HANGER/katie/sequences/fastq/"$i"_1.fq.gz /HANGER/katie/sequences/fastq/"$i"_2.fq.gz
echo "finished_$i" ; 
done




Aaron made a new kraken2 conda environment and install Kraken2 /home/katie/miniconda3/envs/Kraken2

kraken2-build --download-library viral --db /HANGER/katie/kraken/kraken_database


kraken2-build --download-taxonomy --db /HANGER/katie/kraken/kraken_database

 kraken2-build --build --db /HANGER/katie/kraken/kraken_database --threads 16

## to show what is in the database 

kraken2-inspect --db /HANGER/katie/kraken/kraken_database/ --display


###summarise kraken data output #threshold 10 reads (anything mapped less than 10 times will not show up) using bowtie2 through viromescan pipeline attributing a higher read count to reach virus and therefore passing our threshold (10) and showing up in results. Using kracken2, viruses did not surpass threshold therefore do not show in results. kraken has a more stringent threshold for assigning reads to viral sequences. 

cd /HANGER/katie/kraken/

for f in *report; do echo "awk '\$1>0 && \$3>10' "$f" |uniq|sort -k1,1rn  >"$f".summary" ;done >summariser.sh
sh summariser.sh

##code to move all summary files into one folder

find /HANGER/katie/kraken -type f -name "*_kraken2.report.summary" -print0 | xargs -0 -I cheese mv -i cheese /HANGER/katie/kraken/summary


##plotting reults in R 
                                                                                                                                                                                        

setwd("/HANGER/katie/kraken/results")
.libPaths("/home/aaron/R/x86_64-pc-linux-gnu-library/3.6/")
library(knitr)
library(dplyr)
library(ggplot2)

save.image("kraken_data.RData")

file_listB <- list.files(path = "/HANGER/katie/kraken/summary",
                         pattern="*_kraken2.report.summary$",
                         recursive = TRUE,
                         full.names = TRUE)
head(file_listB)
IDs <- as.data.frame(read.csv("/HANGER/katie/kraken/kraken_database/ktaxonomy.tsv", sep="\t", header=F))
colnames(IDs)[1]<-c("tax_code")
colnames(IDs)[9]<-c("name")
IDs<-IDs[,c(1,5,9)]


for (i in  1:length(file_listB)){
  oneFile <- file_listB[i]
  oneFileData <- as.data.frame(read.delim(oneFile, sep="\t", header=F, skip=1))
  
  #extract virus name and aligned columns
  oneFileDataSmall <- oneFileData[,c(5,4,2)]
  oneFileDataSmall<-oneFileDataSmall[oneFileDataSmall$V4=="G",]
  oneFileDataSmall$V4<-NULL
  oneFileName <- unlist(strsplit(oneFile, "/"))
  #get only file name - that last item
  oneFileName <- oneFileName[length(oneFileName)]
  oneFileName <- unlist(strsplit(oneFileName, "_"))
  oneFileName<-oneFileName[1]
  
  
  colnames(oneFileDataSmall) <- c("tax_code",oneFileName)
                            
  
  IDs<-left_join(IDs,oneFileDataSmall, by="tax_code")
  


}
df<-IDs[3:nrow(IDs),]
df[is.na(df)] <- 0
rownames(df)<-paste0(df[,1],sep="_",df[,3])
df[,c(1,2,3)]<-NULL

isnum <- sapply(df, is.numeric)
df[isnum] <- lapply(df[isnum], function(x) ifelse(x == 0, as.numeric(0), as.numeric(1)))

df$sums<-rowSums(df)
df<-df[order(df$sums,decreasing = TRUE),]
head(df)

write.csv(IDs, "taxonomy_aundance_kraken.csv")

cut<-df$sums=="0"
table(cut)
###101  kept and 25728 cut out
df<-df[!cut,]

df$sums<-NULL

write.csv(df, "kraken_frequency.csv", row.names =T)

###seperating mothers 

  metadata<- read.csv("/HANGER/katie/metadata_infant_mothers_meds.csv",header=TRUE)
  df<-read.csv("kraken_frequency.csv")
  rownames(df)<-df$X
  df$X<-NULL
kraken_metadata<-metadata[match(colnames(df),metadata$sample_id),]

infants_metadata<- kraken_metadata[kraken_metadata$PHASE=="INFANT",]
infant_codes<- infants_metadata$sample_id
mothers_metadata<-kraken_metadata[kraken_metadata$PHASE=="MOTHER",]
mother_codes<-mothers_metadata$sample_id


mothers_viruses<-df[,is.na(match(colnames(df),infant_codes))]

match(mothers_metadata$sample_id,colnames(mothers_viruses))
## seperating infants

infant_viruses<-df[,is.na(match(colnames(df),mother_codes))]

infant_viruses$sums<-rowSums(infant_viruses)
mothers_viruses$sums<-rowSums(mothers_viruses)




##########################calculating sample frequency and graphing################
perc2<-NULL
perc<-NULL

infant_viruses<-as.data.frame(infant_viruses)
 perc=c()
  for (i in 1 : (length(rownames(infant_viruses))))
  {
    perc[i]=round(infant_viruses$sums[i]/sum(infant_viruses$sums),4)
  }
 perc2=as.matrix(perc)
  rownames(perc2)=rownames(infant_viruses)
  colnames(perc2)=c('')
  ff=(perc2>0)
  perc2=perc2[ff,]
  perc2=as.matrix(perc2)
  fff=infant_viruses$sums>0
  infant_viruses=infant_viruses[fff,]
  infant_viruses=as.matrix(infant_viruses)
  colnames(perc2)="Detected_frequency"
  write.table(perc2, file='infants_kraken2_percentages%.txt', sep='\t')
  write.table(infant_viruses, file='infants_kraken2_sums.txt', sep='\t')

  filt=(perc2)<0.005
  perc2=as.data.frame(perc2[!filt,])
  table(filt)
  ##33 low abundance were filtered out, 52 were kept###
  
  infant_frequency<-as.data.frame(perc2)
  colnames(infant_frequency)<-"Detected_frequency"
  low_abundance=as.data.frame(round(1-(sum(infant_frequency)),2))
  colnames(low_abundance)<-"Detected_frequency"
  rownames(low_abundance)<-"low_abundance"
  
   infant_frequency=rbind(infant_frequency, low_abundance)
  infant_frequency$Virus<-rownames(infant_frequency)
infant_frequency$Group<-"Infants"
colnames(infant_frequency)<-c("Detected_Frequency","Virus","Group")
rownames(infant_frequency)<-NULL


##same as above for mothers
perc2<-NULL
perc<-NULL

mothers_viruses<-as.data.frame(mothers_viruses)
 perc=c()
  for (i in 1 : (length(rownames(mothers_viruses))))
  {
    perc[i]=round(mothers_viruses$sums[i]/sum(mothers_viruses$sums),4)
  }
 perc2=as.matrix(perc)
  rownames(perc2)=rownames(mothers_viruses)
  colnames(perc2)=c('')
  ff=(perc2>0)
  perc2=perc2[ff,]
  perc2=as.matrix(perc2)
  fff=mothers_viruses$sums>0
  mothers_viruses=mothers_viruses[fff,]
  mothers_viruses=as.matrix(mothers_viruses)
  colnames(perc2)="Detected_frequency"
  write.table(perc2, file='mothers_kraken2_percentages%.txt', sep='\t')
  write.table(mothers_viruses, file='mothers_kraken2_sums.txt', sep='\t')

  filt=(perc2)<0.005
  perc2=as.data.frame(perc2[!filt,])
  table(filt)
  ##176 low abundance were filtered out, 31 were kept###
  
  mothers_frequency<-as.data.frame(perc2)
  colnames(mothers_frequency)<-"Detected_frequency"
  low_abundance=as.data.frame(round(1-(sum(mothers_frequency)),2))
  colnames(low_abundance)<-"Detected_frequency"
  rownames(low_abundance)<-"low_abundance"
  
   mothers_frequency=rbind(mothers_frequency, low_abundance)
  mothers_frequency$Virus<-rownames(mothers_frequency)
mothers_frequency$Group<-"Mothers"
colnames(mothers_frequency)<-c("Detected_Frequency","Virus","Group")
rownames(mothers_frequency)<-NULL




freq_table<-rbind(mothers_frequency, infant_frequency)
###to read back in and plot it@@@
####freq_table<-read.csv("detection_frequency_table_human_virsues.csv")

freq_table$Genus <-unlist(sapply(strsplit(freq_table$Virus,"_"),tail,1))
write.csv(freq_table,file="detection_frequency _kraken.csv")


ggplot(freq_table,aes( fill=Genus,x=Group,y=Detected_Frequency)) + geom_bar(position="stack", stat="identity", colour="white") + xlab(paste("Group")) +
  ylab(paste("Abundance")) + theme(axis.line = element_line(color = 'black'), 
                        panel.background = element_blank() ,
                        panel.grid.major = element_blank() , 
                        panel.grid.minor = element_blank(), 
                        panel.border = element_blank(), 
                        axis.text.x = element_text(color = 'black', size=15), 
                        axis.text.y = element_text(color = 'black', size=15),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())
ggsave(filename ="Kraken_genus_freq.jpeg",dpi = 700, scale = 1,
       limitsize = TRUE)
graphics.off()



###code to match mothers to infants and find evidence of veritcal transmission (link dots by lines)

 total$pair<-metadata2[match(rownames(total),rownames(metadata2)),"Siobhan_code"]


total$duplicate=duplicated(total$pair) | duplicated(total$pair,fromLast=TRUE)
 total_paired<-total[total$duplicate =="TRUE",]
 total_paired
