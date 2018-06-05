############################################################3
#############################################################
### de novo assemblies of environmental vibrio strains.
#############################################################


# Choice of parameters in velvet assembly.
library(plyr)

# Produce different kMERS velveth 

# Modify different exp_cov. parameters


# plotr
dat<-read.csv('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Summary_Assemblies_parameters_DT7.csv',h=T)
summary(dat)
dat$Kmers<-as.factor(dat$Kmers)
plot(dat$Kmers,dat$n50)
library(ggplot2)

ggplot(dat[dat$Strain=="DT7",], aes(x=exp_cov, y=n50, color=Kmers))+
  geom_bar(stat="identity", position=position_dodge())

ggplot(dat[dat$Strain=="372",], aes(x=exp_cov, y=n50, color=Kmers))+
  geom_bar(stat="identity", position=position_dodge())

plot(dat$max~dat$n50)
dat[dat$n50>900000,]




#############################################################
# Estimation by looking multiples assemblies.

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_Thailand/Choice_velvet_parameters/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.delim(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub("ALLK_","",gsub("OUT_","",gsub(".txt","",filesToProcess)    )    )

listOfFiles<-lapply(listOfFiles,function(x) as.character(x[,1]))


Assemblies<-lapply(listOfFiles,function(x) 
  unlist(strsplit(gsub("/"," ",     gsub(" reads","",    gsub(", using","",      gsub(", total","",   gsub(", max","",
    gsub("nodes and n50 of ","",gsub("Final graph has ","",x[grep("Final graph",x)]   )    )  )   )    )    )   )
                  , split = " " )  )
  )

Assemblies2<-list()
for ( i in 1:length(Assemblies)){
Assemblies2[[i]]<-c(unlist(strsplit(names(Assemblies[i]),split="_")),Assemblies[[i]])
tryCatch(names(Assemblies2[[i]])<-c("Strain","K","Exp_cov","Nodes","N50","Max","Total","UsingReads","TotalReads"),
                                  error= function (e) names(Assemblies2[[i]])<-c("Strain","K","Exp_cov"))

}

n.obs <- sapply(Assemblies2, length)
seq.max <- seq_len(max(n.obs))
Assemblies3 <- t(sapply(Assemblies2, "[", i = seq.max))

Assemblies3 <- as.data.frame(Assemblies3)

Assemblies3$Exp_cov<-   factor(Assemblies3$Exp_cov, levels = c("10","20","30","40","50","60","70","80","90",
                                                 "100","110","120","130","140","150","160","170","180","190","200"))
Assemblies3$Nodes<-as.numeric(levels(Assemblies3$Nodes))[Assemblies3$Nodes]
Assemblies3$N50<-as.numeric(levels(Assemblies3$N50))[Assemblies3$N50]
Assemblies3$Max<-as.numeric(levels(Assemblies3$Max))[Assemblies3$Max]
Assemblies3$Total<-as.numeric(levels(Assemblies3$Total))[Assemblies3$Total]
Assemblies3$UsingReads<-as.numeric(levels(Assemblies3$UsingReads))[Assemblies3$UsingReads]
Assemblies3$TotalReads<-as.numeric(levels(Assemblies3$TotalReads))[Assemblies3$TotalReads]



ggplot(Assemblies3[Assemblies3$Strain=="354",], aes(x=Exp_cov, y=N50, color=K))+
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                                               axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  )
Assemblies3[Assemblies3$Strain=="2006",]
Assemblies3[Assemblies3$Strain=="354" & Assemblies3$N50>500000,]

Assemblies3[Assemblies3$Strain=="WKB-T9" & Assemblies3$Exp_cov==200,]



list_plots<-list()
for (X in unique(Assemblies3$Strain)){
 list_plots[[X]]<-ggplot(Assemblies3[Assemblies3$Strain==X,], aes(x=Exp_cov, y=N50, color=K)) + ggtitle(X) +
    geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                                                 axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
    )
}

### PLOTING ALL

pdf("ALL_Estimations_Thailand.pdf",height = 8,width = 10)
for (i in 1:length(unique(Assemblies3$Strain))) {
  print(list_plots[[i]])
}
dev.off()
