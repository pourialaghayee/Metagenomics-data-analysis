# set directory 
setwd("C:/Users/poria/OneDrive/Desktop/project/data/")

# library requirement
library(Biostrings)
library(ggplot2)
library(stringr)

#### read data####
data <- readDNAStringSet("simulated_dataset_uniform.fna")
ground_truth <- read.table("simulated_dataset_uniform.csv",
                           stringsAsFactors = FALSE, 
                           strip.white = TRUE,
                           sep = ' ')
cluster_data <- read.table("contigs.bin.table.tsv",
                           stringsAsFactors = FALSE, 
                           strip.white = TRUE,
                           sep = ' ')

####first part####
storage <- numeric(length = length(data))
for (i in 1:length(data)){
  storage[i] <- length(data[[i]])
}
# parameter in section 1
length(storage)
min(storage)
max(storage)
median(storage)
mean(storage)
sd(storage)
var(storage)

# hist 1
setwd("C:/Users/poria/OneDrive/Desktop/project/result/section 1/")
df <- data.frame(storage)
names(df) <- c('count')
pdf(file = paste("conting lenght diagram",".pdf",sep = "",collapse = NULL))
ggplot(df, aes(x=count)) + geom_histogram()+xlim(0,100000)+labs(y= "conting length", x = "conting #")
dev.off()

####second part####
for (variable in ground_truth) {
  x = (variable)
}
species <- numeric(length = length(data))
contig <- numeric(length = length(data))
for (i in 1:length(x)){
  tmp <- strsplit(x[i],"\t")
  contig[i] <- tmp[[1]][1]
  species[i] <- tmp[[1]][2]
}
df2 = data.frame(contig, species)
names(df2) <- c('ccontig', 'species')

# type  of species
length(unique(df2[["species"]]))

temp <- unique(df2[["species"]])
temp2 <- numeric(length = length(temp))
for (i in 1:length(temp)){
  temp2[i] <- length(which(grepl(temp[i], df2$species)))
}
df3 <- data.frame(temp2, temp)
names(df3) <- c('numberofconting', 'species')

setwd("C:/Users/poria/OneDrive/Desktop/project/result/section2/")
# plot number of contig of a species
pdf(file = paste("# contig of a species",".pdf",sep = "",collapse = NULL))
ggplot(data=df3, aes(x=numberofconting, y=species))+geom_bar(stat="identity")
dev.off()

uniq_species = unique(species)
#### A C T G for species ####
for (i in 1:length(uniq_species)){
  print(i)
  num_of_A <- 0
  num_of_C <- 0
  num_of_G <- 0
  num_of_T <- 0
  for (j in 1:length(contig)){
    if(uniq_species[i]==species[j]){
      temp <- match(contig[j],names(data))
      num_of_A = str_count(data[[temp]], "A") + num_of_A
      num_of_C = str_count(data[[temp]], "C") + num_of_C
      num_of_G = str_count(data[[temp]], "G") + num_of_G
      num_of_T = str_count(data[[temp]], "T") + num_of_T
    }
  }
  name = paste("spacies",i)
  name = paste(name,".pdf",sep='')
  pdf(name)
  plt = (c(num_of_A, num_of_C, num_of_G, num_of_T))
  barplot(plt, main="Num of nucleotides",
          xlab="A   C   G   T")
  dev.off()
}


#### 3 ####
uniq_cluster = sort(unique(cluster_data[[2]]))
num_of_each_cluster =  numeric(length = length(uniq_cluster))
zero <- 0
for (i in 1:length(cluster_data[[2]])){
  tmp = cluster_data[[2]][i]
  if (tmp==0){
    zero = zero +1
  }
  num_of_each_cluster[tmp] = num_of_each_cluster[tmp]+1
}
setwd("C:/Users/poria/OneDrive/Desktop/project/result/section3/")
pdf("number of each cluster.pdf")
barplot(num_of_each_cluster, main="number of each cluster",
        xlab="from 1 to 38")
dev.off()
# number of zerro cluster
print(zero)

#### four section ####
setwd("C:/Users/poria/OneDrive/Desktop/project/result/section4/")
str_f = "@Version:0.9.1\n@SampleID:gsa\n\n@@SEQUENCEID\tBINID\tTAXID\tLENGTH"
write(str_f,file="GR.binning",append=FALSE)

for (i in 1:length(contig)){
  line0=paste("c",contig[i],sep='')
  line1=paste(line0,species[i],sep="\t")
  line2=paste(species[i],storage[i],sep="\t")
  line=paste(line1,line2,sep="\t")
  
  write(line,file="GR.binning",append=TRUE)
}

str_f2 = "@Version:0.9.1\n@SampleID:gsa\n\n@@SEQUENCEID\tBINID"
write(str_f2,file="BB.binning",append=FALSE)

for (i in 1:length(contig)){
  line0=paste("c",cluster_data[[1]][i],sep='')
  line=paste(line0,cluster_data[[2]][i],sep='\t')
  
  write(line,file="BB.binning",append=TRUE)
}

setwd("C:/Users/poria/OneDrive/Desktop/project/result/section5/")
#### section 5 ####
num_of_A <- numeric(length = length(contig))
num_of_C <- numeric(length = length(contig))
num_of_G <- numeric(length = length(contig))
num_of_T <- numeric(length = length(contig))

for (i in 1:length(contig)){
  num_of_A[i] = str_count(data[[i]], "A")/storage[i]
  num_of_C[i] = str_count(data[[i]], "C")/storage[i]
  num_of_G[i] = str_count(data[[i]], "G")/storage[i]
  num_of_T[i] = str_count(data[[i]], "T")/storage[i]
}

df_me <- data.frame(num_of_A, num_of_C, num_of_G, num_of_T)
rownames(df_me) <- contig

x = kmeans(df_me,centers = 47)

center <-x$cluster

str_f3 = "@Version:0.9.1\n@SampleID:gsa\n\n@@SEQUENCEID\tBINID"
write(str_f3,file="ME.binning",append=FALSE)

for (i in 1:length(contig)){
  line0=paste("c",contig[i],sep='')
  line1=paste(line0,center[i],sep="\t")
  write(line1,file="ME.binning",append=TRUE)
}

