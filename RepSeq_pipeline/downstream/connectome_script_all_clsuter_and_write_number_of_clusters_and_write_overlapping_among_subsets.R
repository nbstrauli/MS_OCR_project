#!/usr/bin/ Rscript

library(stringr)
library(reshape)
library(plyr)
library(dplyr)
library(combinat)

ptm <- proc.time()

source("init_shared.R")

# register multicore for foreach
library(doMC)
registerDoMC(cores=12)


data<-read.csv(all_clusters_csv_file,stringsAsFactors = FALSE)
#remove first column, 
# TODO trt data_without_rowname<-read.csv(all_clusters_csv_file, stringsAsFactors = FALSE, row.names=FALSE)
# NOTE: above does not work because read.csv does not expect row.name config
data<-data[,-1]
subset_list<-unique(data$subset)
if (length(subset_list) == 1) {
    print("only one subset")
}
####### PART I ########
# get number of clusters by susbset

#DONE: better ways to create a dataframe?
#count_cluster_per_subset<-as.data.frame(matrix(nrow=length(subset_list),ncol=2))
#colnames(count_cluster_per_subset)<-c('subset','number_of_clusters')

##DONE: use group function to do this --Hao
## select count(1) from data group by cluster_membership
#for (i in 1:length(subset_list)){
#  data_sub<-data[data$subset==subset_list[i],]
#  count_cluster_per_subset[i,]<-c(subset_list[i],length(unique(data_sub$cluster_membership)))
#}


# better version for the code segment above
count_cluster_per_subset <- data %>% group_by(subset) %>% 
  summarise(number_of_clusters = n_distinct(cluster_membership)) %>% collect()


print("number of clusters by subset")
print(count_cluster_per_subset)


count_cluster_per_subset<-count_cluster_per_subset[order(count_cluster_per_subset$subset),]
write.csv(count_cluster_per_subset,file.path(result_dir,paste0(patient, "_count_cluster_per_subset.csv")))



#only keep cluster with node from at least 2 subset
cluster<-unique(data$cluster_membership)
final_clus<-data.frame(subset1=as.character(),subset2=as.character())
for (data_clu in split(data, data$cluster_membership)) {
  subset_count_in_cluster_member<-as.character(unique(data_clu$subset))
  if (length(subset_count_in_cluster_member)>1){
    sub<-t(data.frame(combn(subset_count_in_cluster_member, 2,simplify = T)))
    final_clus<-rbind(final_clus,sub)
  }

}


#  TODO: use the parallel version
# filter_clusters_with_same_membership <- function(data_in_cluster) {
#   subset_count_in_cluster_member<-as.character(unique(data_in_cluster$subset))
#   if (length(subset_count_in_cluster_member)>1){
#     result_sub<-t(data.frame(combn(subset_count_in_cluster_member, 2,simplify = T)))
#     return(result_sub)
#   }
# }
# 
# final_clus<-foreach(data_in_subset = split(data,data$cluster_membership), .combine='rbind') %dopar% {
#    filter_clusters_with_same_membership(data_in_subset)
# }
#   
# colnames(final_clus) <-c('V1', 'V2')
# final_clus <- as.data.frame(final_clus)
  
# rewrite this part,. too slow the logic is bad, should be done in igraph



#final_clus
# count_of_subset_per_cluster <- data %>% group_by(cluster_membership)  %>% 
#   summarise(number_of_subset = n_distinct(subset)) %>%
#   filter( number_of_subset > 1)
# final_clus <- 
#   
  
write.table(final_clus,file.path(result_dir,paste0(patient, '_connectivity_rawlist')),sep='\t',row.names=F)
#print(final_clus)
if (nrow(final_clus) == 0) {
    stop("no connectivity exists")

} 

data<-final_clus[order(final_clus$V1,final_clus$V2),]

data2<-data.frame(table(data))
data2<-data2[data2$Freq!=0,]
data2$V1<-as.character(data2$V1)
data2$V2<-as.character(data2$V2)

data2$member<-ifelse(data2$V1<data2$V2,paste(data2$V1,data2$V2,sep=';'),paste(data2$V2,data2$V1,sep=';'))
conn<-ddply(data2[,c(3,4)], .(member), summarize,count = sum(Freq))

write.csv(data2,file.path(result_dir,paste0(patient, '_connectivity_rawlist_summary.csv')),row.names=F)
conn$source<-str_replace(str_extract(conn$member,'^.*;'),';','')
conn$target<-str_replace(str_extract(conn$member,';.*$'),';','')


conn<-conn[,-1]
write.table(conn,file.path(result_dir,paste0(patient, '_subset_connectome.txt')),
            row.names=F,col.names=F,quote =F)

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
