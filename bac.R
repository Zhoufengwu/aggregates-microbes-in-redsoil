bacterial<- function(){setwd("./16s")

map = read.csv("group.csv",row.names = 1,header = T)
map=unite(map, "yeartre", year,tre, sep = ".", remove = FALSE)
asv = read.table("otu_table.tsv",header = T,row.names = 1)
asv=asv[,colnames(asv)%in%row.names(map)]
c=asv
tree =read.tree("tree.nwk")
tax = read.delim("taxonomy.tsv",row.names = 1,header = T)
tax=subset(tax,select = -c(Confidence))
b=str_split_fixed(tax$Taxon, ";", 7)
tax=cbind(tax,b)
tax=subset(tax,select = -c(Taxon))
names(tax)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
# remove the sequences we don't need, such as Unassigned sequences, here we will drop out few information.
# drop out the sequences Unassigned on phylum level. We don't use the ASVs which unassigned on phylum level, cause there are no ecologic means with those sequences. 
# selecting Archaea
taxarch<- tax[which(tax$Kingdom == "d__Archaea"), ]
taxarch<- taxarch[!(taxarch$Phylum == ""), ]
write.csv(taxarch,"taxarch.csv")
taxarch$row_ID <- row.names(taxarch)
asv$row_ID <- row.names(asv)
new=semi_join(asv,taxarch, by = "row_ID")
asvarch= subset(new,select = -c(row_ID))
write.csv(asvarch,"asvarch.csv")


tax<- tax[!(tax$Kingdom == "d__Archaea"), ]
tax<- tax[!(tax$Phylum == ""), ]
tax$row_ID <- row.names(tax)
asv$row_ID <- row.names(asv)
new16s=semi_join(asv,tax, by = "row_ID")
asv= subset(new16s,select = -c(row_ID))
print("number of drop out 16S ASVs",str(nrow(c)-nrow(asv)))
print("the min count of sequences",str(colSums(asv)[which.min(colSums(asv))]))
# Using function rrarefy of vegan to subsample asv-table
asv_sub16s = as.data.frame(t(rrarefy(t(asv), min(colSums(asv)))))
print(colSums(asv_sub16s))    
# note: raring is a method to estimate alpha diversity.
# NOW, U MIGHT HAVE A GOOD ASV TABLE TO ANALYSIS MCDF.
# Next, we will change the counts into percentages, usually we analysis relative abundance more than 0.1%,0.5% and 1%.
# This time we also show U a method published on FEMS Microbiology Ecology, doi: 10.1093/femsec/fiw203
asv_abundacne = asv_sub16s/colSums(asv_sub16s)
n=ncol(asv_abundacne)
asv_abundacne0.01=asv_abundacne[which(rowSums(asv_abundacne)/n>0.001),]
asv_abundacne0.05=asv_abundacne[which(rowSums(asv_abundacne)/n>0.005),]
asv_abundacne0.1=asv_abundacne[which(rowSums(asv_abundacne)/n>0.01),]
otu16s=list(asv_sub16s,map,tax,asv_abundacne,asv_abundacne0.01,asv_abundacne0.05,asv_abundacne0.1)

}

