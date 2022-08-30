#input ASV table
#Constructing network. you can choose the optimal method which you like
#but we still suggest to use SPARCC bulid your co-networks
#install packages
#install.packages('devtools')
#devtools::install_github('zdk123/SpiecEasi')
#install.packages('igraph')
#install.packages('randomForest')
library('randomForest')
library('igraph')
library('SpiecEasi')
dir.create('test', recursive = TRUE)
otu=read.csv("test/test.csv",header = T,row.names = 1)
map=read.csv("test/group.csv",header = T,row.names = 1)
#Extract a subset which ASVs' relative abundance greater than 0.05%
otud<-otu[which(rowSums(otu/colSums(otu))/ncol(otu) >= 0.0005),]

  otu_n<- cbind(t(otud), map[c("tre")])
  obb=NULL
  
  for (i in 5:nrow(otu_n)) {
    #i=14
    set.seed(123)
    train_n <- sample(nrow(otu_n), nrow(otu_n)*(i/15))
    otu_train_n <- otu_n[train_n, ]
    otu_test_n <- otu_n[-train_n, ]
    set.seed(123)
    otu_train.forest <- randomForest(as.factor(tre)~., data = otu_train_n, importance = TRUE,ntree=50000)
    #summary(otu_train.forest)
    otu_train.forest
    x=data.frame(otu_train.forest[["err.rate"]])$OOB
    z=data.frame(table(x))
    obbnum=as.vector(z$x)[which.max(z$Freq)]
    oob_i=c(i,as.numeric(obbnum)*100)
    obb=rbind(obb,oob_i)}
  print(obb[,1][which.min(obb[,2])])
  #Possible parameter
  pp=obb[,1][which.min(obb[,2])]
  set.seed(123)
  train_n <- sample(nrow(otu_n), nrow(otu_n)*(pp/15))
  otu_train_n <- otu_n[train_n, ]
  otu_test_n <- otu_n[-train_n, ]

  set.seed(123)
  otu_train.forest <- randomForest(as.factor(tre)~., data = otu_train_n, importance = TRUE,ntree=50000)
  otu_train.forest
  x1=data.frame(otu_train.forest[["err.rate"]])$OOB
  z1=data.frame(table(x1))
  obbnum1=as.vector(z1$x)[which.max(z1$Freq)]
  oob_i1=c(i,as.numeric(obbnum1))
  oob1=cbind("OOB estimate of  error rate",oob_i1[2])
  set.seed(123)
  otu_train.cv <- replicate(10, rfcv(otu_n[-ncol(otu_n)], as.factor(otu_n$tre), cv.fold = 10,step = 1.5), simplify = FALSE)
  otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
  otu_train.cv$otus <- rownames(otu_train.cv)
  otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
  otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
  otu_train.cv<-otu_train.cv[order(otu_train.cv$otus),]
  test=aggregate(subset(otu_train.cv,select=-c(variable)),by = list(otu_train.cv$otus),FUN =mean) 
  
  library(ggplot2)
  library(splines)
  
  p <- ggplot(otu_train.cv, aes(otus, value)) +
    geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
  p=p + geom_vline(xintercept = 50)
  print(p)
  
  print(varImpPlot(otu_train.forest, n.var = min(50, nrow(otu_train.forest$importance)), 
                   main = paste('Top',50, '- variable importance',sep='') ))
  importance_otu <- data.frame(importance(otu_train.forest))
  importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
  #importance_otu<-cbind(importance_otu,taxonomy[row.names(importance_otu),]);if you hava taxonomy table
  
  otu_select <- rownames(importance_otu)[1:50]
  otu_train_top <- otu_n[ ,c(otu_select, 'tre')]
  write.csv(rbind(subset(test,select = -c(Group.1)),oob1[1,]),paste('test/RFoob',".csv",sep=''))
  write.csv(otu_train_top,paste('test/','biomaker',".csv",sep=''))
  write.csv(importance_otu[1:50,], paste('test/','top','improtance',".csv",sep=''))


#network
NET=function(otu){
  #otu=otud
set.seed(123)
#a slow step
amgut1.filt.sparcc <- sparcc(t(otu), iter = 20, inner_iter = 10, th = 0.1)

sparcc0 <- amgut1.filt.sparcc$Cor

colnames(sparcc0) <- colnames(t(otu))
rownames(sparcc0) <- colnames(t(otu))
#output the correlation matrix
write.table(sparcc0, 'test/sparcc0.txt', sep = '\t', 
col.names = NA, quote = FALSE)

#Get random cor-matrix by 100 bootstrap samples
#slow......
set.seed(123)
n = 100

for (i in 1:n) {
  amgut1.filt.boot <- sample(otu, replace = TRUE)  #bootstrap
  amgut1.filt.sparcc_boot <- sparcc(amgut1.filt.boot, iter = 20, 
                                    inner_iter = 10, th = 0.1)  
  sparcc_boot <- amgut1.filt.sparcc_boot$Cor
  colnames(sparcc_boot) <- colnames(amgut1.filt.boot)
  rownames(sparcc_boot) <- colnames(amgut1.filt.boot)
  write.table(sparcc_boot, paste('test/testsparcc_boot', i, '.txt', sep = ''), 
              sep = '\t', col.names = NA, quote = FALSE)  
}
# p-value
p <- sparcc0
p[p!=0] <- 0

for (i in 1:n) {
  p_boot <- read.delim(paste('test/testsparcc_boot', i, '.txt', sep = ''), sep = '\t', 
                       row.names = 1)
  p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
}

p <- p / n
write.table(p, 'test/pvals.two_sided.txt', sep = '\t', col.names = NA, quote = FALSE)

#network
cor_sparcc <- read.delim('test/sparcc0.txt', row.names = 1, sep = '\t', check.names = FALSE)
pvals <- read.delim('test/pvals.two_sided.txt', row.names = 1, sep = '\t', check.names = FALSE)
cor_sparcc[abs(cor_sparcc)<=0.3 | pvals>=0.05] <- 0
diag(cor_sparcc) <- 0
write.table(cor_sparcc, 'test/neetwork.adj.txt', col.names = NA, sep = '\t', quote = FALSE)
}

NET(otud)
library(igraph)
adjacency_unweight <- read.delim('test/neetwork.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(adjacency_unweight)[1:6] 
g <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = TRUE, diag = FALSE)

adjacency_unweight <-abs(adjacency_unweight)

adjacency_unweight[adjacency_unweight >0] = 1

g <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = TRUE, diag = FALSE)


V(g)$degree <- igraph::degree(g,normalized = T)    
nodes_list <- data.frame(
  nodes_id = V(g)$name, 
  degree = V(g)$degree
 
)
head(nodes_list)   
write.table(nodes_list,paste('test/nodepro','.txt',sep=''), sep = '\t', row.names = F, quote = FALSE)


nodes_list <- read.delim(paste('test/nodepro','.txt',sep=''),row.names = 1, sep = '\t', check.names = FALSE)
nam=row.names(nodes_list)[order(nodes_list$degree, decreasing = F)]
nodes_list  <- data.frame(nodes_list[order(nodes_list$degree, decreasing = F),])
row.names(nodes_list)=nam
names(nodes_list)=c('ND')
nodes_list$no<-order(-nodes_list[,1])
#Targeted removal
#removal with replacement
clustering_coefficient=NULL
for (i in 1:length(row.names(nodes_list))) {
  
  adjacency_unweight_i=adjacency_unweight[row.names(nodes_list),row.names(nodes_list)]
  adjacency_unweight_i=adjacency_unweight_i[-i,-i]
  g_i <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight_i), mode = 'undirected', weighted = TRUE, diag = FALSE)
  clustering_coefficient_i <- transitivity(g_i)
  clustering_coefficient=rbind(clustering_coefficient,clustering_coefficient_i)
}
clustering_coefficient=data.frame(clustering_coefficient)

nodes_list=cbind(nodes_list,clustering_coefficient)
library(ggplot2)

p1=ggplot(nodes_list, aes(x=ND, y=clustering_coefficient)) + 
  geom_line()+
  geom_point(size=0.2)+
  xlab("ND of ASVs removed")+
  ylab("clustering coefficient of species remained")+
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6))+
  geom_vline(xintercept = 0.125)
p1
#removal with no-replacement
nodes_list2 <- read.delim(paste('test/nodepro','.txt',sep=''),row.names = 1, sep = '\t', check.names = FALSE)
nam2=row.names(nodes_list2)[order(nodes_list2$degree, decreasing = F)]
nodes_list2  <- data.frame(nodes_list2[order(nodes_list2$degree, decreasing = F),])
row.names(nodes_list2)=nam
names(nodes_list2)=c('ND')
nodes_list2$no<-order(-nodes_list2[,1])
clustering_coefficient2=NULL
for (i in 1:(length(row.names(nodes_list2))-30)) {
  
  adjacency_unweight_i=adjacency_unweight[row.names(nodes_list2),row.names(nodes_list2)]
  adjacency_unweight_i=adjacency_unweight_i[-(1:i),-(1:i)]
  g_i <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight_i), mode = 'undirected', weighted = TRUE, diag = FALSE)
  clustering_coefficient_i <- transitivity(g_i)
  clustering_coefficient2=rbind(clustering_coefficient2,clustering_coefficient_i)
}
nodes_list2=cbind(nodes_list2[1:(nrow(nodes_list2)-30),],clustering_coefficient2[,1])
names(nodes_list2)=c('ND','NO.','CC')
p2 =ggplot(nodes_list2, aes(x=ND, y=CC)) + 
  geom_line()+
  geom_point(size=0.2)+
  xlab("ND of ASVs removed")+
  ylab("clustering coefficient of species remained")+
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6))+
  geom_vline(xintercept = 0.125)
p2
nodes_list[which(nodes_list$ND>0.125),'roles'] <- 'coremicrobiome'

write.csv(nodes_list, paste('test/cormi','.csv',sep=''), quote = FALSE)

#overlap
rfotu=read.csv(paste('test/','biomaker',".csv",sep=''),header = T,row.names = 1)
rfotu=subset(rfotu,select = -c(tre))
rfotu=as.data.frame(t(rfotu))
netotu=read.csv(paste('test/cormi','.csv',sep=''),header = T)
netotu=netotu[netotu$roles=="coremicrobiome",]
reim=rfotu[intersect(row.names(rfotu),netotu$X),]
write.csv(reim,"keystonetaxa.csv")
#done



