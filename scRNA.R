library(readr)
library(Seurat)
library( tidyverse)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(ggsignif)
library(CellChat)
library(patchwork)
#options(digits=5)


DEG2RNK <- function(DEG,p_hold = 1.1, log2fc_hold = 0,name){
  
  
  up <- DEG[which(DEG$avg_log2FC > 0 & abs(DEG$avg_log2FC) > log2fc_hold & DEG$p_val_adj < p_hold) ,]
  
  down <- DEG[which(DEG$avg_log2FC < 0 & abs(DEG$avg_log2FC) > log2fc_hold & DEG$p_val_adj < p_hold) ,]
  
  
  up <- up[order(up$avg_log2FC,decreasing = T),]
  
  down <- down[order(down$avg_log2FC,decreasing = F),]
  
  
  rnk <- DEG[order(DEG$avg_log2FC,decreasing = T),]
  write_tsv(data.frame(row=rownames(rnk),log2fc=rnk[,'avg_log2FC']) ,paste0(name,'.rnk'),col_names = F)
}

initialization <- function(){
  dt <- read.csv('GSE115978_counts.csv',row.names = 1)
  meta <- read.csv('GSE115978_cell.annotations.csv',row.names = 1)
  all <- CreateSeuratObject(dt,
                            project='GSE115978',
                            meta.data = meta)
  write_rds(all,'GSE115978.rds')
  return(all)
}

standard_process <- function(object){
  
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  
  object <- ScaleData(object)
  
  object <- RunPCA(object, features = VariableFeatures(object = object))
  
  
  object <- FindNeighbors(object, dims = 1:10)
  object <- FindClusters(object, resolution = 0.75) 
  object <- RunUMAP(object, dims = 1:10)
  
  
  return(object)
}

cellsubset_discover_DEG <- function(object,cell,name,high_sample=high_samples,low_sample=low_samples){
  
  ident.1 = colnames(object)[which(object$samples %in% high_sample & object$cell.types == cell)]
  
  ident.2 = colnames(object)[which(object$samples %in% low_sample & object$cell.types == cell)]
  
  Makers <- FindMarkers(object,ident.1 =  ident.1,ident.2 = ident.2)         ## 1 is positive, 2 is negetive (2 is control relative to 1)
  
  DEG2RNK(Makers,name=name)
  
}

#all <- initialization()     run once to create project


###### data processing

{
  all <- all[,which(all$cell.types!='?')]
  
  all <- standard_process(all)
  
  table(all$samples[all$cell.types=='Mal']) > 50
  
  samples <- c('78','79','88','71','81','80','89','194','102','110','103','106','98','129pa')
  
  samples <- paste0('Mel',samples)
  
  # samples <- unique(all$samples)
  
  Mal <- standard_process(all[,which(all$cell.types=='Mal' & all$samples %in% samples)])
  
  Mal <- RunTSNE(Mal, dims = 1:50)
  
}


###### General plotting          load Rdata and continue here

{
  ##########################       all cells
  
  Gene <- 'UBA3'
  
  unique(all$Cohort)
  
  DimPlot(all,group.by = 'cell.types')

  DimPlot(all)
  
  DimPlot(all,group.by = 'Cohort')
  
  DimPlot(all,cells = all$cell.types == 'Mal')
  
  FeaturePlot(all,features = Gene)
  
  DotPlot(all,features = Gene,group.by = 'cell.types')

  ############################          Mal cells
  
  DimPlot(Mal,reduction='umap',group.by = 'samples')

  TSNEPlot(Mal,group.by = 'samples')

  DimPlot(Mal,reduction='tsne',group.by = 'samples')

  VlnPlot(Mal,features = Gene,group.by = 'cell.types')
  
  #VlnPlot(Mal,features = c('SERPING1','SLC1A5'),group.by = 'status')
  
}






Expression <- as.matrix(Mal[Gene,]@assays$RNA@data)

half <- quantile(Expression, probs = seq(0, 1, 1/4))[3]

Mal_high <- colnames(Mal)[which(Expression >= half )]

Mal_low <- colnames(Mal)[which(Expression < half )]

Mal$status <- NA

Mal$status[Mal_high] <- paste0(Gene,'_high')

Mal$status[Mal_low] <- paste0(Gene,'_low')

################################################## optional DEG

DEG <- FindMarkers(Mal,ident.1 = Mal_high,ident.2 = Mal_low,features = rownames(Mal)) ## 1 is positive, 2 is negetive

DEG2RNK(DEG,name='tumor_cell_DEG')

######################################################

low <- paste0(Gene,'_low_Mal')
high <- paste0(Gene,'_high_Mal')

TSNEPlot(Mal,group.by = 'status')

Marker_Tb <- cbind(table(Mal[,Mal_high]$samples) ,
             table(Mal[,Mal_low]$samples))

colnames(Marker_Tb) <- c(high,low)

assay <- all[,which(all$samples %in% samples)]

cell_stastic <- table(assay$samples,assay$cell.types,dnn=c('samples','cell.types'))

df <- cbind(cell_stastic[rownames(Marker_Tb),which(colnames(cell_stastic) !='Mal')],Marker_Tb)
 
df_proportion_general <- df/rowSums(df)         ######  rowSums(df/rowSums(df))

df_proportion_specific <- df/rowSums(Marker_Tb)

plot_dt <- df_proportion_specific   # ratio of all cell types = number of a specific  cell type / total Mal cell number

#plot_dt <- df_proportion_general
#############                                                         
                            # Higher CD8 T in UAB low group, consistent with immune response pathway enrichment

plot_dt <- plot_dt[order(plot_dt[,high]),]

ggplot(melt(plot_dt), aes(Var1, Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

##################################################          mean expression as group threshold
 
low_samples <- rownames(plot_dt )[which(plot_dt[,high]<0.5) ]

high_samples <- rownames(plot_dt )[which(plot_dt[,high]>=0.5) ]

box_plot_dt <- data.frame(cbind( plot_dt,status=rep(NA,length(rownames(plot_dt)))))

box_plot_dt[low_samples,'status'] <- low 
box_plot_dt[high_samples,'status'] <- high

ggplot(melt(box_plot_dt),aes(x=status,y=value,fill=status))+
  facet_wrap(~ variable)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1) 

#### REMOVE OUTLIER SAMPLE

Samples_to_remove <- c('Mel98')

plot_dt_cleaned <- plot_dt[-which(rownames(plot_dt) %in% Samples_to_remove ),]

ggplot(melt(plot_dt_cleaned), aes(Var1, Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

box_plot_dt_cleaned <- box_plot_dt[-which(rownames(box_plot_dt) %in% Samples_to_remove),]

ggplot(melt(box_plot_dt_cleaned),aes(x=status,y=value,fill=status))+
  facet_wrap(~ variable)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1) +geom_signif(comparisons = list(c(low,high)),test='t.test')


###########

df_melted <- melt(df)

ggplot(df_melted, aes(fill=Var2,x=Var1,y=value)) + 
  geom_bar(position='stack', stat='identity')

################# DEG discovery

cellsubset_discover_DEG(all,'T.CD8','CD8T_in_UBA_HIGH_VS_UBA_LOW')
cellsubset_discover_DEG(all,'Macrophage','Macrophage_in_UBA_HIGH_VS_UBA_LOW')
cellsubset_discover_DEG(all,'B.cell','Bcell_in_UBA_HIGH_VS_UBA_LOW')
cellsubset_discover_DEG(all,'T.CD4','CD4T_in_UBA_HIGH_VS_UBA_LOW')

################ cell chat

run_cellchat <- function(obj,cell_label,cell_info,group_label,group_info){

  obj_list <- sapply(group_info, function(x) obj[,which(unlist(obj[[group_label]]) %in% x
                                                        & unlist(obj[[cell_label]]) %in% cell_info) ] )
  
  run_ <- function(obj){
    
    data.input <- obj[["RNA"]]@data
    labels <- Idents(obj)
    meta <- data.frame(labels = labels, row.names = names(labels))
    
    cellchat <- createCellChat(object = obj, group.by = cell_label)
    CellChatDB <- CellChatDB.human
    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use
    
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    #> The number of highly variable ligand-receptor pairs used for signaling inference is 692
    
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    return(cellchat)
    
  }
  cellchat_list<- sapply(obj_list,function(x) run_(x))
  return(cellchat_list)
  #return(obj_list)
}

cell_chat_visualization <- function(cellchat){
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

all.cleaned <- all[,which(all$samples %in% c('Mel98'))]

cell_subset <- c('B.cell','Macrophage','Mal','T.CD4','T.CD8','T.cell')
cellchat_list <- run_cellchat(all,'cell.types',cell_subset,'samples',list(low=low_samples,high=high_samples))
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

save(cellchat_list, file = "cellchat_object.UBA3.RData")
save(cellchat, file = "cellchat_merged_UBA3.RData")



#################################### visualization

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#######(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
## red-increase / blue decrease, in second one compared with first one
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


#(C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets

weight.max <- getMaxWeight(cellchat_list , attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cellchat_list)[i]))
}


#(D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types
group.cellType <- c('B.cell', 'Macrophage', 'Mal', rep('T',3))
group.cellType <- factor(group.cellType, levels = c('B.cell', 'Macrophage', 'Mal', 'T'))
object.list <- lapply(cellchat_list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}




#(A) Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)




#(B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T.CD8", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mal", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))




#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2






gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2





#   1  Identify dysfunctional signaling by comparing the communication probabities
pdf('differential_ligand_receptor.pdf',height = 50,width = 8)
netVisual_bubble(cellchat, sources.use = 'T.CD8', targets.use = c('Mal','Macrophage'),  comparison = c(1, 2), angle.x = 45,remove.isolate = TRUE)
dev.off()
### 2
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
