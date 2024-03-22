# ScRNAseq/Tumor-immune-infiltration
## Description
scRNA analysis for studying how the expression of oncogenes affects tumor-immune cell interplay

## Overview

![Picture2](https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/f734f0cc-8047-4846-95db-835323550cd2)


## Software dependencies :
R 4.3.1


## Material :
### Sample dataset : GSE115978   Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes   https://doi.org/10.1016/j.cell.2018.09.006
### Use FN1 as a target gene

## Steps :
### 1. Generate Seurat object with count matrix & cell annotation

### 2. Exclude tumor samples that contain less than 50 malignant cells

### 3. Run TSNE to check tumor diversity 
<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/b60b3835-5891-4d98-ad20-764e09cfce9b" width="200" />

### 4. Assign tumor cell identities based on FN1 expression (default uses mean FN1 expression as threshold)
```
half <- quantile(Expression, probs = seq(0, 1, 1/4))[3]

Mal_high <- colnames(Mal)[which(Expression >= half )]

Mal_low <- colnames(Mal)[which(Expression < half )]
```
<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/eb430e6a-3cc6-484a-a0b9-58c863c9668b" width="200" />

### 5. Run DEG discovery then generate rnk file for pre-ranked GSEA  
```
DEG <- FindMarkers(Mal,ident.1 = Mal_high,ident.2 = Mal_low,features = rownames(Mal)) ## 1 is positive, 2 is negetive
DEG2RNK(DEG,name='tumor_cell_DEG')
```

### 6. Generate immune-infiltration matrix (score of each line = number of corresponding cell in that sample / total malignant cell in that sample)
```
ggplot(melt(plot_dt), aes(Var1, Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
```
<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/249389bb-ad5f-4729-acab-b6c8c1719bc1" width="400" />

<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/16957c15-c59a-43ed-9caa-ffd46ecdfb12" width="200" />




### 7. Assign tumor sample identities based on proportion of high-level marker-expressing cells (default use 50% as threshold)
```
low_samples <- rownames(plot_dt )[which(plot_dt[,high]<0.5) ]
high_samples <- rownames(plot_dt )[which(plot_dt[,high]>=0.5) ]

```
Visualize the infiltration difference between tumor subgroups

<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/c01df1a7-40b1-4c9f-9f7d-407cd2a6e4c1" width="400" />

The result may help indicate the immune cells that are affected by the expression of FN1

### 8. Run immune cell functionality analysis with DEG discovery (for immune cell subsets in high-FN1-expression tumor vs low-FN1-expression tumor)
```
cellsubset_discover_DEG(all,'T.CD8','CD8T_in_FN1_HIGH_VS_FN1_LOW')
cellsubset_discover_DEG(all,'Macrophage','Macrophage_in_FN1_HIGH_VS_FN1_LOW')
cellsubset_discover_DEG(all,'B.cell','Bcell_in_FN1_HIGH_VS_FN1_LOW')
cellsubset_discover_DEG(all,'T.CD4','CD4T_in_FN1_HIGH_VS_FN1_LOW')
```
Create .rnk files ready for pre-ranked GSEA

### 9. Cell chat analysis
Official tutorial : https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

<img src="https://github.com/Gico1941/Tumor-Immune-infiltration-discovery/assets/127346166/8453cfff-4ffe-4ad7-bef1-101bcc7a387f" width="400" />


