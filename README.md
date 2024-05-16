# TCMNP-Traditional Chinese Medicine Network Pharmacology

### author: Jinkun Liu, Min Ying
### time: 2024-5-16

# TCMNP <img src= https://github.com/tcmlab/image/blob/main/TCMNP.png align="right" height="200" />

## Introduction

TCMR包名称变更为TCMNP包，以避免与现有的其他TCMR包的冲突，感谢大家的使用与反馈！

Hello, here is the TCMNP package for Traditional Chinese Medicine Network Pharmacology research.
DATABASE:TCMNP
https://tcmlab.top/tcmr/

## Installation

You can install the development version of TCMNP like so:

### 1. online

```{r}
if(!require(devtools))install.packages("devtools")
if(!require(TCMNP))devtools::install_github("tcmlab/TCMNP",upgrade = FALSE,dependencies = TRUE)
```

### 2. local

Click the green button "code" on this page, then click "Download ZIP" to download it to your working directory. Install it with 

```{r}
devtools::install_local("TCMNP-main.zip", upgrade = F, dependencies = T)                                               
```

## functions
<img src= https://github.com/tcmlab/image/blob/main/TCMR-4-1-01.png height="600" />

## 0. loaded R package 

```{r}
library(TCMNP)
library(dplyr)
library(ggplot2)
library(ggraph)
library(clusterProfiler, quietly= TRUE)
library(org.Hs.eg.db, quietly= TRUE)
library(DOSE, quietly= TRUE)
```

## 1. tcm_comp

```{r}
xfbdf.compostion <- data.frame(
  herb = c(
    "mahuang", "kuxingren", "shengshigao",
    "shengyiren", "maocangzhu", "guanghuoxiang",
    "qinghaocao", "mabiancao", "ganlugen", "tinglizi",
    "huajuhong", "shenggancao", "huzhang"
  ),
  weight = c(6, 15, 30, 30, 10, 15, 12, 30, 30, 15, 15, 10, 20)
)
tcm_comp(xfbdf.compostion)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_comp.png height="400" />

```{r}
herb = c("ma huang", "ku xing ren", "hua ju hong",
     "cang zhu", "guang huo xiang", "yi yi ren",
     "hu zhang", "qing hao", "ma bian cao", "lu gen",
     "ting li zi", "shi gao", "gan cao")
 xfbdf.compostion <- data.frame(
   herb = herb,
   weight = c(6, 15, 15, 10, 15, 30, 20, 12, 30, 30, 15, 30, 10),
   property = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Property,
   flavor = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Flavor,
   meridian = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Meridian
 )
 tcm_comp_plus(xfbdf.compostion)
```
<img src= https://github.com/tcmlab/image/blob/main/tcm_comp_plus2.png height="400" />

## 2. herb_target

```{r}
herbs<-c('麻黄', '甘草','苦杏仁','石膏',
          '薏苡仁', '苍术', '青蒿', '猪苓',
          '马鞭草', '葶苈子','化橘红', 
          '虎杖', '广藿香','芦根')
xfbdf <- herb_target(herbs, type = "Herb_cn_name")
head(xfbdf)
```

```{r}
herb      molecule_id molecule target
cang zhu   MOL000173  wogonin   NOS2
cang zhu   MOL000173  wogonin  PTGS1
cang zhu   MOL000173  wogonin   ESR1
cang zhu   MOL000173  wogonin     AR
cang zhu   MOL000173  wogonin  SCN5A
cang zhu   MOL000173  wogonin  PPARG
```

```{r}
herbs2 <- c("ma huang", "ku xing ren")
fufang2 <- herb_target(herbs2, type ="Herb_name_pin_yin")
head(fufang2)
```

```{r}
herb         molecule_id molecule target
ku xing ren   MOL010921  estrone    D1R
ku xing ren   MOL010921  estrone   DRD1
ku xing ren   MOL010921  estrone  CHRM3
ku xing ren   MOL010921  estrone     F2
ku xing ren   MOL010921  estrone  CHRM1
ku xing ren   MOL010921  estrone  CHRM5
```

## 3. target_herb

```{r}
gene <- c("MAPK1", "JUN", "FOS", "RAC1", "IL1", "IL6")
herb.data <- target_herb(gene)
head(herb.data)
```

```{r}
   herb    molecule_id      molecule   target
ai di cha   MOL000422      kaempferol    JUN
ai di cha   MOL000098       quercetin    FOS
ai di cha   MOL000098       quercetin  MAPK1
ai di cha   MOL000098       quercetin    JUN
ai di cha   MOL000098       quercetin    IL6
    ai ye   MOL000358 beta-sitosterol    JUN
```
## 4.tcm_prescription 

```{r}
#通过疾病靶点基因寻找中药及其处方
#寻找到的中药
tcm_prescription(disease_data)[[1]]
# A tibble: 273 × 2
# Groups:   Herb_cn_name [273]
   Herb_cn_name  freq
   <chr>        <int>
 1 人参           437
 2 甘草           433
 3 胡芦巴         325
 4 皂角刺         288
 5 地榆           276
# ℹ 263 more rows
# ℹ Use `print(n = ...)` to see more rows
```

```{r}
#寻找到的处方
library(formattable)
formattable(tcm_prescription(disease_data)[[2]], list(Count = color_bar("lightblue")))
```
<img src= https://github.com/tcmlab/image/blob/main/tcm_pre2.png height="400" />

```{r}
# draw plot
tcm_pre_plot(tcm_prescription(disease_data)[[2]],color = 'Spectral')
```
<img src= https://github.com/tcmlab/image/blob/main/tcm_pre_plot.png height="400" />

## 5. tcm_net

```{r}
data("xfbdf", package = "TCMR")
network.data <- xfbdf %>%
  dplyr::select(herb, molecule, target) %>%
  sample_n(100, replace = FALSE) %>%
  as.data.frame()
tcm_net(network.data, 
        label.degree = 3,
        rem.dis.inter = TRUE)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_net.png height="400" />

## 6. degree_plot

```{r}
degree_plot(xfbdf,plot.set='horizontal')
```

<img src= https://github.com/tcmlab/image/blob/main/degree_plot2.png height="400" />

## 7. tcm_sankey

```{r}
sankey.data <- xfbdf[sample(nrow(xfbdf), 30), ]
 tcm_sankey(sankey.data,
   text.size = 3,
   text.position = 1
 )
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_sankey.png height="400" />

## 8. tcm_alluvial

```{r}
alluvial.data <- xfbdf[sample(nrow(xfbdf), 30), ]
 tcm_alluvial(alluvial.data,
   text.size = 3,
   text.position = 1
 )
 
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_alluvial.png height="400" />

## 9. bar_plot

```{r}
data(xfbdf, package = "TCMR")
eg <- bitr(unique(xfbdf$target), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
KK <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
bar_plot(KK,title = "KEGG")
```

<img src= https://github.com/tcmlab/image/blob/main/barplot_kk.png height="400" />

```{r}
BP <- enrichGO(
  gene = eg$ENTREZID,
  "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)
bar_plot(BP,title = "biological process")
```

<img src= https://github.com/tcmlab/image/blob/main/barplot_bp.png height="400" />

## 10. dot_plot

```{r}
dot_plot(KK, title = "KEGG")
```

<img src= https://github.com/tcmlab/image/blob/main/dot_plot_KK.png height="400" />

## 11. lollipop_plot

```{r}
lollipop_plot(KK, title = "KEGG")
```

<img src= https://github.com/tcmlab/image/blob/main/lollipop_plot_KK.png height="400" />

## 12. cir_plot

```{r}
cir_plot(KK)
```

<img src= https://github.com/tcmlab/image/blob/main/cir_plot_kk.png height="400" />



## 13. pathway_cirplot

```{r}
#Filter the pathways to display
data(KK, package = "TCMR")
newdata <- KK %>% dplyr::slice(11, 15, 17, 33, 53, .preserve = .preserve)
pathway_cirplot(newdata)
```

<img src= https://github.com/tcmlab/image/blob/main/pathway_cirplot.png height="400" />


## 14. pathway_ccplot

```{r}
pathway_ccplot(KK, root = "KEGG")
```

<img src= https://github.com/tcmlab/image/blob/main/pathway_ccplot.png height="400" />

## 15. bubble_plot

```{r}
bubble_plot(KK)
```

<img src= https://github.com/tcmlab/image/blob/main/bubble_plot.png height="400" />

## 16. dot_sankey

```{r}
dot_sankey(newdata, dot.x = 0.35, dot.y = 0.25)
```

<img src= https://github.com/tcmlab/image/blob/main/dot_sankey.png height="400" />

```r
data(kegg.filter, package = "TCMR")
#Because not all pathways and genes are what we want to display, and displaying too many genes at the same time will cause text overlap.
#The data format must be an S4 object or data frame.
#"kegg.filter" was derived from the data frame after kegg-enriched data is screened for pathways and genes.
dot_sankey(kegg.filter)
```

<img src= https://github.com/tcmlab/image/blob/main/dot_sankey2.png height="400" />

## 17. go_barplot

```{r}
eg <- bitr(unique(xfbdf$target), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go.diff <- enrichGO(
  gene = eg$ENTREZID,
  org.Hs.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  ont = "all",
  readable = T
)
go_barplot(go.diff)
```

<img src= https://github.com/tcmlab/image/blob/main/go_barplot(go.diff).png height="400" />

## 18. go_dotplot

```{r}
go_dotplot(go.diff)
```

<img src= https://github.com/tcmlab/image/blob/main/go_dotplot(go.diff).png height="400" />

## 19. go_lollipop

```{r}
go_lollipop(go.diff)
```
<img src= https://github.com/tcmlab/image/blob/main/go_lollipop(go.diff).png height="400" />

## 20. go_cir

```{r}
go_cir(go.diff)
```

<img src= https://github.com/tcmlab/image/blob/main/go_cir(go.diff).png height="400" />

## 21.  tcm_sankey_dot

```{r}
 KK2 <- KK %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
 path <- separate_rows(KK2@result, geneID, sep = "/")
 data_sankey <- left_join(xfbdf, path,
   by = c("target" = "geneID"),
   relationship = "many-to-many"
 ) %>%
   distinct() %>%
   drop_na() %>%
   sample_n(30, replace = FALSE) %>%
   as.data.frame()
 tcm_sankey_dot(data_sankey)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_sankey_dot.png height="400" />

```{r}
data_sankey2<- data_sankey %>% 
               dplyr::select(herb, molecule, target, Description)
tcm_sankey(data_sankey2,text.position = 1)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_sankey2.png height="400" />

## 22. tcm_alluvial_dot

```{r}
tcm_alluvial_dot(data_sankey)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_alluvial_dot.png height="400" />

```{r}
tcm_alluvial(data_sankey2, text.position = 1)
```

<img src= https://github.com/tcmlab/image/blob/main/tcm_alluvial2.png height="400" />

## 23. ppi_plot

```{r}
data(string, package = "TCMR")
ppi_plot(string,
    label.degree = 1,
    nodes.color = "Spectral",
    label.repel = TRUE)
```

<img src= https://github.com/tcmlab/image/blob/main/ppi_plot_2.png height="400" />

```{r}

 ppi_plot(string, nodes.color = "Spectral",
                label.degree = 1,
                label.size = 3,
                label.repel = TRUE,
                graph.layout = 'circle')
```

<img src= https://github.com/tcmlab/image/blob/main/ppi_plot3.png height="400" />

## 24. dock_plot

```{r}
data <- matrix(rnorm(81), 9, 9)
data[1:9, seq(1, 9, 2)] <- data[1:9, seq(1, 9, 2)] - 4
colnames(data) <- paste("molecule", 1:9, sep = "")
rownames(data) <- paste("target", 1:9, sep = "")
data <- round(data, 2)
dock_plot(data)
```

<img src= https://github.com/tcmlab/image/blob/main/dock_plot.png height="400" />

```{r}
dock_plot(data, shape = "circle", legend.height = 3)
```

<img src= https://github.com/tcmlab/image/blob/main/dock_plot2.png height="400" />

## 25. venn_plot

```{r}
data(venn.data, package = "TCMR")
venn_plot(venn.data, type = 1)
```

<img src= https://github.com/tcmlab/image/blob/main/venn_plot(venn%2C%20type%3D1).png height="400" />

```{r}
venn_plot(venn.data, type = 2)
```

<img src= https://github.com/tcmlab/image/blob/main/venn_plot(venn%2C%20type%3D3).png height="400" />

## 26. venn_net

```{r}
data(venn.data, package = "TCMR")
gene <- names(sort(table(venn.data$gene), decreasing = TRUE))[1:50]
data <- venn.data[venn.data$gene %in% gene, ]
data2 <- sample_n(venn.data, 100) %>% rbind(data)
venn_net(data2, label.degree = 1)
```
<img src= https://github.com/tcmlab/image/blob/main/venn_net.png height="400" />

```{r}
#targets shared by four datasets
venn_net(data2, label.degree = 4)
```
<img src= https://github.com/tcmlab/image/blob/main/venn_net4.png height="400" />

```{r}
venn_net(data2, edge.type = "hive", label.degree = 4)
```
<img src= https://github.com/tcmlab/image/blob/main/venn_net4-2.png height="400" />

## 27. tf_filter

```{r}
# Finding transcription factors and their target genes
data(venn_data, package = "TCMNP")
head(venn_result(venn_data))
```
```{r}
    Detail RNA_seq MS_protein DisGeNET XFBDT SharedSets
1          TIMP1       1          1        1     1          4
2           SOD2       1          1        1     1          4
3           MMP2       1          1        1     1          4
4            MPO       1          1        1     1          4
5           MMP9       1          1        1     1          4
6       SERPINE1       1          1        1     1          4
```

## 28. tf_filter

```{r}
# Finding transcription factors and their target genes
data(xfbdf, package = "TCMR")
newdata <- tf_filter(xfbdf$target)
head(newdata)
```
```{r}
  TF Target Mode_of_Regulation References(PMID) Species Database
  AHR  ABCG2            Unknown         20460431   Human   TRRUST
  AHR   AHRR         Activation         18848529   Human   TRRUST
  AHR   ARNT            Unknown 19255421;8631989   Human   TRRUST
  AHR  BRCA1            Unknown         18259752   Human   TRRUST
  AHR    CA9         Repression         19154183   Human   TRRUST
  AHR  CCND1            Unknown         24380854   Human   TRRUST
```

## 29. tf_cirplot

```{r}
# Visualization of filtered transcription factor data
data(xfbdf, package = "TCMR")
tf_data <- tf_filter(xfbdf$target)
set.seed(1234)
data <- tf_data[, 1:2] %>%
            distinct() %>%
            sample_n(100)
tf_cirplot(data, color = "Spectral")
```
<img src = https://github.com/tcmlab/image/blob/main/tf_cirplot_2.png height="400" />

## 30. etcm

```{r}
data("mahuang", package = "TCMR")
data <- etcm(mahuang, herb = "ma huang")
head(data)
```

```{r}
     herb   molecule target    QED
 ma huang Kaempferol   ACTB 0.9643
 ma huang Kaempferol    AHR 0.9643
 ma huang Kaempferol AKR1C1 0.8621
 ma huang Kaempferol   AKT1  0.963
 ma huang Kaempferol ATP5A1 0.9643
 ma huang Kaempferol  ATP5B 0.9643
```
