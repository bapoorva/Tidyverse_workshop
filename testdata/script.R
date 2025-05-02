install.packages("tidyverse")
library(tidyverse)

setwd('~/Desktop/testdata/')
test=read.csv("~/Desktop/testdata/test_meta.csv")

#Lets subset the data so that its easier to view the data

df=test[1:11,]
df
#### Select: Subset by variable (columns) ####
#single column
select(df,Sample)
#multiple column
select(df,Sample,CRS)
#Range of cols
select(df,Group:possible_anno)

select(df,Sample,rpca_clusters:possible_anno)
#remove a col
select(df,-CRS)
select(df,-nCount_RNA:-Group)

select(df,starts_with('S'))

select(df,starts_with('r') & ends_with('s'))
select(df,contains('RNA'))
select(df,contains('RNA') | ends_with('ay'))
select(df,matches("[abcd]a"))
select(df,num_range("Day", 0:5))


#### filter: Subset observations (rows) ####

filter(df,possible_anno=='CD16')

filter(df,possible_anno=='CD16' & rpca_clusters==8)

filter(df,possible_anno=='CD16' | rpca_clusters==0)

#### Pipes ####
# %>% is called a pipe 
#Pipes send the output of one function as the first argument to the next function. 

df1= select(df,Sample,CRS,possible_anno)
df1
df2= filter(df1,possible_anno=='CD16')
df2

df %>% select(Sample,CRS,possible_anno) %>% filter(possible_anno=='CD16')


#### Mutate ####
df %>% mutate(new_col=paste0("Grade",CRS,sep=""))

df %>% mutate(new_col=Day+5)

df %>% mutate(new_col=nCount_RNA- mean(nCount_RNA, na.rm = TRUE))

#### Arrange Sort column by observations ####

df %>% arrange(rpca_clusters)
df %>% arrange(-rpca_clusters)
df %>% arrange(rpca_clusters,-nFeature_RNA) 
df %>% arrange(rpca_clusters,-nFeature_RNA) %>% head(n=5)
df %>% arrange(rpca_clusters,-nFeature_RNA) %>% tail(n=5)

#### Rename ####
df %>% rename("celltype"="possible_anno")
df %>% rename("celltype"="possible_anno","clusters"="rpca_clusters")


#### summarize #####
df %>% summarise(mean(nCount_RNA),sd(nCount_RNA))
df %>% summarise(n=n())

# group_by ----------------------------------------------------------------

df %>% group_by(rpca_clusters) %>% summarise(n=n())
test %>% group_by(rpca_clusters) %>% summarise(n=n())
test %>% group_by(possible_anno) %>% summarise(n=n())

test %>% group_by(Sample,possible_anno) %>% summarise(n())

df2=test %>% group_by(Sample) %>% summarise(mean(nCount_RNA),sd(nCount_RNA))
df2
test %>% summarise(n=n_distinct(possible_anno))
test %>% summarise(n=n_distinct(Group))
test %>% group_by(Sample) %>% summarise(min(nFeature_RNA))


#### join ####
#Open a 2nd file and "merge' results. 
anno=read.csv("~/Desktop/testdata/test_anno.csv")

dim(anno)
dim(df)

#using join, there's multiple types of join.We'll look at inner_join and left/right join 
inner_join(df,anno,by='rpca_clusters')
left_join(df,anno,by='rpca_clusters')
right_join(df,anno,by='rpca_clusters')
full_join(df,anno,by='rpca_clusters')

#### pivot ####
# use to be called gather and spread

#Go from long to wide
test %>% group_by(Sample,possible_anno) %>% summarise(n())
df2=test %>% group_by(Sample,possible_anno) %>% summarise(count=n()) %>% pivot_wider(names_from = possible_anno, values_from = count)
df2
#Go from wide to long
df2 %>% pivot_longer(!Sample, names_to = "celltype", values_to = "count")

#### Rownames to column and vice versa####
df$Sample=make.names(df$Sample,unique=T)
df

df %>% column_to_rownames('Sample')
df %>% rownames_to_column('Sample')

## Slice family  function to subset rows by "position"

df %>% slice(1)
df %>% slice(6:8)

df %>% slice_max(nCount_RNA,n=2)
df %>% slice_min(nFeature_RNA,n=3)

df %>% slice_head(n=5)
df %>% slice_tail(n=4)

df %>% group_by(possible_anno) %>% slice_max(nCount_RNA,n=1) 

#We can even pipe to ggplot as well. 

test %>% ggplot(.,aes(x=possible_anno,y=nCount_RNA,fill=possible_anno)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
+ facet_grid(~Group) 

test %>% group_by(Sample,possible_anno) %>% summarise(count=n())  %>% ggplot(.,aes(x=possible_anno,y=count,fill=possible_anno)) + geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#### Real use examples ####
#Example 1 : Looking at celltype proportions per group
ct=as.data.frame(table(test$Group))
total_min=ct$Freq[ct$Var1=="Minimal"]
total_sev=ct$Freq[ct$Var1=="Severe"]

test1=test %>% group_by(possible_anno,Group) %>% count() 
test1=test1 %>% filter(possible_anno %in% c("CD8+ T Cells","NK Cells","CD4+ T Cells","NKT","CD16","CD14","B Cells","intermediate monocytes"))
test1=test1 %>% mutate(Proportion=ifelse(Group=="Minimal",(n/total_min)*100,(n/total_sev)*100)) 
test1=test1 %>%   mutate(Group_prop=Group) %>% pivot_wider(.,id_cols = possible_anno, names_from = Group, values_from = c("n", "Proportion")) 


cnt_prop_table=test %>% group_by(possible_anno,Group) %>% count() %>% 
  filter(possible_anno %in% c("CD8+ T Cells","NK Cells","CD4+ T Cells","NKT","CD16","CD14","B Cells","intermediate monocytes"))%>% 
  mutate(Proportion=ifelse(Group=="Minimal",(n/total_min)*100,(n/total_sev)*100)) %>% 
  mutate(Group_prop=Group) %>% pivot_wider(.,id_cols = possible_anno, names_from = Group, values_from = c("n", "Proportion")) 

#Example 2 : Calculation CD8:CD4 ratio in both groups
#
test%>% group_by(possible_anno,Sample,Group) %>% filter(possible_anno %in% c("CD8+ T Cells","CD4+ T Cells")) %>% count()  %>% spread(possible_anno,n) %>% pivot_wider(names_from='Group',values_from=c("CD8+ T Cells","CD4+ T Cells")) %>% mutate(min.ratio=`CD4+ T Cells_Minimal`/`CD8+ T Cells_Minimal`,sev.ratio=`CD4+ T Cells_Severe`/`CD8+ T Cells_Severe`)

#Example 3 : Make a composition stacked bar plot
scrna=readRDS('Seurat_downsampled.RDS')
scrna@meta.data %>% 
  group_by(orig.ident,final_anno) %>%
  summarise(n = n()) %>%
  mutate(pct = n / sum(n)) %>% 
  ggplot(.,aes(x=final_anno,y=pct,fill=orig.ident)) +  
  geom_bar(position="fill", stat="identity")  + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Paired")+ theme_bw() + theme(legend.position = 'bottom') 

#Example 4 : Make a volcano plot with top 10 genes labeled
load("../../Downloads/NGSViewer-master/data/MILR1_bulk.RData")
diff_df=results$limma$NALM6_withMILR1_vs_NALM6_noMILR1
cols <- c("Upregulated with MILR1" = "red4", "Upregulated without MILR1" = "royalblue3", "NotSignificant" = "grey") 
diff_df$group <- "NotSignificant"
diff_df[diff_df$adj.P.Val < 0.05 & diff_df$fc < -2 ,"group"] <- "Upregulated without MILR1"
diff_df[diff_df$adj.P.Val < 0.05 & diff_df$fc > 2,"group"] <- "Upregulated with MILR1"

#n=20 #Top genes
top_genes <- diff_df %>% filter(logFC> 2 & adj.P.Val < 0.05) %>% drop_na(any_of('SYMBOL'))%>%arrange(adj.P.Val) %>%head(10)
bottom_genes <- diff_df %>% filter(logFC < -2 & adj.P.Val < 0.05) %>% drop_na(any_of('SYMBOL'))%>%arrange(adj.P.Val) %>%head(10)
top_peaks <- rbind(top_genes, bottom_genes)

#png(paste0("plots/",samp,"_volcano.png",sep=""),width = 10,height = 10,units="in",res =250)
gg=ggplot(data = diff_df, aes(x = logFC, y = -log10(adj.P.Val),col=group)) + geom_point() +
  ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") + scale_colour_manual(values = cols)+
  geom_label_repel(data = top_peaks,max.overlaps = 20, # Add labels last to appear as the top layer  
                   aes(label = SYMBOL),size=5,
                   label.size=NA)+ theme_classic()+ggtitle("NALM6 +/- MILR1")

#### Exercise ####
#Use the same test data and perform the following actions
# 1.Rename the column CRS to CRS_Grade (Hint:rename)
# 2.Get all rows for Sample9 (Hint:filter)
# 3.Day column has negative values. Make that column all positive by adding 1 to it(Hint:mutate)
# 4.Find the mean nCount_RNA for Sample3(Hint:filter+summarize)
# 5.Save a new data frame with cluster numbers and corresponding celltype similar to the test_anno file provided(Hint:select)
# 6.Make a barplot of number of minimal and severe cell count in test(Hint:group_by+summarize+ggplot)
# 7.Which sample in Severe group has the cell with highest nCount_RNA?(Hint:filter+group+Arrange)

#Tad bit more complex (ask me for dummy data if you want to attempt)
# 1.Take the DEG list, compute delta values for pct.1 and pct.2 (mutate), arrange by +ve fold change and delta (Arrange) and pick top 10 genes
# 2.You have metadata for a single cell project. Read the csv file, make X into rownames. You have some new labeling information after further analysis. Join the new labels to 
# your metadata with a different column name and create new dataframe with just barcodes, clusters at all resolutions and annotation
# 3.Take the bulk DEG data, filter for p.value < 0.05, convert fold change to Log2 FC,map the ENSEMBL id's to gene symbols and select 
#for only protein coding gene located in chromosome 12