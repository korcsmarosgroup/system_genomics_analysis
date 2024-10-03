#CD
#Perform GO analysis of TaMMA dataset for CD 

#First perform GO analysis of DEGs from CD colon vs healthy control colon
CD_colon_tamma<-read.csv("Colon_CD-vs-Colon_Control.diffexp.tsv",  sep="\t")
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(GO.db)
library(GOSemSim)

# we want the log2 fold change 
original_gene_list_CD_colon_tamma<- CD_colon_tamma$log2FoldChange
length(original_gene_list_CD_colon_tamma)
# name the vector
names(original_gene_list_CD_colon_tamma) <-  CD_colon_tamma$Gene

# omit any NA values 
gene_list_CD_colon_tamma<-na.omit(original_gene_list_CD_colon_tamma)
length(gene_list_CD_colon_tamma)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_CD_colon_tamma = sort(gene_list_CD_colon_tamma, decreasing = TRUE)

#GSEA
gse_CD_colon_tamma <- gseGO(geneList=gene_list_CD_colon_tamma, 
                            ont ="BP", 
                            keyType = "SYMBOL", 
                            minGSSize = 3, 
                            maxGSSize = 800, 
                            pvalueCutoff = 0.05, 
                            verbose = TRUE,  
                            OrgDb = organism, 
                            pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse_CD_colon_tamma, showCategory=20, split=".sign", font.size=6) + facet_grid(.~.sign)
dotplot(gse_CD_colon_tamma, showCategory=20,font.size=6)

#GO analysis of TaMMA colon
fc_genes_CD_colon_tamma<-(CD_colon_tamma[CD_colon_tamma$log2FoldChange>1,])
fc_genes_CD_colon_tamma<-fc_genes_CD_colon_tamma[fc_genes_CD_colon_tamma$padj<0.05,]
fc_genelist_CD_colon_tamma<-fc_genes_CD_colon_tamma$Gene
GO_results_CD_colon_tamma<-enrichGO(gene = fc_genelist_CD_colon_tamma, universe = names(gene_list_CD_colon_tamma), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "BP", pAdjustMethod = "fdr")

GO_barplot_CD_colon_tamma<- plot(barplot(GO_results_CD_colon_tamma, showCategory = 20, font.size = 6))
#GO_barplot_CD_colon_tamma

#Extract significant GO terms of TaMMA colon
cluster_summary_CD_colon<-data.frame(GO_results_CD_colon_tamma)
signficant_GO_terms_CD_colon_tamma<-cluster_summary_CD_colon[cluster_summary_CD_colon$p.adjust<0.05,]
#signficant_GO_terms_CD_colon_tamma$ID

print("First done")

#Second perform GO analysis of DEGs from CD rectum vs healthy control rectum
CD_rectum_tamma<-read.csv("Rectum_CD-vs-Rectum_Control.diffexp.tsv",  sep="\t")

# we want the log2 fold change 
original_gene_list_CD_rectum_tamma<- CD_rectum_tamma$log2FoldChange
length(original_gene_list_CD_rectum_tamma)

# name the vector
names(original_gene_list_CD_rectum_tamma) <-  CD_rectum_tamma$Gene

# omit any NA values 
gene_list_CD_rectum_tamma<-na.omit(original_gene_list_CD_rectum_tamma)
length(gene_list_CD_rectum_tamma)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_CD_rectum_tamma = sort(gene_list_CD_rectum_tamma, decreasing = TRUE)

#GSEA
gse_CD_rectum_tamma <- gseGO(geneList=gene_list_CD_rectum_tamma, 
                             ont ="BP", 
                             keyType = "SYMBOL", 
                             minGSSize = 3, 
                             maxGSSize = 800, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE,  
                             OrgDb = organism, 
                             pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse_CD_rectum_tamma, showCategory=20, split=".sign", font.size=6) + facet_grid(.~.sign)
dotplot(gse_CD_rectum_tamma, showCategory=20,font.size=6)

#GO analysis of TaMMA rectum
fc_genes_CD_rectum_tamma<-(CD_rectum_tamma[CD_rectum_tamma$log2FoldChange>1,])
fc_genes_CD_rectum_tamma<-fc_genes_CD_colon_tamma[fc_genes_CD_rectum_tamma$padj<0.05,]
fc_genelist_CD_rectum_tamma<-fc_genes_CD_rectum_tamma$Gene
GO_results_CD_rectum_tamma<-enrichGO(gene = fc_genelist_CD_rectum_tamma, universe = names(gene_list_CD_rectum_tamma), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "BP", pAdjustMethod = "fdr")

GO_barplot_CD_rectum_tamma<- plot(barplot(GO_results_CD_rectum_tamma, showCategory = 20, font.size = 6))
#GO_barplot_CD_rectum_tamma

#Extract significant GO terms of TaMMA rectum
cluster_summary_rectum<-data.frame(GO_results_CD_rectum_tamma)
signficant_GO_terms_CD_rectum_tamma<-cluster_summary_rectum[cluster_summary_rectum$p.adjust<0.05,]
#signficant_GO_terms_CD_rectum_tamma$ID

print("Second done")

#Third perform GO analysis of DEGs from CD ileum vs healthy control ileum
CD_ileum_tamma<-read.csv("Ileum_CD-vs-Ileum_Control.diffexp.tsv",  sep="\t")

# we want the log2 fold change 
original_gene_list_CD_ileum_tamma<- CD_ileum_tamma$log2FoldChange
length(original_gene_list_CD_ileum_tamma)

# name the vector
names(original_gene_list_CD_ileum_tamma) <-  CD_ileum_tamma$Gene

# omit any NA values 
gene_list_CD_ileum_tamma<-na.omit(original_gene_list_CD_ileum_tamma)
length(gene_list_CD_ileum_tamma)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_CD_ileum_tamma = sort(gene_list_CD_ileum_tamma, decreasing = TRUE)

#GSEA
gse_CD_ileum_tamma <- gseGO(geneList=gene_list_CD_ileum_tamma, 
                            ont ="BP", 
                            keyType = "SYMBOL", 
                            minGSSize = 3, 
                            maxGSSize = 800, 
                            pvalueCutoff = 0.05, 
                            verbose = TRUE,  
                            OrgDb = organism, 
                            pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse_CD_ileum_tamma, showCategory=20, split=".sign", font.size=6) + facet_grid(.~.sign)
dotplot(gse_CD_ileum_tamma, showCategory=20,font.size=6)

#GO analysis of TaMMA ileum
fc_genes_CD_ileum_tamma<-(CD_ileum_tamma[CD_ileum_tamma$log2FoldChange>1,])
fc_genes_CD_ileum_tamma<-fc_genes_CD_ileum_tamma[fc_genes_CD_ileum_tamma$padj<0.05,]
fc_genelist_CD_ileum_tamma<-fc_genes_CD_ileum_tamma$Gene
GO_results_CD_ileum_tamma<-enrichGO(gene = fc_genelist_CD_ileum_tamma, universe = names(gene_list_CD_ileum_tamma), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "BP", pAdjustMethod = "fdr")

GO_barplot_CD_ileum_tamma<- plot(barplot(GO_results_CD_ileum_tamma, showCategory = 20, font.size = 6))
GO_barplot_CD_ileum_tamma

#Extract significant GO terms of TaMMA ileum
cluster_summary_ileum<-data.frame(GO_results_CD_ileum_tamma)
signficant_GO_terms_CD_ileum_tamma<-cluster_summary_ileum[cluster_summary_ileum$p.adjust<0.05,]
signficant_GO_terms_CD_ileum_tamma$ID


#Combine GO terms from CD TaMMA rectum and colon and remove duplicates
significant_GO_terms_CD_tamma<-c(signficant_GO_terms_CD_rectum_tamma$ID, signficant_GO_terms_CD_colon_tamma$ID, signficant_GO_terms_CD_ileum_tamma$ID)
significant_GO_terms_CD_tamma_final<-significant_GO_terms_CD_tamma[!duplicated(significant_GO_terms_CD_tamma)]
significant_GO_terms_CD_tamma_final<-as.data.frame(significant_GO_terms_CD_tamma_final)
colnames(significant_GO_terms_CD_tamma_final)<-"GOID"

print("Third done")

#Read GO switches file
CD_GO_switches<-read.csv("17_11_2021_CDLeuven_network_all_snps_23_07_2024.tsv",  sep="\t")


#Switch GO terms to GO IDs in CD GO switches file

#Load libraries

goIdToTerm <- function(x, names = TRUE, keepNA = TRUE) {
  stopifnot(requireNamespace("GO.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ans <- rep(NA_character_, length(x))
  names(ans) <- x
  ids <- AnnotationDbi::GOID(GO.db::GOTERM)
  i <- match(x, ids)
  k <- which(!is.na(i))
  res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
  ans[k] <- res
  if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
  if (!names) names(ans) <- NULL
  return(ans)
}

goTermToId <- function(x, names = TRUE, keepNA = TRUE) {
  stopifnot(requireNamespace("GO.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ans <- rep(NA_character_, length(x))
  names(ans) <- x
  terms <- AnnotationDbi::Term(GO.db::GOTERM)
  i <- match(x, terms)
  k <- which(!is.na(i))
  res <- AnnotationDbi::GOID(GO.db::GOTERM[i[k]])
  ans[k] <- res
  if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
  if (!names) names(ans) <- NULL
  return(ans)
}


##' @rdname goIdToTerm
flipGoTermId <- function(x, names = TRUE, keepNA = TRUE) {
  isId <- grepl("GO:", x)
  if (any(isId)) ans <- goIdToTerm(x, names, keepNA)
  else ans <- goTermToId(x, names, keepNA)
  return(ans)
}

##' @rdname goIdToTerm
prettyGoTermId <- function(x) {
  y <- flipGoTermId(x)
  if (any(grepl("GO:", x))) ans <- paste0(y, " (", x, ")")
  else ans <- paste0(x, " (", y, ")")
  return(ans)
}



#New GO switches file (post-Revigo analysis) and then clean 
CD_GO_switches<-read.csv("17_11_2021_CDLeuven_network_all_snps_23_07_2024.tsv",  sep="\t")

length(CD_GO_switches$Target)
CD_GO_switches$Source <- gsub("healthy", "", CD_GO_switches$Source)
CD_GO_switches$Source <- gsub("_", "", CD_GO_switches$Source)
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

CD_GO_switches$Source<-trim(CD_GO_switches$Source)

CD_GO_switches$Target <- gsub("nonhealthy", "", CD_GO_switches$Target)
CD_GO_switches$Target <- gsub("_", " ", CD_GO_switches$Target)
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}
CD_GO_switches$Target<-trim(CD_GO_switches$Target)

Healthy_GOID<-flipGoTermId(CD_GO_switches$Source)
CD_GO_switches_2<-cbind(CD_GO_switches, Healthy_GOID)
IBD_GOID<-flipGoTermId(CD_GO_switches_2$Target)
CD_GO_switches_3<-cbind(CD_GO_switches_2, IBD_GOID)


#Find TaMMA CD GO terms in Target GO switches
#CD_GO_switches_tamma_1<-CD_GO_switches_3[CD_GO_switches_3$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]

#Find TaMMA CD GO terms in Source GO switches
#CD_GO_switches_tamma_2<-CD_GO_switches_tamma_1[CD_GO_switches_tamma_1$Healthy_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]


#GO similarity analysis for each pair of GO switch terms
CD_GO_switches_3$Healthy_GOID
CD_GO_switches_3$IBD_GOID

end_number <- length(CD_GO_switches_3$Healthy_GOID)
my_range<-1:end_number

e <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
e


GO_sim_score<-NULL

for(i in my_range) {
  source <- CD_GO_switches_3$Healthy_GOID[i]
  target <- CD_GO_switches_3$IBD_GOID[i]
  output<-goSim(source, target, semData=e, measure="Lin")
  GO_sim_score<-rbind(GO_sim_score, output)
}

GO_sim_score<-as.data.frame(GO_sim_score)
CD_GO_switches_4<-cbind(CD_GO_switches_3, GO_sim_score$V1)
colnames(CD_GO_switches_4)[7]="GO_sim_score"
colnames(CD_GO_switches_4)[1]="Healthy_GO_term"
colnames(CD_GO_switches_4)[2]="IBD_GO_term"


#Filter out GO terms that have a GO_sim_score of >=0.5, >=0.4, >=0.3 and >=0.2 respectively
CD_GO_switches_filtered<-CD_GO_switches_4[CD_GO_switches_4$GO_sim_score<0.5,]
# CD_GO_switches_filtered_2<-CD_GO_switches_4[CD_GO_switches_4$GO_sim_score<0.4,]
# CD_GO_switches_filtered_3<-CD_GO_switches_4[CD_GO_switches_4$GO_sim_score<0.3,]
# CD_GO_switches_filtered_4<-CD_GO_switches_4[CD_GO_switches_4$GO_sim_score<0.2,]
# CD_GO_switches_filtered_5<-CD_GO_switches_4[CD_GO_switches_4$GO_sim_score<0.15,]



#Come back to CD GO switches file filtered by GO similarity score (choose level of filter), and filter out IBD GO terms that are not present in above CD TaMMA GO term list
CD_GO_switches_filtered_tamma<-CD_GO_switches_filtered[CD_GO_switches_filtered$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]
# CD_GO_switches_filtered_tamma_2<-CD_GO_switches_filtered_2[CD_GO_switches_filtered_2$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]
# CD_GO_switches_filtered_tamma_3<-CD_GO_switches_filtered_3[CD_GO_switches_filtered_3$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]
# CD_GO_switches_filtered_tamma_4<-CD_GO_switches_filtered_4[CD_GO_switches_filtered_4$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]
# CD_GO_switches_filtered_tamma_5<-CD_GO_switches_filtered_5[CD_GO_switches_filtered_5$IBD_GOID %in% significant_GO_terms_CD_tamma_final$GOID,]

#Save filtered dataframe as csv file 
write.csv(CD_GO_switches_filtered_tamma, "CD_GO_switches_filtered_tamma_23_07_2024.csv", row.names=FALSE)
