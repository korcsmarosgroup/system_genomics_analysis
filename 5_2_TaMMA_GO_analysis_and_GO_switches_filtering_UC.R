#UC
#Perform GO analysis of TaMMA dataset for UC 

#First perform GO analysis of DEGs from UC colon vs healthy control colon
UC_colon_tamma<-read.csv("Colon_UC-vs-Colon_Control.diffexp.tsv",  sep="\t")
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(GO.db)
library(GOSemSim)

# we want the log2 fold change 
original_gene_list_UC_colon_tamma<- UC_colon_tamma$log2FoldChange
length(original_gene_list_UC_colon_tamma)
# name the vector
names(original_gene_list_UC_colon_tamma) <-  UC_colon_tamma$Gene

# omit any NA values 
gene_list_UC_colon_tamma<-na.omit(original_gene_list_UC_colon_tamma)
length(gene_list_UC_colon_tamma)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_UC_colon_tamma = sort(gene_list_UC_colon_tamma, decreasing = TRUE)

#GSEA
gse_UC_colon_tamma <- gseGO(geneList=gene_list_UC_colon_tamma, 
                            ont ="BP", 
                            keyType = "SYMBOL", 
                            minGSSize = 3, 
                            maxGSSize = 800, 
                            pvalueCutoff = 0.05, 
                            verbose = TRUE,  
                            OrgDb = organism, 
                            pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse_UC_colon_tamma, showCategory=20, split=".sign", font.size=6) + facet_grid(.~.sign)
dotplot(gse_UC_colon_tamma, showCategory=20,font.size=6)

#GO analysis of TaMMA colon upregulated in UC
fc_genes_UC_colon_tamma<-(UC_colon_tamma[UC_colon_tamma$log2FoldChange>1,])
fc_genes_UC_colon_tamma<-fc_genes_UC_colon_tamma[fc_genes_UC_colon_tamma$padj<0.05,]
fc_genelist_UC_colon_tamma<-fc_genes_UC_colon_tamma$Gene
GO_results_UC_colon_tamma<-enrichGO(gene = fc_genelist_UC_colon_tamma, universe = names(gene_list_UC_colon_tamma), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "BP", pAdjustMethod = "fdr")

GO_barplot_UC_colon_tamma<- plot(barplot(GO_results_UC_colon_tamma, showCategory = 20, font.size = 6))
#GO_barplot_UC_colon_tamma

#Extract significant GO terms of TaMMA colon upregulated in UC 
cluster_summary_UC_colon<-data.frame(GO_results_UC_colon_tamma)
signficant_GO_terms_UC_colon_tamma<-cluster_summary_UC_colon[cluster_summary_UC_colon$p.adjust<0.05,]
signficant_GO_terms_UC_colon_tamma$ID

#Second perform GO analysis of DEGs from UC rectum vs healthy control rectum
UC_rectum_tamma<-read.csv("Rectum_UC-vs-Rectum_Control.diffexp.tsv",  sep="\t")

# we want the log2 fold change 
original_gene_list_UC_rectum_tamma<- UC_rectum_tamma$log2FoldChange
length(original_gene_list_UC_rectum_tamma)

# name the vector
names(original_gene_list_UC_rectum_tamma) <-  UC_rectum_tamma$Gene

# omit any NA values 
gene_list_UC_rectum_tamma<-na.omit(original_gene_list_UC_rectum_tamma)
length(gene_list_UC_rectum_tamma)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_UC_rectum_tamma = sort(gene_list_UC_rectum_tamma, decreasing = TRUE)

#GSEA
gse_UC_rectum_tamma <- gseGO(geneList=gene_list_UC_rectum_tamma, 
                             ont ="BP", 
                             keyType = "SYMBOL", 
                             minGSSize = 3, 
                             maxGSSize = 800, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE,  
                             OrgDb = organism, 
                             pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse_UC_rectum_tamma, showCategory=20, split=".sign", font.size=6) + facet_grid(.~.sign)
dotplot(gse_UC_rectum_tamma, showCategory=20,font.size=6)

#GO analysis of TaMMA rectum
fc_genes_UC_rectum_tamma<-(UC_rectum_tamma[UC_rectum_tamma$log2FoldChange>1,])
fc_genes_UC_rectum_tamma<-fc_genes_UC_rectum_tamma[fc_genes_UC_rectum_tamma$padj<0.05,]
fc_genelist_UC_rectum_tamma<-fc_genes_UC_rectum_tamma$Gene
GO_results_UC_rectum_tamma<-enrichGO(gene = fc_genelist_UC_rectum_tamma, universe = names(gene_list_UC_rectum_tamma), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "BP", pAdjustMethod = "fdr")

GO_barplot_UC_rectum_tamma<- plot(barplot(GO_results_UC_rectum_tamma, showCategory = 20, font.size = 6))
GO_barplot_UC_rectum_tamma

#Extract significant GO terms of TaMMA UC rectum
cluster_summary_UC_rectum<-data.frame(GO_results_UC_rectum_tamma)
signficant_GO_terms_UC_rectum_tamma<-cluster_summary_UC_rectum[cluster_summary_UC_rectum$p.adjust<0.05,]
signficant_GO_terms_UC_rectum_tamma$ID

# session_info()

#Combine GO terms from UC TaMMA rectum and colon and remove duplicates
significant_GO_terms_UC_tamma<-c(signficant_GO_terms_UC_rectum_tamma$ID, signficant_GO_terms_UC_colon_tamma$ID)
significant_GO_terms_UC_tamma_final<-significant_GO_terms_UC_tamma[!duplicated(significant_GO_terms_UC_tamma)]
significant_GO_terms_UC_tamma_final<-as.data.frame(significant_GO_terms_UC_tamma_final)
colnames(significant_GO_terms_UC_tamma_final)<-"GOID"



#New GO switches file (post-Revigo analysis) and then clean 
UC_GO_switches<-read.csv("05_11_2021_UCLeuven_network_all_snps_23_07_2024.tsv",  sep="\t")

length(UC_GO_switches$Target)
UC_GO_switches$Source <- gsub("healthy", "", UC_GO_switches$Source)
UC_GO_switches$Source <- gsub("_", "", UC_GO_switches$Source)
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

UC_GO_switches$Source<-trim(UC_GO_switches$Source)

UC_GO_switches$Target <- gsub("nonhealthy", "", UC_GO_switches$Target)
UC_GO_switches$Target <- gsub("_", " ", UC_GO_switches$Target)
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}
UC_GO_switches$Target<-trim(UC_GO_switches$Target)


#Switch GO terms to GO IDs in new GO switches file

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

Healthy_GOID<-flipGoTermId(UC_GO_switches$Source)
UC_GO_switches_2<-cbind(UC_GO_switches, Healthy_GOID)
IBD_GOID<-flipGoTermId(UC_GO_switches_2$Target)
UC_GO_switches_3<-cbind(UC_GO_switches_2, IBD_GOID)


#Next perform semantic similarity analysis
UC_GO_switches_3$Healthy_GOID
UC_GO_switches_3$IBD_GOID

end_number <- length(UC_GO_switches_3$Healthy_GOID)
my_range<-1:end_number

d <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
d


GO_sim_score_2<-NULL

for(i in my_range) {
  source <- UC_GO_switches_3$Healthy_GOID[i]
  target <- UC_GO_switches_3$IBD_GOID[i]
  output<-goSim(source, target, semData=d, measure="Lin")
  GO_sim_score_2<-rbind(GO_sim_score_2, output)
}

GO_sim_score_2<-as.data.frame(GO_sim_score_2)
UC_GO_switches_4<-cbind(UC_GO_switches_3, GO_sim_score_2$V1)
colnames(UC_GO_switches_4)[7]="GO_sim_score"
colnames(UC_GO_switches_4)[1]="Healthy_GO_term"
colnames(UC_GO_switches_4)[2]="IBD_GO_term"

#Filter out GO terms that have a GO_sim_score of >=0.5, >=0.4, >=0.3 and >=0.2 respectively
UC_GO_switches_filtered<-UC_GO_switches_4[UC_GO_switches_4$GO_sim_score<0.5,]
# UC_GO_switches_filtered_2<-UC_GO_switches_4[UC_GO_switches_4$GO_sim_score<0.4,]
# UC_GO_switches_filtered_3<-UC_GO_switches_4[UC_GO_switches_4$GO_sim_score<0.3,]
# UC_GO_switches_filtered_4<-UC_GO_switches_4[UC_GO_switches_4$GO_sim_score<0.2,]
# UC_GO_switches_filtered_5<-UC_GO_switches_4[UC_GO_switches_4$GO_sim_score<0.15,]



#Come back to UC GO switches file filtered by GO similarity score (choose level of filter), and filter out IBD GO terms that are not present in above UC TaMMA GO term list
#Selected UC_GO_switches_filtered file (i.e. similarity score<0.5)
UC_GO_switches_filtered_tamma<-UC_GO_switches_filtered[UC_GO_switches_filtered$IBD_GOID %in% significant_GO_terms_UC_tamma_final$GOID,]
# UC_GO_switches_filtered_tamma_2<-UC_GO_switches_filtered_2[UC_GO_switches_filtered_2$IBD_GOID %in% significant_GO_terms_UC_tamma_final$GOID,]
# UC_GO_switches_filtered_tamma_3<-UC_GO_switches_filtered_3[UC_GO_switches_filtered_3$IBD_GOID %in% significant_GO_terms_UC_tamma_final$GOID,]
# UC_GO_switches_filtered_tamma_4<-UC_GO_switches_filtered_4[UC_GO_switches_filtered_4$IBD_GOID %in% significant_GO_terms_UC_tamma_final$GOID,]
# UC_GO_switches_filtered_tamma_5<-UC_GO_switches_filtered_5[UC_GO_switches_filtered_5$IBD_GOID %in% significant_GO_terms_UC_tamma_final$GOID,]

#Save filtered dataframe as csv file 
write.csv(UC_GO_switches_filtered_tamma, "UC_GO_switches_filtered_tamma_23_07_2024.csv", row.names=FALSE)
