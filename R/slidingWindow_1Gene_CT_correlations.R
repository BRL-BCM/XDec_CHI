#' Identify all genes that correlate to a gene of interest in a cell type of interest
#'
#' This function takes bulk gene expression and esimated cell type proportions
#' and performs sliding window deconvolution. The output is all genes that correlate.
#' to a gene of interest in a specific cell type. It assumes that both input
#' matrixes are numeric with features (expression) or samples (proportion) as row names.
#' Columns of the expression matrix are samples and rows of the proportion matrix are
#' samples. Only samples present in both matrixes will be used. There must be more
#' than 40 samples total.
#'
#' @param GeneExp numeric matrix with samples as columns and features as rows
#' @param CTproportions numeric matrix with samples as rows and cell type names as columns
#' @param GOI_anchor gene to order samples by
#' @param CT_anchor cell type of anchor gene
#' @param CT_cor cell type where you want to identify correlated genes
#' @param Type_pValue p.value cut off FDR or p.value
#' @param Significance_cutoff significance must be less than this to be reports
#' @param Correlation_cutoff correlation must be less than this to be reports
#' @return reports p.value, FDR, and correlation for genes passing thresholds in a given cell type
#' @export
slidingWindow_1Gene_CT_correlations <- function(GeneExp,
                                                CTproportions,
                                                GOI_anchor,
                                                CT_anchor,
                                                CT_cor,
                                                Type_pValue="p.value",
                                                Significance_cutoff=0.05,
                                                Correlation_cutoff=0.5){

  CT_anchor=gsub(" ",".",CT_anchor)
  CT_cor=gsub(" ",".",CT_anchor)
  colnames(CTproportions)=gsub(" ",".",colnames(CTproportions))


  CTproportions=CTproportions[colnames(GeneExp),]
  GeneExp=GeneExp[,rownames(CTproportions)]

  stopifnot(length(colnames(GeneExp))>40)

  stopifnot(GOI_anchor %in% rownames(GeneExp))
  stopifnot(CT_anchor %in% colnames(CTproportions))
  stopifnot(CT_cor %in% colnames(CTproportions))

  stopifnot(Type_pValue=="p.value" | Type_pValue=="FDR")
  stopifnot(Significance_cutoff > 0 && Significance_cutoff < 1)
  stopifnot(Correlation_cutoff > -1 && Correlation_cutoff < 1)

  get_p_val_gene_2CT <- function(...,list_CT,list_CT2) {
    CT_df = as.numeric(list_CT)
    CT_df2 = as.matrix(list_CT2)
    p_Val = c()
    correlation = c()
    #For every gene in second cell type, correlated GOI in first cell type
    for (gene in 1:dim(CT_df2)[1]) {
      if (sd(CT_df2[gene,])==0) {
        p_Val=c(p_Val,1)
        correlation = c(correlation,suppressWarnings(cor(CT_df,CT_df2[gene,])))
      } else {
        #Use pearson correlation
        p_Val=c(p_Val,suppressWarnings(cor.test(CT_df,CT_df2[gene,]))$p.value)
        correlation = c(correlation,cor(CT_df,CT_df2[gene,]))
      }
    }
    p_val_adjusted = p.adjust(p_Val, method = "bonferroni")
    all_info = data.frame(genes = rownames(CT_df2),
                          p.val = p_Val,
                          FDR = p_val_adjusted,
                          Correlation = correlation)
    return(all_info)
  }

  GOI_order = GeneExp[,order(GeneExp[GOI_anchor,])]
  order_cancer = colnames(GOI_order)
  num_use=length(order_cancer)-39

  for (j in 1:num_use) {
    start = j
    end=j+39
    assign(paste("group_",j,sep=""),order_cancer[start:end])
  }

  ##2 myeloid, 2 cancer
  input.Stage2 = CTproportions[order_cancer,]

  props.names = list()
  props.matrix = list()
  Expression.matrix = list()
  Stage2.matrix = list()

  means_names = colnames(CTproportions)

  CT1_get_means = c()
  CT2_get_means = data.frame(rep(0,dim(GOI_order)[1]))

  #########
  ###2 cancer
  ########
  for(i in 1:num_use){
    # for(i in 1:3){
    props.names = rownames(input.Stage2[rownames(input.Stage2) %in% get(paste("group_",i,sep="")),])
    props.matrix = as.matrix(input.Stage2[props.names,])
    Expression.matrix = as.matrix(GOI_order[,props.names])
    ##Run Stage 2 on each subset
    Stage2.matrix = EDec::run_edec_stage_2(gene_exp_bulk_samples = Expression.matrix, cell_type_props = props.matrix)

    CT1_get_means = c(CT1_get_means,Stage2.matrix$means[GOI_anchor,which(means_names==CT_anchor)])
    CT2_get_means = cbind.data.frame(CT2_get_means,Stage2.matrix$means[,which(means_names==CT_cor)])

    print(paste("performing",i,"of",num_use,"windows"))
  }

  CT2_get_means=CT2_get_means[,-1]


  print("performing correlations")
  GOI_CT_cor = get_p_val_gene_2CT(list_CT = CT1_get_means,
                                  list_CT2 = CT2_get_means)


  if (Type_pValue=="p.value") {
    GOI_CT_cor_sig = GOI_CT_cor[GOI_CT_cor[,2] < Significance_cutoff,]
  } else {
    GOI_CT_cor_sig = GOI_CT_cor[GOI_CT_cor[,3] < Significance_cutoff,]
  }
  GOI_CT_cor_sig_cor = GOI_CT_cor_sig[abs(GOI_CT_cor_sig[,4]) > Correlation_cutoff,]


  Cor_table <- GOI_CT_cor_sig_cor[order(GOI_CT_cor_sig_cor[,4],decreasing = T),]
  Cor_table_up <- GOI_CT_cor_sig_cor[GOI_CT_cor_sig_cor[,4]>Correlation_cutoff,]
  Cor_table_up=Cor_table_up[order(Cor_table_up[,4],decreasing = T),]
  Cor_table_down <- GOI_CT_cor_sig_cor[GOI_CT_cor_sig_cor[,4]<Correlation_cutoff,]
  Cor_table_down=Cor_table_down[order(Cor_table_down[,4],decreasing = F),]

  my_list = list("all_correlations"=Cor_table,
                 "Positive_correlations"=Cor_table_up,
                 "Negative_correlations"=Cor_table_down)

  return(my_list)
}
