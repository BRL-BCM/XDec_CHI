#' Correlate 2 genes in a cell type specific manner
#'
#' This function takes bulk gene expression and esimated cell type proportions
#' and performs sliding window deconvolution. The output is correlation of 2
#' features of interest in a cell type specific manner. It assumes that both input
#' matrixes are numeric with features (expression) or samples (proportion) as row names.
#' Columns of the expression matrix are samples and rows of the proportion matrix are
#' samples. Only samples present in both matrixes will be used. There must be more
#' than 40 samples total.
#'
#' @param GeneExp numeric matrix with samples as columns and features as rows
#' @param CTproportions numeric matrix with samples as rows and cell type names as columns
#' @param GOI_anchor gene to order samples by
#' @param GOI_cor gene you want to correlate to
#' @param CT_anchor cell type of anchor gene
#' @param CT_cor cell type of gene you want to correlate to
#' @return Correlation, p.value, and figures of bulk and cell type specific expression of genes
#'
#' @import EDec
#' @import ggplot2
#' @import ggpubr
#'
#' @export
slidingWindow_2Gene_2CT <- function(GeneExp,
                                    CTproportions,
                                    GOI_anchor,
                                    GOI_cor,
                                    CT_anchor,
                                    CT_cor){

  stopifnot(GOI_anchor %in% rownames(GeneExp))
  stopifnot(GOI_cor %in% rownames(GeneExp))

  CT_anchor=gsub(" ",".",CT_anchor)
  CT_cor=gsub(" ",".",CT_anchor)
  colnames(CTproportions)=gsub(" ",".",colnames(CTproportions))

  stopifnot(CT_anchor %in% colnames(CTproportions))
  stopifnot(CT_cor %in% colnames(CTproportions))

  CTproportions=CTproportions[colnames(GeneExp),]
  GeneExp=GeneExp[,rownames(CTproportions)]

  stopifnot(length(colnames(GeneExp))>40)


  GOI_order = GeneExp[c(GOI_anchor,GOI_cor),]
  GOI_order = GOI_order[,order(GOI_order[GOI_anchor,])]
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
  CT2_get_means = c()

  CT1_get_SD = c()
  CT2_get_SD = c()

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
    CT2_get_means = c(CT2_get_means,Stage2.matrix$means[GOI_cor,which(means_names==CT_cor)])

    CT1_get_SD = c(CT1_get_SD,Stage2.matrix$std.errors[GOI_anchor,which(means_names==CT_anchor)])
    CT2_get_SD = c(CT2_get_SD,Stage2.matrix$std.errors[GOI_cor,which(means_names==CT_cor)])

    print(paste("performing",i,"of",num_use,"windows"))
  }


  #Correlation
  cor_return = cor(as.numeric(CT1_get_means),as.numeric(CT2_get_means),method = "pearson")
  p_return = cor.test(as.numeric(CT1_get_means),as.numeric(CT2_get_means),method = "pearson")
  #Correlation plot
  Cor_plot = cbind.data.frame(as.numeric(CT1_get_means),
                              as.numeric(CT2_get_means),
                              1:num_use)
  colnames(Cor_plot) = c("CT1","CT2","Order")
  p1=ggplot2::ggplot(Cor_plot, aes(x=CT1, y=CT2)) + geom_point(aes(colour=Order)) +
    xlab(paste(CT_anchor,GOI_anchor)) +
    ylab(paste(CT_cor,GOI_cor)) +
    ggpubr::stat_cor(method = "pearson", label.x = mean(Cor_plot$CT1), label.y = (max(Cor_plot$CT2)+5)) +
    scale_colour_viridis_c()

  #Expression plot 1st gene
  Means_EDec = data.frame(means=CT1_get_means,
                          variable=1:length(CT1_get_means))
  EDec_std = data.frame(std.error=CT1_get_SD,
                        variable=1:length(CT1_get_means))

  x_mean_std = merge(Means_EDec,EDec_std,by = "variable")
  x_mean_std$variable=as.factor(x_mean_std$variable)

  p2 = ggplot2::ggplot(x_mean_std,aes(x = variable , y = means, fill=variable)) +
    geom_bar(stat="identity", position = "dodge")+
    scale_fill_viridis_d()+
    geom_errorbar(aes(ymax= means + std.error,
                      ymin=ifelse(means-std.error < 0,0,means-std.error)),
                  position = position_dodge(0.95), width = 0.25) +
    ggtitle(paste(CT_anchor,GOI_anchor)) +
    theme(plot.title = element_text(size=10,hjust = 0.5, face="bold")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=5))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) +
    theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    xlab("Profiles") + ylab("Predicted Expression") +
    theme(legend.position = "none")

  #Expression plot 2nd gene
  Means_EDec = data.frame(means=CT2_get_means,
                          variable=1:length(CT2_get_means))
  EDec_std = data.frame(std.error=CT2_get_SD,
                        variable=1:length(CT2_get_means))



  x_mean_std = merge(Means_EDec,EDec_std,by = "variable")
  x_mean_std$variable=as.factor(x_mean_std$variable)

  p3 = ggplot2::ggplot(x_mean_std,aes(x = variable , y = means, fill=variable)) +
    geom_bar(stat="identity", position = "dodge")+
    scale_fill_viridis_d()+
    geom_errorbar(aes(ymax= means + std.error,
                      ymin=ifelse(means-std.error < 0,0,means-std.error)),
                  position = position_dodge(0.95), width = 0.25) +
    ggtitle(paste(CT_cor,GOI_cor)) +
    theme(plot.title = element_text(size=10,hjust = 0.5, face="bold")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=5))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) +
    theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    xlab("Profiles") + ylab("Predicted Expression") +
    theme(legend.position = "none")


  Means_EDec = data.frame(means=as.numeric(GOI_order[GOI_anchor,]),
                          variable=1:length(GOI_order[GOI_anchor,]))
  EDec_std = data.frame(std.error=rep(sd(GOI_order[GOI_anchor,]),length(GOI_order[GOI_anchor,])),
                        variable=1:length(GOI_order[GOI_anchor,]))

  x_mean_std = merge(Means_EDec,EDec_std,by = "variable")
  x_mean_std$variable=as.factor(x_mean_std$variable)

  p4 = ggplot2::ggplot(x_mean_std,aes(x = variable , y = means, fill=variable)) +
    geom_bar(stat="identity", position = "dodge")+
    scale_fill_viridis_d()+
    geom_errorbar(aes(ymax= means + std.error,
                      ymin=ifelse(means-std.error < 0,0,means-std.error)),
                  position = position_dodge(0.95), width = 0.25) +
    ggtitle(paste("Bulk",GOI_anchor)) +
    theme(plot.title = element_text(size=10,hjust = 0.5, face="bold")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=5))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) +
    theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    xlab("Samples") + ylab("Expression") +
    theme(legend.position = "none")


  Means_EDec = data.frame(means=as.numeric(GOI_order[GOI_cor,]),
                          variable=1:length(GOI_order[GOI_cor,]))
  EDec_std = data.frame(std.error=rep(sd(GOI_order[GOI_cor,]),length(GOI_order[GOI_cor,])),
                        variable=1:length(GOI_order[GOI_cor,]))

  x_mean_std = merge(Means_EDec,EDec_std,by = "variable")
  x_mean_std$variable=as.factor(x_mean_std$variable)

  p5 = ggplot2::ggplot(x_mean_std,aes(x = variable , y = means, fill=variable)) +
    geom_bar(stat="identity", position = "dodge")+
    scale_fill_viridis_d()+
    geom_errorbar(aes(ymax= means + std.error,
                      ymin=ifelse(means-std.error < 0,0,means-std.error)),
                  position = position_dodge(0.95), width = 0.25) +
    ggtitle(paste("Bulk",GOI_cor)) +
    theme(plot.title = element_text(size=10,hjust = 0.5, face="bold")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=5))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) +
    theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    xlab("Samples") + ylab("Expression") +
    theme(legend.position = "none")

  my_list <- list("Correlation" = cor_return,
                  "p.value" = p_return,
                  "Cor_plot" = p1,
                  "CT_anchor_expression" = p2,
                  "CT_cor_expression" = p3,
                  "Bulk_anchor_expression" = p4,
                  "Bulk_cor_expression" = p5)
  return(my_list)
}
