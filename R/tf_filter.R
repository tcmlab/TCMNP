#' Finding transcription factors and their target genes
#'
#' @param data The data must be character as human or mouse gene symbol
#'
#' @return data.frame
#' @export
#' @examples
#' data(xfbdf, package = "TCMNP")
#' newdata <- tf_filter(xfbdf$target)
#' head(newdata)
tf_filter <- function(data) {
  data<-as.data.frame(data)
  colnames(data)<-'gene'
  if (is.character(data$gene)) {
    tf.data <- trrust[trrust$TF %in% data$gene, ] %>% as.data.frame()
    return(tf.data)
  } else {
    print("The data must be character as human or mouse gene symbol.")
  }
}
