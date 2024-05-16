#' Intersection of datasets
#'
#' @param data data.frame
#' column names must be two columns of data for gene and source
#'
#' @return data.frame
#'
#' @importFrom dplyr filter
#' @importFrom methods new
#' @importFrom VennDetail make.subset
#' @importFrom VennDetail venndetail
#' @importFrom VennDetail result
#' @export
#'
#' @examples
#'  \dontrun{
#' data(venn_data, package = "TCMNP")
#' venn_result(venn_data)
#' }
#'
venn_result <- function(data) {
  # data processing
  data <- as.data.frame(data)
  if(is.data.frame(data)) {
    data <- data %>% dplyr::distinct()
    if (is.data.frame(data)) {
      if (length(unique(data[, 1])) > length(unique(data[, 2]))) {
        colnames(data) <- c("gene", "source")
      } else {
        colnames(data) <- c("source", "gene")
      }
    } else {
      print("The data must be a data frame with two columns.")
    }
  }
  database <- unique(data$source)
  df <- list()
  df2 <- list()
  for (i in database) {
      df[[i]] <- data %>% dplyr::filter(source == i)
      df2[[i]] <- df[[i]]$gene
    }
names(df2) <- database

results <- VennDetail::venndetail(df2) %>%
          VennDetail::result(wide = TRUE) %>%
          as.data.frame()
rownames(results) <- 1:nrow(results)
return(results)
}

