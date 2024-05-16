#' Find Chinese medicines and their active compounds through the target
#'
#' @param gene a vector of gene name
#' @param OB.value oral bioavailability
#' @param DL.value drug-likeness
#' @param ... additional parameters
#'
#' @return data.frame
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' gene <- c("MAPK1", "JUN", "FOS", "RAC1", "IL1", "IL6")
#' target_herbs <- target_herb(gene)
target_herb <- function(gene, OB.value = 30, DL.value = 0.18, ...) {
  if (all(gene %in% (unique(tcmsp$Target))) == FALSE) {
    gene_deficiency <- setdiff(gene, unique(tcmsp$Target))
    print(paste0(gene_deficiency, " is/are not in the datasets."))
    gene <- gene[-match(gene_deficiency, gene)]
  }
  data_herbs <- tcmsp[which(tcmsp$Target %in% gene), ] %>%
    dplyr::filter(OB >= OB.value & DL >= DL.value) %>%
    dplyr::select(c(3:5, 9)) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()
  colnames(data_herbs) <- c("herb", "molecule_id", "molecule", "target")
  return(data_herbs)
}
