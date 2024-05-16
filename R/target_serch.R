#' Screen target based on ingredients from Chinese herbal medicine
#'
#' @param target The name of the target filter
#' @param ... additional parameters
#'
#' @return data.frame
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' target <- c("PRSS3", "ADH1C")
#' target_filter<- target_serch(target, type = "Herb_name_pin_yin")
#' head(target_filter)
target_serch <- function(target, ...){
    if (all(target %in% (unique(tcmsp$Target))) == FALSE) {
      target_deficiency <- setdiff(target, unique(tcmsp$Target))
      print(paste0(target_deficiency, " is/are not in the datasets."))
      target <- target[-match(target_deficiency, target)]
    }
    fufang <- tcmsp[which((tcmsp$Target) %in% target), ] %>%
      dplyr::select(c(3:5, 9)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(fufang) <- c("herb", "molecule_id", "molecule", "target")
    return(fufang)
}

