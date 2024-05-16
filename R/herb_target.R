#' Screen targets based on ingredients from Chinese herbal medicine
#'
#' @param herb The name of the compound Chinese herbal medicine
#' @param OB.value oral bioavailability
#' @param DL.value drug-likeness
#' @param type Writing format of herb name:
#' "Herb_name_pin_yin" or "Herb_cn_name"
#' @param ... additional parameters
#'
#' @return data.frame
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' herbs <- c("ma huang", "ku xing ren")
#' fufang <- herb_target(herbs, type = "Herb_name_pin_yin")
#' head(fufang)
herb_target <- function(herb,
                        OB.value = 30,
                        DL.value = 0.18,
                        type = "Herb_cn_name",
                        ...) {
  if (type == "Herb_cn_name") {
    if (all(herb %in% (unique(tcmsp$Herb_cn_name))) == FALSE) {
      herb_deficiency <- setdiff(herb, unique(tcmsp$Herb_cn_name))
      print(paste0(herb_deficiency, " is/are not in the datasets."))
      herb <- herb[-match(herb_deficiency, herb)]
    }
    fufang <- tcmsp[which((tcmsp$Herb_cn_name) %in% herb), ] %>%
      dplyr::filter(OB >= OB.value & DL >= DL.value) %>%
      dplyr::select(c(3:5, 9)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(fufang) <- c("herb", "molecule_id", "molecule", "target")
    return(fufang)
  } else if (type == "Herb_name_pin_yin") {
    if (all(herb %in% (unique(tcmsp$Herb_name_pin_yin))) == FALSE) {
      herb_deficiency <- setdiff(herb, unique(tcmsp$Herb_name_pin_yin))
      print(paste0(herb_deficiency, " is/are not in the datasets."))
      herb <- herb[-match(herb_deficiency, herb)]
    }
    fufang <- tcmsp[which((tcmsp$Herb_name_pin_yin) %in% herb), ] %>%
      dplyr::filter(OB >= OB.value & DL >= DL.value) %>%
      dplyr::select(c(3:5, 9)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(fufang) <- c("herb", "molecule_id", "molecule", "target")
    return(fufang)
  }
}
