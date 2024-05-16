#' Screen molecules based on ingredients from Chinese herbal medicine
#'
#' @param molecule The name of the molecule filter
#' @param ... additional parameters
#'
#' @return data.frame
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' molecule <- c("Hexyl acetate", "quercetin 3-o-rhamnopyranosyl")
#' molecule_filter<- molecule_serch(molecule)
#' head(molecule_filter)
molecule_serch <- function(molecule, ...){
    if (all(molecule %in% (unique(tcmsp$Molecule_name))) == FALSE) {
      molecule_deficiency <- setdiff(molecule, unique(tcmsp$Molecule_name))
      print(paste0(molecule_deficiency, " is/are not in the datasets."))
      molecule <- molecule[-match(molecule_deficiency, molecule)]
    }
  molecule2 <- tcmsp[which((tcmsp$Molecule_name) %in% molecule), ] %>%
      dplyr::select(c(3:5, 9)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(molecule2) <- c("herb", "molecule_id", "molecule", "target")
    return(molecule2)
}
