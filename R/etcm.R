#' ETCM database download data collation
#'
#' @param data data downloaded from ETCM database
#' @param herb herb herb name
#'
#' @return data.frame
#' @export
#'
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom tidyr drop_na
#' @importFrom stringr str_extract
#' @importFrom stringr str_trim
#' @importFrom utils stack
#' @examples
#' data("mahuang", package = "TCMNP")
#' data <- etcm(mahuang, herb = "ma huang")
#' head(data)
etcm <- function(data, herb = NULL) {
  data <- data %>% tidyr::drop_na()
  # data cleaning
  if ("Chemical Component" %in% colnames(data)) {
    df2 <- lapply(
      data$`Candidate Target genes`,
      function(x) base::strsplit(x, ",")[[1]]
    )
    names(df2) <- data$`Chemical Component`
    df3 <- utils::stack(df2) %>%
      dplyr::rename(target = "values", molecule = "ind")
    # add new column
    df4 <- df3 %>%
      dplyr::mutate(
        QED = stringr::str_extract(df3$target, "(?<=\\().+?(?=\\))") %>%
          stringr::str_trim(),
        target = stringr::str_extract(df3$target, "^.*(?=\\()") %>%
          stringr::str_trim(),
        herb = herb
      ) %>%
      dplyr::select(herb, molecule, target, QED)
    return(df4)
  } else if ("Chemical.Component" %in% colnames(data)) {
    df2 <- lapply(
      data$Candidate.Target.genes,
      function(x) base::strsplit(x, ",")[[1]]
    )
    names(df2) <- data$Chemical.Component
    df3 <- utils::stack(df2) %>%
      dplyr::rename(target = "values", molecule = "ind")
    # add new column
    df4 <- df3 %>%
      dplyr::mutate(
        QED = stringr::str_extract(df3$target, "(?<=\\().+?(?=\\))") %>%
          stringr::str_trim(),
        target = stringr::str_extract(df3$target, "^.*(?=\\()") %>%
          stringr::str_trim(),
        herb = herb
      ) %>%
      dplyr::select(herb, molecule, target, QED)
    return(df4)
  }
}
