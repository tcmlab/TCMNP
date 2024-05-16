#' Screen herbs and compound based on disease gene
#'
#' @param data disease gene
#' @param OB.value oral bioavailability
#' @param DL.value drug-likeness
#' @param compound.source Writing format of herb name:
#' Chinese medicine prescription(tcmp) and Chinese patent medicine(cpm)
#' "tcmp" or "pcm"
#' @param freq.value herbs freqence
#' @param digits.value Keep the number of decimal places
#' @param ... additional parameters
#'
#' @return data.frame
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' data(disease_data, package = "TCMNP")
#' herbs <- tcm_prescription(disease_data)[[1]]
#' compound <- tcm_prescription(disease_data)[[2]]

tcm_prescription <- function(data,
                             compound.source = "tcmp",
                             OB.value = 30,
                             DL.value = 0.18,
                             freq.value = 50,
                             digits.value = 2, ...) {
  # Search for the right herbal medicine
  data <- as.data.frame(data)
  colnames(data)<- 'disease_gene'

  tcmsp2 <- tcmsp[tcmsp$Target %in% data$disease_gene, ] %>%
    tidyr::drop_na() %>%
    dplyr::distinct(., .keep_all = TRUE) %>%
    dplyr::filter(., OB >= OB.value & DL >= DL.value)

  # Count the frequency of occurrence of herbs
  herb_frq <- tcmsp2 %>%
    dplyr::select(Herb_cn_name) %>%
    dplyr::group_by(Herb_cn_name) %>%
    dplyr::mutate(freq = dplyr::n()) %>%
    arrange(desc(freq)) %>%
    dplyr::distinct() %>%
    dplyr::filter(freq > freq.value)

  # Construct a data set of traditional Chinese medicine prescriptions
  tmp <- compoundfinal %>%
    dplyr::filter(type == compound.source) %>%
    dplyr::select(herb, compound) %>%
    as.data.frame()

  # Hypergeometric distribution calculation
  # 背景中药数目(m），就是所有复方注释到的所有草药的数目
  # The number of background Chinese medicines (m) is the
  # number of all herbal medicines annotated in all compounds.

  # 特定复方上的草药数目（n）
  # Number of herbs on a specific recipe (n)

  # 目的复方中草药与背景草药m重合的草药数目（k）
  # The number of herbs (k) that overlap the purpose compound
  # Chinese herbal medicine and the background herbal medicine m

  # k中落在特定复方中的草药数（a）
  # The number of herbs in k that fall into a specific compound (a)
  herb_compound <- tapply(tmp[, 2], as.factor(tmp[, 1]), function(x) x)
  compound_herb <- tapply(tmp[, 1], as.factor(tmp[, 2]), function(x) x)
  herb_name <- herb_frq$Herb_cn_name
  herb_name_has_compound <- intersect(herb_name, names(herb_compound))
  n <- length(herb_name)
  universeherb <- unique(tcmsp$Herb_cn_name)
  N <- length(universeherb)

  results <- c()
  for (i in names(compound_herb)) {
    M <- length(intersect(compound_herb[[i]], universeherb))
    if (M < 5) next
    exp_count <- n * M / N
    k <- 0
    for (j in herb_name_has_compound) {
      if (i %in% herb_compound[[j]]) k <- k + 1
    }
    Oddratio <- k / exp_count
    if (k == 0) next
    p <- phyper(k - 1, M, N - M, n, lower.tail = F)
    results <- rbind(results, c(i, p, Oddratio, exp_count, k, M))
  }
  results <- as.data.frame(results, stringsAsFactors = F)
  colnames(results) <- c("CompoundId", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size")
  results$p.adjust <- stats::p.adjust(results$Pvalue, method = "BH") %>%
    base::format(digits = digits.value, nsmall = digits.value) %>%
    as.numeric()
  results <- results[order(results$Count, decreasing = TRUE), ]
  rownames(results) <- 1:nrow(results)
  results$Pvalue <- as.numeric(results$Pvalue) %>%
    base::format(digits = digits.value, nsmall = digits.value) %>%
    as.numeric()
  results$OddsRatio <- as.numeric(results$OddsRatio) %>%
    base::format(digits = digits.value, nsmall = digits.value) %>%
    as.numeric()
  results$ExpCount <- as.numeric(results$ExpCount) %>%
    base::format(digits = digits.value, nsmall = digits.value) %>%
    as.numeric()
  results$Count <- as.numeric(results$Count)
  results$Size <- as.numeric(results$Size)
  # Add the intersection of compound herbal medicine and background herbal medicine
  herb_list <- list()
  for (i in seq(dplyr::row_number(results))) {
    herb_list[[i]] <- intersect(tmp[tmp$compound %in% results$CompoundId[i], ]$herb, herb_name) %>%
      paste0(sep = "", collapse = "/")
  }
  results$herb <- unlist(herb_list)
  return(list(herb_frq, results))
}

