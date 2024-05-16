#' Lollipop chart display of KEGG/GO results
#'
#' @param data R clusterprofiler package for KEGG and GO results
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param color color color see "RColorBrewer::display.brewer.all()"
#' @param title title
#' @param text.size text size
#' @param text.width text width
#' @param ... additional parameters
#'
#' @return lollipop plot
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 margin
#' @importFrom dplyr mutate
#' @importFrom dplyr top_n
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @importFrom stats reorder
#'
#' @examples
#' \dontrun{
#' data(xfbdf)
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' eg <- bitr(unique(xfbdf$target),
#'   fromType = "SYMBOL",
#'   toType = "ENTREZID",
#'   OrgDb = "org.Hs.eg.db"
#' )
#' KK <- enrichKEGG(
#'   gene = eg$ENTREZID,
#'   organism = "hsa",
#'   pvalueCutoff = 0.05
#' )
#' MF<- enrichGO(
#'   gene = eg$ENTREZID,
#'   "org.Hs.eg.db",
#'   ont = "MF",
#'   pvalueCutoff = 0.05,
#'   readable = TRUE
#' )
#' lollipop_plot(MF, title = "Molecular Function", color = "Spectral")
#' lollipop_plot(KK)
#' }
lollipop_plot <- function(data,
                          top = 15,
                          color = "RdBu",
                          title = NULL,
                          text.size = 10,
                          text.width = 35, ...) {
  # data processing
  if (isS4(data)) {
    data <- data@result %>% tidyr::drop_na()
  } else if (is.data.frame(data)) {
    data <- data %>% tidyr::drop_na()
  } else {
    print("The data format must be S4 object or data frame.")
  }
  data2 <- dplyr::mutate(data,
    richFactor = Count / as.numeric(sub(
      "/\\d+", "",
      BgRatio
    ))
  ) %>% dplyr::slice(1:top)
  data2$richFactor <- data2$richFactor %>% round(digits = 2)
  # ggplot2 plotting
  p <- ggplot(
    data2,
    aes(richFactor, stats::reorder(Description, richFactor))
  ) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = text.size,
        colour = "black",
        vjust = 1
      ),
      axis.text.y = element_text(
        size = text.size,
        colour = "black",
        hjust = 1
      )
    ) +
    theme(axis.title = element_text(
      margin = margin(10, 5, 0, 0),
      color = "black",
      size = text.size
    )) +
    xlab("Rich Factor") +
    ylab(NULL) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width))
  return(p)
}
