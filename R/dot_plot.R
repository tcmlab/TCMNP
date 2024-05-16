#'  Bubble graphs showed the results of GO and KEGG analysis
#'
#' @param data R clusterprofiler package for KEGG or GO results
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param top number of categories to show
#' @param text.size text size
#' @param text.width y-axis label length
#' @param title  title
#' @param ... additional parameters
#'
#' @return dotplot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guide_colorbar
#' @importFrom dplyr mutate
#' @importFrom stats reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMNP")
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
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
#' KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
#' BP <- enrichGO(
#'   gene = eg$ENTREZID,
#'   "org.Hs.eg.db",
#'   ont = "BP",
#'   pvalueCutoff = 0.05,
#'   readable = TRUE
#' )
#' dot_plot(KK, title = "KEGG")
#' dot_plot(BP, title = "biological process")
#' }
dot_plot <- function(data,
                     color = "RdBu",
                     top = 15,
                     text.size = 10,
                     text.width = 35,
                     title = NULL, ...) {
  # Rich Factor
  data2 <- dplyr::mutate(data,
                         richFactor = Count / as.numeric(sub(
                           "/\\d+", "",
                           BgRatio
                         ))
  )
  # ggplot2 plotting
  p <- ggplot(data2,
              showCategory = top,
              aes(richFactor, stats::reorder(Description, richFactor))
  ) +
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
    ggtitle(title) +
    xlab("Rich Factor") +
    ylab(NULL) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width))
  return(p)
}
