#' Dot plot showed the results of GO results
#'
#' @param go.diff R clusterprofiler package for enrichGO results
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param text.size text size
#' @param text.width text width
#' @param title title
#' @param ... additional parameters
#'
#' @return dotplot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 facet_grid
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom tidyr drop_na
#' @importFrom stringr str_wrap
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
#' data(xfbdf, package = "TCMNP")
#' eg <- bitr(unique(xfbdf$target),
#'   fromType = "SYMBOL",
#'   toType = "ENTREZID",
#'   OrgDb = "org.Hs.eg.db"
#' )
#' go.diff <- enrichGO(
#'   gene = eg$ENTREZID,
#'   org.Hs.eg.db,
#'   pAdjustMethod = "BH",
#'   pvalueCutoff = 0.01,
#'   qvalueCutoff = 0.05,
#'   ont = "all",
#'   readable = TRUE
#' )
#' go_dotplot(go.diff, color = "Spectral")
#' }
go_dotplot <- function(go.diff,
                       top = 5,
                       color = "RdBu",
                       text.size = 12,
                       text.width = 35,
                       title = "GO enrichment of cluster", ...) {
  # data processing
  if (isS4(go.diff)) {
    go.diff <- go.diff@result
  } else if (is.data.frame(go.diff)) {
    go.diff <- go.diff
  } else {
    print("The data format must be S4 object or data frame.")
  }

  go_enrich <- go.diff %>%
    tidyr::drop_na() %>%
    dplyr::group_by(ONTOLOGY) %>%
    top_n(n = top) %>%
    dplyr::mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  go_enrich <- go_enrich[order(go_enrich$richFactor,
    decreasing = FALSE
  ), ]
  go_enrich$Description <- factor(go_enrich$Description,
    levels = go_enrich$Description
  )
  # ggplot2 plotting
  p <- ggplot(go_enrich, aes(richFactor, Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 8)) +
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
    xlab("Rich Factor") +
    ylab(NULL) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width)) +
    theme(strip.text = element_text(face = "bold", size = text.size))
  return(p)
}
