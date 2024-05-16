#' sankey and dot graphs showed the results of GO and KEGG analysis
#'
#' @param data the results after KEGG/GO analysis
#' @param top according to the order of p adjust value from small to large
#' the number of categories to show
#' @param sankey.text.size Sankey diagram text size
#' @param sankey.x.axis.text.position Location of x-axis labels in sankey plot
#' @param sankey.x.axis.text.size Sankey diagram x-axis text size
#' @param dot.color dotplot color: see "RColorBrewer::display.brewer.all()"
#' @param dot.position dotplot position
#' @param dot.lable "ID" or "Description"
#' @param dot.text.size dotplot text size
#' @param dot.scale dotplot scale
#' @param dot.x dotplot horizontal position
#' @param dot.y dotplot up and down position
#' @param dot.width dotplot width
#' @param dot.height dotplot height
#' @param ... additional parameters
#'
#' @return figure
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_distiller
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom dplyr slice
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom ggsankey make_long
#' @importFrom ggsankey geom_sankey
#' @importFrom ggsankey geom_sankey_text
#' @importFrom ggsankey theme_sankey
#' @importFrom cols4all c4a
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @importFrom stats reorder
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMNP")
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
#' library(tidyverse)
#' eg <- bitr(unique(xfbdf$target),
#'   fromType = "SYMBOL",
#'   toType = "ENTREZID",
#'   OrgDb = "org.Hs.eg.db"
#' )
#' KK <- enrichKEGG(
#'   eg$ENTREZID,
#'   organism = "hsa",
#'   pvalueCutoff = 0.05
#' )
#' KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
#' dot_sankey(KK)
#' dot_sankey(kegg.filter)
#' }
dot_sankey <- function(data,
                       top = 8,
                       sankey.text.size = 3,
                       sankey.x.axis.text.position = 1,
                       sankey.x.axis.text.size = 12,
                       dot.color = "RdBu",
                       dot.position = 8,
                       dot.lable = "Description",
                       dot.text.size = 12,
                       dot.scale = 0.9,
                       dot.x = 0.45,
                       dot.y = 0.17,
                       dot.width = 0.5,
                       dot.height = 0.63,
                       ...) {
  # data processing
  if (isS4(data)) {
    data <- data@result %>% tidyr::drop_na()
  } else if (is.data.frame(data)) {
    data <- data %>% tidyr::drop_na()
  } else {
    print("The data format must be S4 object or data frame.")
  }
  data_dot <- data %>%
    dplyr::slice(1:top) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
    dplyr::arrange(richFactor, descreasing = TRUE) %>%
    dplyr::select("ID", "Description", "p.adjust", "Count", "richFactor")
  data_dot$richFactor <- data_dot$richFactor %>% round(digits = 2)
  # Start drawing Bubble Chart
  if (dot.lable == "Description") {
    p1 <- ggplot(
      data = data_dot,
      aes(
        x = richFactor,
        y = reorder(Description, richFactor)
      )
    ) +
      geom_point(aes(size = Count, color = -log10(p.adjust))) +
      scale_size_continuous(range = c(2, 6)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 0, hjust = 0.5, vjust = 0.5,
          size = dot.text.size, colour = "black"
        ),
        axis.text.y = element_text(
          angle = 0, size = dot.text.size,
          face = "plain", colour = "black"
        )
      ) +
      theme(axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = dot.text.size
      )) +
      scale_colour_distiller(palette = dot.color, direction = -1) +
      labs(x = "richFactor", y = "") +
      theme(
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )
  } else if (dot.lable == "ID") {
    p1 <- ggplot(
      data = data_dot,
      aes(
        x = richFactor,
        y = reorder(Description, richFactor)
      )
    ) +
      geom_point(aes(size = Count, color = -log10(p.adjust))) +
      scale_size_continuous(range = c(2, 6)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = dot.text.size,
          colour = "black"
        ),
        axis.text.y = element_text(
          angle = 0,
          size = dot.text.size,
          face = "plain",
          colour = "black"
        )
      ) +
      theme(axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = dot.text.size
      )) +
      scale_colour_distiller(palette = dot.color, direction = -1) +
      labs(x = "richFactor", y = "") +
      theme(
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )
  }

  data_sankey <- data %>%
    dplyr::slice(1:top) %>%
    dplyr::mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
    dplyr::arrange(richFactor, descreasing = TRUE) %>%
    dplyr::select(geneID, dot.lable) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    as.data.frame() %>%
    dplyr::distinct()

  df <- data_sankey %>% ggsankey::make_long(colnames(data_sankey))
  if (dot.lable == "Description") {
    df$node <- factor(df$node, levels = c(
      data_sankey$Description %>% unique(),
      data_sankey$geneID %>% unique()
    ))
  } else if (dot.lable == "ID") {
    df$node <- factor(df$node, levels = c(
      data_sankey$ID %>% unique(),
      data_sankey$geneID %>% unique()
    ))
  } else {
    (print('The dot.lable must be "Description" or "ID".'))
  }

  mycol <- cols4all::c4a("rainbow_wh_rd", length(unique(df$node)))
  mycol2 <- sample(mycol, length(mycol), replace = FALSE)

  p2 <- ggplot(df, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = node
  )) +
    ggsankey::geom_sankey(
      flow.alpha = 0.5,
      node.color = NA,
      show.legend = FALSE
    ) +
    ggsankey::geom_sankey_text(
      size = sankey.text.size,
      color = "black",
      hjust = 1,
      position = position_nudge(x = -0.05)
    ) +
    ggsankey::theme_sankey(base_size = 20) +
    scale_fill_manual(values = mycol2) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(
      size = sankey.x.axis.text.size,
      hjust = sankey.x.axis.text.position,
      vjust = 10, colour = "black"
    ))
  p3 <- p2 + theme(plot.margin = unit(c(0, dot.position, 0, 0),
    units = "cm"
  ))
  p4 <- cowplot::ggdraw() + cowplot::draw_plot(p3) +
    cowplot::draw_plot(p1,
      scale = dot.scale,
      x = dot.x,
      y = dot.y,
      width = dot.width,
      height = dot.height
    )
  return(p4)
}
