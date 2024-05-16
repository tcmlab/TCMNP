#' Herb-ingredient-target-function enrichment analysis
#'
#' @param data data.frame
#' @param sankey.text.size Sankey diagram text size
#' @param sankey.x.axis.text.position Location of x-axis labels in sankey plot
#' @param sankey.x.axis.text.size Sankey diagram x-axis text size
#' @param molecule.lable "molecule" or ""molecule_id"
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
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_colour_distiller
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom tidyr drop_na
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
#' KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID") %>%
#'   dplyr::mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#' path <- separate_rows(KK@result, geneID, sep = "/")
#' data_sankey <- left_join(xfbdf, path,
#'   by = c("target" = "geneID"),
#'   relationship = "many-to-many"
#' ) %>%
#'   dplyr::distinct() %>%
#'   tidyr::drop_na() %>%
#'   dplyr::sample_n(30, replace = FALSE) %>%
#'   as.data.frame()
#' tcm_sankey_dot(data_sankey)
#' }
tcm_sankey_dot <- function(data,
                           sankey.text.size = 3,
                           sankey.x.axis.text.position = 1,
                           sankey.x.axis.text.size = 12,
                           molecule.lable = "molecule_id",
                           dot.color = "RdBu",
                           dot.position = 6,
                           dot.lable = "Description",
                           dot.text.size = 12,
                           dot.scale = 0.72,
                           dot.x = 0.59,
                           dot.y = -0.173,
                           dot.width = 0.48,
                           dot.height = 1.33,
                           ...) {
  # Data collation after KEGG or GO enrichment analysis
  data3 <- data %>%
    dplyr::select("ID", "Description", "p.adjust", "Count", "richFactor") %>%
    dplyr::distinct() %>%
    dplyr::arrange(p.adjust, descreasing = TRUE)
  data3$richFactor <- data3$richFactor %>% round(digits = 2)

  # bubble chart
  if (dot.lable == "Description") {
    p1 <- ggplot(
      data = data3,
      aes(
        x = richFactor,
        y = stats::reorder(Description, richFactor)
      )
    ) +
      geom_point(aes(size = Count, color = -log10(p.adjust))) +
      scale_size_continuous(range = c(2, 6)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          size = dot.text.size,
          colour = "black", vjust = 1
        ),
        axis.text.y = element_text(
          size = dot.text.size,
          colour = "black", hjust = 1
        )
      ) +
      scale_colour_distiller(palette = dot.color, direction = -1) +
      labs(x = "richFactor", y = "") +
      theme_bw() +
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
      data = data3,
      aes(
        x = richFactor,
        y = stats::reorder(ID, richFactor)
      )
    ) +
      geom_point(aes(size = Count, color = -log10(p.adjust))) +
      scale_size_continuous(range = c(2, 6)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          size = dot.text.size,
          colour = "black", vjust = 1
        ),
        axis.text.y = element_text(
          size = dot.text.size,
          colour = "black", hjust = 1
        )
      ) +
      scale_colour_distiller(palette = dot.color, direction = -1) +
      labs(x = "richFactor", y = "") +
      theme_bw() +
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
    dplyr::arrange(richFactor, descreasing = FALSE) %>%
    dplyr::select(herb, molecule_id, target, dot.lable) %>%
    as.data.frame() %>%
    dplyr::distinct()

  df <- data_sankey %>% ggsankey::make_long(colnames(data_sankey))

  if (dot.lable == "Description" && molecule.lable == "molecule_id") {
    df$node <- factor(df$node, levels = c(
      data_sankey$Description %>% unique(),
      data_sankey$target %>% unique(),
      data_sankey$molecule_id %>% unique(),
      data_sankey$herb %>% unique()
    ))
  } else if (dot.lable == "ID" && molecule.lable == "molecule_id") {
    df$node <- factor(df$node, levels = c(
      data_sankey$ID %>% unique(),
      data_sankey$target %>% unique(),
      data_sankey$molecule_id %>% unique(),
      data_sankey$herb %>% unique()
    ))
  } else if (dot.lable == "Description" && molecule.lable == "molecule") {
    df$node <- factor(df$node, levels = c(
      data_sankey$Description %>% unique(),
      data_sankey$target %>% unique(),
      data_sankey$molecule_id %>% unique(),
      data_sankey$herb %>% unique()
    ))
  } else if (dot.lable == "ID" && molecule.lable == "molecule") {
    df$node <- factor(df$node, levels = c(
      data_sankey$ID %>% unique(),
      data_sankey$target %>% unique(),
      data_sankey$molecule %>% unique(),
      data_sankey$herb %>% unique()
    ))
  }

  # color settings
  mycol <- cols4all::c4a("rainbow_wh_rd", length(unique(df$node)))
  mycol2 <- sample(mycol, length(mycol))

  # Sankey diagram
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
    ggsankey::theme_sankey(base_size = 20, base_family = "sans") +
    scale_fill_manual(values = mycol2) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(
      size = sankey.x.axis.text.size,
      hjust = sankey.x.axis.text.position,
      vjust = 10, colour = "black"
    ))
  p3 <- p2 + theme(plot.margin = unit(c(0, dot.position, 0, 0), units = "cm"))
  p4 <- cowplot::ggdraw() +
    cowplot::draw_plot(p3) +
    cowplot::draw_plot(p1,
      scale = dot.scale,
      x = dot.x, y = dot.y,
      width = dot.width,
      height = dot.height
    )
  return(p4)
}
