#' Herb-ingredient-target-function enrichment analysis: alluvial plot
#'
#' @param data data.frame
#' @param alluvial.text.size alluvial diagram text size
#' @param alluvial.x.axis.text.position
#' Location of x-axis labels in alluvial plot
#' @param alluvial.x.axis.text.size alluvial diagram x-axis text size
#' @param molecule.lable "molecule" or ""molecule_id"
#' @param dot.color dotplot color: see "RColorBrewer::display.brewer.all()"
#' @param dot.position dotplot position
#' @param dot.lable "ID" or "Description"
#' @param dot.text.size  dotplot text size
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
#' @import dplyr
#' @import ggplot2
#' @importFrom ggsankey make_long
#' @importFrom ggsankey geom_alluvial
#' @importFrom ggsankey geom_alluvial_text
#' @importFrom ggsankey theme_sankey
#' @importFrom cols4all c4a
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @importFrom stats reorder
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
#' library(tidyverse)
#' data(xfbdf)
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
#'   mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#' path <- separate_rows(KK@result, geneID, sep = "/")
#' data_sankey <- left_join(xfbdf, path,
#'   by = c("target" = "geneID"),
#'   relationship = "many-to-many"
#' ) %>%
#'   distinct() %>%
#'   drop_na() %>%
#'   sample_n(30, replace = FALSE) %>%
#'   as.data.frame()
#' alluvial_dot(data_sankey)
#' }
alluvial_dot <- function(data,
                         alluvial.text.size = 3,
                         alluvial.x.axis.text.position = 1,
                         alluvial.x.axis.text.size = 12,
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
  options(digits = 2) # Set to two digits after the decimal point
  data_kegg <- data %>%
    dplyr::select("ID", "Description", "p.adjust", "Count", "richFactor") %>%
    distinct() %>%
    arrange(p.adjust, descreasing = TRUE)
  data_kegg$richFactor <- data_kegg$richFactor %>% round(2)
  # bubble chart
  if (dot.lable == "Description") {
    p1 <- ggplot(
      data = data_kegg,
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
          colour = "black",
          vjust = 1
        ),
        axis.text.y = element_text(
          size = dot.text.size,
          colour = "black",
          hjust = 1
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
      data = data_kegg,
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
          colour = "black",
          vjust = 1
        ),
        axis.text.y = element_text(
          size = dot.text.size,
          colour = "black",
          hjust = 1
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
    arrange(richFactor, descreasing = FALSE) %>%
    dplyr::select(herb, molecule_id, target, dot.lable) %>%
    as.data.frame() %>%
    distinct()

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
    ggsankey::geom_alluvial(
      flow.alpha = 0.5,
      node.color = NA,
      show.legend = FALSE
    ) +
    ggsankey::geom_alluvial_text(
      size = alluvial.text.size,
      color = "black",
      hjust = 1,
      position = position_nudge(x = -0.05)
    ) +
    theme_sankey(base_size = 20) +
    scale_fill_manual(values = mycol2) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(
      size = alluvial.x.axis.text.size,
      hjust = alluvial.x.axis.text.position,
      vjust = 10, colour = "black"
    ))
  p3 <- p2 + theme(plot.margin = unit(c(0, dot.position, 0, 0), units = "cm"))
  cowplot::ggdraw() +
    cowplot::draw_plot(p3) +
    cowplot::draw_plot(p1,
      scale = dot.scale,
      x = dot.x, y = dot.y,
      width = dot.width,
      height = dot.height
    )
}
