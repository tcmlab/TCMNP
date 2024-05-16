#' Bar graphs showed the results of GO results
#'
#' @param go.diff R clusterprofiler package for enrichGO results
#' @param top  According to the order of P adjust value from small to large
#' the number of categories to show
#' @param plot.set "horizontal" or "vertical"
#' @param facet facet
#' @param color color, Recommended color:
#' colour1 <- c("#c0583d", "#daa152", "#b3c6dc")
#' colour2 <- c("#8fa1c7", "#ef916a", "#7dbfa5")
#' colour3 <- c("#7dbfa5", "#487eb3", "#d2342b")
#' colour4 <- c("#6176b5", "#DFC86B", "#87BD5B")
#' colour5 <- c("#852f88", "#eb990c", "#0f8096")
#' colour6 <- c("#3c4b8d", "#da2e20", "#3c884c")
#' colour7 <- c("#7aa354", "#d67d44", "#6594c4")
#' colour8 <- c("#501d8a", "#e55709", "#1c8041")
#' colour9 <- c("#A40545", "#7FCBA4", "#4B65AF")
#' @param bar.width  bar width
#' @param text.size  text size
#' @param text.width text width
#' @param title title
#' @param ... additional parameters
#'
#' @return bar plot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 facet_grid
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom tidyr drop_na
#' @importFrom stats reorder
#' @importFrom stringr str_wrap
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
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
#'   readable = T
#' )
#' go_barplot(go.diff)
#' }
go_barplot <- function(go.diff,
                       top = 5,
                       plot.set = "vertical",
                       facet = TRUE,
                       color = c("#A40545", "#7FCBA4", "#4B65AF"),
                       bar.width = 0.8,
                       text.size = 12,
                       text.width = 30,
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
    dplyr::slice(1:top)
  go_enrich <- go_enrich[order(go_enrich$ONTOLOGY,
    go_enrich$p.adjust,
    decreasing = TRUE
  ), ]
  go_enrich$Description <- factor(go_enrich$Description,
    levels = go_enrich$Description
  )
  # ggplot2 plotting
  if (plot.set == "vertical" && facet == FALSE) {
    p1 <- ggplot(go_enrich, aes(Description, -log10(p.adjust))) +
      geom_bar(aes(fill = ONTOLOGY),
        width = bar.width,
        stat = "identity"
      ) +
      geom_text(
        aes(
          label = round(-log10(p.adjust), digits = 0),
          y = round(-log10(p.adjust), digits = 0) + 1.5
        ),
        size = text.size / 3
      ) +
      coord_flip() +
      labs(x = "", y = "-log10 p.adjust", title = title) +
      labs(fill = "type") +
      scale_fill_manual(values = color) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          size = text.size,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text.size,
          colour = "black"
        )
      ) +
      theme(axis.title.x = element_text(size = text.size)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    p1
  } else if (plot.set == "vertical" && facet == TRUE) {
    p2 <- ggplot(go_enrich, aes(Description, -log10(p.adjust))) +
      geom_bar(aes(fill = ONTOLOGY),
        width = bar.width,
        stat = "identity"
      ) +
      geom_text(
        aes(
          label = round(-log10(p.adjust), digits = 0),
          y = round(-log10(p.adjust), digits = 0) + 1.5
        ),
        size = text.size / 3
      ) +
      coord_flip() +
      labs(x = "", y = "-log10 p.adjust", title = title) +
      labs(fill = "type") +
      scale_fill_manual(values = color) +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = text.size, colour = "black"),
        axis.text.y = element_text(size = text.size, colour = "black")
      ) +
      theme(axis.title.x = element_text(size = text.size)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme(strip.text = element_text(face = "bold", size = text.size)) +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      theme(legend.position = "none")
    p2
  } else if (plot.set == "horizontal" && facet == FALSE) {
    p3 <- ggplot(go_enrich, aes(Description, -log10(p.adjust))) +
      geom_col(aes(fill = ONTOLOGY), width = bar.width) +
      geom_text(
        aes(
          label = round(-log10(p.adjust), digits = 0),
          y = round(-log10(p.adjust), digits = 0) + 1
        ),
        size = text.size / 3
      ) +
      scale_fill_manual(values = color) +
      theme(
        axis.text.x = element_text(
          size = text.size,
          colour = "black",
          angle = 60,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = element_text(
          size = text.size,
          colour = "black"
        )
      ) +
      theme(axis.title.y = element_text(size = text.size)) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(
          color = "black",
          fill = "transparent"
        )
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "", y = "-Log10 p.adjust") +
      labs(fill = "type") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width))
    p3
  } else if (plot.set == "horizontal" && facet == TRUE) {
    p4 <- ggplot(go_enrich, aes(Description, -log10(p.adjust))) +
      geom_col(aes(fill = ONTOLOGY), width = bar.width) +
      geom_text(
        aes(
          label = round(-log10(p.adjust), digits = 0),
          y = round(-log10(p.adjust), digits = 0) + 1
        ),
        size = text.size / 3
      ) +
      scale_fill_manual(values = color) +
      theme(
        axis.text.x = element_text(
          size = text.size,
          colour = "black",
          angle = 60,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = element_text(size = text.size, colour = "black")
      ) +
      theme(axis.title.y = element_text(size = text.size)) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(
          color = "black",
          fill = "transparent"
        )
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "", y = "-Log10 p.adjust") +
      labs(fill = "type") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width)) +
      facet_grid(. ~ ONTOLOGY, scales = "free_x", space = "free_x") +
      theme(strip.text = element_text(
        face = "bold",
        size = text.size
      )) +
      theme(legend.position = "none")
    p4
  } else {
    print("The plot.set must be 'horizontal' or 'vertical'.")
  }
}
