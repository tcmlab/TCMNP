#' Display the relationship between genes and pathways
#' enriched by KEGG or GO in the form of a chord diagram
#'
#' @param data  R clusterprofiler package for KEGG and GO results
#' The data format must be S4 object with KEGG and GO results or
#' a data frame with KEGG and GO results or
#' a data frame with two columns.
#'
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param start.degree the angle at which the circle plot starts
#' @param label.name "ID" or "Description"
#' @param text.width label length
#' @param text.size the font size of each layer circle
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param ... additional parameters
#'
#' @return chord diagram
#' @export
#'
#' @importFrom circlize circos.par
#' @importFrom circlize chordDiagram
#' @importFrom circlize circos.track
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom circlize colorRamp2
#' @importFrom circlize get.all.sector.index
#' @importFrom circlize circos.clear
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot ggdraw
#' @importFrom tidyr separate_rows
#'
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMNP")
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' library(DOSE)
#' library(tidyverse)
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
#' pathway_cirplot(KK)
#' }
pathway_cirplot <- function(data,
                            top = 5,
                            start.degree = -90,
                            label.name = "Description",
                            text.width = 15,
                            text.size = 0.6,
                            color = "RdBu", ...) {
  # data processing
  if (isS4(data)) {
    data <- data@result %>% tidyr::drop_na()
  } else if (is.data.frame(data)) {
    data <- data %>% tidyr::drop_na()
  } else {
    print("The data format must be S4 object or data frame.")
  }
  if (all(c("ID", "Description", "geneID") %in% colnames(data))) {
    df <- separate_rows(data[1:top, ], geneID, sep = "/") %>%
      dplyr::mutate(col = "grey") %>%
      dplyr::distinct()
    gene_col <- df$col %>% head(length(unique(df$geneID)))
    names(gene_col) <- unique(df$geneID)

    if (label.name == "ID") {
      df <- df %>% dplyr::select(ID, geneID)
      cols_number <- length(unique(df$ID))
      cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
      pathway_col <- colorRampPalette(cols)(length(unique(df$ID)))
      names(pathway_col) <- unique(df$ID)
    } else if (label.name == "Description") {
      df <- df %>% dplyr::select(Description, geneID)
      df$Description <- lapply(
        strwrap(df$Description,
          width = text.width,
          simplify = FALSE
        ),
        paste,
        collapse = "\n"
      )
      cols_number <- length(unique(df$Description))
      cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
      pathway_col <- colorRampPalette(cols)(length(unique(df$Description)))
      names(pathway_col) <- unique(df$Description)
    } else {
      print("The label.name is ID' or 'Description'. ")
    }
  } else {
    df <- separate_rows(data[1:top, ], 2, sep = "/")
    colnames(df)[1:2] <- c("Description", "geneID")
    df <- df %>%
      dplyr::mutate(col = "grey")%>%
      dplyr::distinct()
    gene_col <- df$col %>% head(length(unique(df$geneID)))
    names(gene_col) <- unique(df$geneID)

    df$Description <- lapply(
      strwrap(df$Description,
        width = text.width,
        simplify = FALSE
      ),
      paste,
      collapse = "\n"
    )
    cols_number <- length(unique(df$Description))
    cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
    pathway_col <- colorRampPalette(cols)(length(unique(df$Description)))
    names(pathway_col) <- unique(df$Description)
  }

  grid.col <- c(pathway_col, gene_col)
  set.seed(1234)
  circlize::circos.par(
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    start.degree = start.degree
  )
  chordDiagram(df,
    grid.col = grid.col,
    link.decreasing = TRUE,
    transparency = 0.1,
    big.gap = 1,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = .1)
  )

  for (si in get.all.sector.index()) {
    xlim <- get.cell.meta.data("xlim",
      sector.index = si,
      track.index = 1
    )
    ylim <- get.cell.meta.data("ylim",
      sector.index = si,
      track.index = 1
    )
    circos.text(mean(xlim),
      ylim[1],
      labels = si,
      sector.index = si,
      track.index = 1,
      facing = "clockwise",
      cex = text.size,
      adj = c(0, 0.5),
      niceFacing = TRUE
    )
  }
  circos.clear()
}
