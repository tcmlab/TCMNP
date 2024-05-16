#' Chord diagram display of GO results
#'
#' @param go.diff R clusterprofiler package for enrichGO results
#' @param top  According to the order of P adjust value from small to large
#' the number of categories to show
#' @param id.color The color of the first layer ID of the path
#' Recommended color:
#' colour1 <- c("#c0583d", "#daa152", "#b3c6dc")
#' colour2 <- c("#8fa1c7", "#ef916a", "#7dbfa5")
#' colour3 <- c("#7dbfa5", "#487eb3", "#d2342b")
#' colour4 <- c("#6176b5", "#DFC86B", "#87BD5B")
#' colour5 <- c("#852f88", "#eb990c", "#0f8096")
#' colour6 <- c("#3c4b8d", "#da2e20", "#3c884c")
#' colour7 <- c("#7aa354", "#d67d44", "#6594c4")
#' colour8 <- c("#501d8a", "#e55709", "#1c8041")
#' colour9 <- c("#A40545", "#7FCBA4", "#4B65AF")
#' @param id.cex The text size of the first layer ID of the path
#' @param gene.color The color of the number of all genes in this pathway
#' @param gene.cex The font size of the number of all genes in this pathway
#' @param ratio.color
#' （1）The color of the number of genes enriched this time on this pathway
#' （2）Other genes involved in this enrichment analysis
#' @param ratio.cex
#' The text size of the numbers involved in this enrichment analysis
#' @param p.color p adjust color
#' @param p.cex The text size of padjust
#' @param richFactor.cex The text size of richFactor
#' @param legend.position legend position
#' @param text.width legend text width
#' @param ... additional parameters
#'
#' @return chord diagram
#' @export
#'
#' @importFrom circlize circos.par
#' @importFrom circlize circos.genomicInitialize
#' @importFrom circlize circos.genomicRect
#' @importFrom circlize circos.track
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom circlize colorRamp2
#' @importFrom circlize circos.clear
#' @importFrom circlize circos.genomicTrackPlotRegion
#' @importFrom circlize circos.lines
#' @importFrom circlize circos.genomicTrack
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr arrange
#' @importFrom grDevices colorRampPalette
#' @importFrom stringr str_split
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap packLegend
#' @importFrom stringr str_wrap
#' @importFrom stringr str_split
#' @importFrom tidyr drop_na
#' @importFrom grid gpar
#' @importFrom grid pushViewport
#' @importFrom grid grid.draw
#' @importFrom grid viewport
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
#'   readable = T
#' )
#' go_cir(go.diff)
#' }
go_cir <- function(
    go.diff,
    top = 5,
    id.color = c("#852f88", "#eb990c", "#0f8096"),
    id.cex = 0.8,
    gene.color = "#CFB0D4",
    gene.cex = 0.8,
    ratio.color = c("#6291bd", "#B1D1E7"),
    ratio.cex = 0.8,
    p.color = c("#eda9aa", "#c25254"),
    p.cex = 0.8,
    richFactor.cex = 0.7,
    legend.position = c(0.8, 0.45),
    text.width = 35,
    ...) {
  # data processing
  if (isS4(go.diff)) {
    go.diff <- go.diff@result %>% tidyr::drop_na()
  } else if (is.data.frame(go.diff)) {
    go.diff <- go.diff %>% tidyr::drop_na()
  } else {
    print("The data format must be S4 object or data frame.")
  }
  go_enrich <- go.diff %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::slice(1:top) %>%
    dplyr::mutate(RichFactor = Count / as.numeric(sub(
      "/\\d+", "",
      BgRatio
    ))) %>%
    dplyr::distinct(p.adjust, RichFactor, .keep_all = TRUE) %>%
    as.data.frame()
  go_enrich$RichFactor <- go_enrich$RichFactor %>% round(digits = 2)
  go_enrich$left <- 0
  go_enrich$this_pathway_gene <- go_enrich$BgRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][1]
    }) %>%
    as.numeric()
  go_enrich$all_pathway_gene <- go_enrich$BgRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][2]
    }) %>%
    as.numeric()
  go_enrich$DEGnum <- go_enrich$GeneRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][2]
    }) %>%
    as.numeric()
  go_enrich$right <- max(go_enrich$this_pathway_gene)
  go_enrich$ONTOLOGY <- factor(go_enrich$ONTOLOGY,
    levels = sort(unique(go_enrich$ONTOLOGY))
  )
  go_enrich <- go_enrich %>% dplyr::arrange(ONTOLOGY, desc(RichFactor))
  rownames(go_enrich) <- go_enrich$ID

  # start drawing
  plotdata <- go_enrich %>% drop_na()
  circos.clear()
  circos.par(
    canvas.xlim = c(-0.2, 1.2),
    canvas.ylim = c(-1.0, 1.2),
    clock.wise = TRUE,
    start.degree = 90,
    gap.degree = 0.8,
    xaxis.clock.wise = TRUE
  )
  ### The first circle: classification information
  color.selected <- id.color[1:3]
  color.order <- rep(color.selected, each = top)
  # Initialize circular plot with any genomic data
  circos.genomicInitialize(
    plotdata[, c(
      "ID", "left",
      "right", "ONTOLOGY"
    )],
    plotType = NULL
  )

  circos.track(
    ylim = c(0, 1), track.height = 0.08,
    bg.border = NA, bg.col = color.order,
    panel.fun = function(x, y) {
      ylim <- get.cell.meta.data("ycenter")
      xlim <- get.cell.meta.data("xcenter")
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(xlim, ylim, sector.name,
        cex = id.cex, niceFacing = TRUE
      )
    }
  )
  ### The second circle: How many genes are there in total in this pathway
  plotdata2 <- plotdata[, c("ID", "left", "this_pathway_gene")]
  circos.genomicTrackPlotRegion(
    plotdata2,
    track.height = 0.08, bg.border = NA,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = gene.color, border = NA, ...)
      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata2[get.cell.meta.data("sector.index"), 3] / 2
      sector.name <- plotdata2[get.cell.meta.data("sector.index"), 3]
      circos.text(xlim, ylim, sector.name,
        cex = gene.cex, niceFacing = TRUE
      )
    }
  )
  ### The third circle: generation related
  plotdata3 <- plotdata[, c("ID", "left", "Count", "DEGnum", "right")]
  plotdata3$ratio <- plotdata3$Count / plotdata3$DEGnum
  plotdata3$len <- plotdata3$ratio * plotdata3$right
  plotdata3$len2 <- plotdata3$right - plotdata3$len

  tmpdf1 <- plotdata3[, c("ID", "left", "len")]
  colnames(tmpdf1) <- c("ID", "start", "end")
  tmpdf1$type <- 1
  tmpdf2 <- plotdata3[, c("ID", "len", "right")]
  colnames(tmpdf2) <- c("ID", "start", "end")
  tmpdf2$type <- 2
  tmpdata <- tmpdf1 %>% rbind(tmpdf2)
  color_assign <- colorRamp2(breaks = c(1, 2), colors = ratio.color)

  circos.genomicTrackPlotRegion(
    tmpdata,
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        col = color_assign(value),
        border = NA, ...
      )

      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata3[get.cell.meta.data("sector.index"), "len"] / 2
      sector.name <- plotdata3[get.cell.meta.data("sector.index"), "Count"]
      circos.text(xlim, ylim, sector.name, cex = ratio.cex, niceFacing = TRUE)

      xlim <- plotdata3[get.cell.meta.data("sector.index"), "len"] +
        plotdata3[get.cell.meta.data("sector.index"), "len2"] / 2
      sector.name <- plotdata3[get.cell.meta.data("sector.index"), "DEGnum"] -
        plotdata3[get.cell.meta.data("sector.index"), "Count"]
      circos.text(xlim, ylim, sector.name,
        cex = ratio.cex, niceFacing = TRUE
      )
    }
  )
  ### Fourth circle: p-value
  plotdata4 <- plotdata[, c("ID", "left", "right", "p.adjust")]
  total.len <- unique(plotdata4$right)
  plotdata4$p.adjust_neg <- -log10(plotdata4$p.adjust)
  plotdata4$relative_value <-
    plotdata4$p.adjust_neg / max(plotdata4$p.adjust_neg) * total.len

  #round up integer-valued function
  my_ceiling <- function(x) {
    if (x > 0 & x - as.integer(x) != 0) {
      return(as.integer(x) + 1)
    } else {
      return(as.integer(x))
    }
  }
  p_max <- max(plotdata4$p.adjust_neg) %>% my_ceiling()
  colorsChoice <- colorRampPalette(p.color)
  color_assign <- colorRamp2(
    breaks = 0:p_max,
    colors = colorsChoice(p_max + 1)
  )

  plotdata4 <- plotdata4[, c(
    "ID", "left",
    "relative_value",
    "p.adjust_neg"
  )]

  circos.genomicTrackPlotRegion(
    plotdata4,
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        col = color_assign(value),
        border = NA, ...
      )

      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata4[
        get.cell.meta.data("sector.index"),
        "relative_value"
      ] / 2
      sector.name <- plotdata4[
        get.cell.meta.data("sector.index"),
        "p.adjust_neg"
      ] %>% round(digits = 2)
      circos.text(xlim, ylim, sector.name,
        cex = p.cex, niceFacing = TRUE
      )
    }
  )
  ### Fifth Circle: Enrichment Factor
  plotdata5 <- plotdata[, c("ID", "left", "right", "RichFactor", "ONTOLOGY")]
  color_assign <- color.selected
  names(color_assign) <- as.character(unique(plotdata5$ONTOLOGY))

  circos.genomicTrack(
    plotdata5,
    ylim = c(0, max(plotdata5$RichFactor)),
    track.height = 0.45, bg.col = "gray95", bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.genomicRect(region, value,
        col = color_assign[plotdata5[sector.name, 5]],
        border = NA,
        ytop.column = 1,
        ybottom = 0, ...
      )
      ylim <- plotdata5[
        get.cell.meta.data("sector.index"),
        "RichFactor"
      ] / 2
      xlim <- get.cell.meta.data("xcenter")
      sector.name <- plotdata5[
        get.cell.meta.data("sector.index"),
        "RichFactor"
      ]
      circos.text(xlim, ylim, sector.name, cex = richFactor.cex, niceFacing = TRUE)
    }
  )

  circos.clear()
  # draw legend
  category_legend1 <- ComplexHeatmap::Legend(
    labels = as.character(unique(
      paste0(
        plotdata$ONTOLOGY, ":",
        rep(
          c(
            "Biological process",
            "Cellular component",
            "Molecular function"
          ),
          each = top
        )
      )
    )),
    background = color.selected,
    type = "points", pch = NA,
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  label.text <- as.character(unique(paste0(
    plotdata$ID, ":",
    plotdata$Description
  )))
  category_legend2 <- ComplexHeatmap::Legend(
    labels = stringr::str_wrap(label.text, width = text.width),
    background = color.order,
    type = "points",
    pch = NA,
    labels_gp = grid::gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )


  thispathway_gene_legend <- Legend(
    labels = c("gene number of all this pathway"),
    background = gene.color,
    type = "points", pch = NA,
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  generatio_legend <- Legend(
    labels = c("gene number in this pathway", "other genes"),
    background = ratio.color,
    type = "points", pch = NA,
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  pvalue_legend <- Legend(
    col_fun = colorRamp2(
      breaks = 0:p_max,
      colors = colorRampPalette(p.color)(p_max + 1)
    ),
    legend_height = unit(3, "cm"),
    labels_gp = gpar(fontsize = 8),
    title = "-log10(p.adjust)",
    title_gp = gpar(fontsize = 9),
    title_position = "lefttop",
    direction = "horizontal"
  )

  pack_Legend <- packLegend(
    category_legend1,
    category_legend2,
    thispathway_gene_legend,
    generatio_legend,
    pvalue_legend,
    row_gap = unit(0.2, "cm")
  )

  grid::pushViewport(grid::viewport(
    x = legend.position[1],
    y = legend.position[2]
  ))
  grid::grid.draw(pack_Legend)
}
