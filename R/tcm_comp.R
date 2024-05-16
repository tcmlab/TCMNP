#' @title Drug composition in prescription
#'
#' @param data data must consist of herb and weight columns
#' @param start.degree the angle at which the circle plot starts
#' @param color the color of circle plot.
#' see "RColorBrewer::display.brewer.all()"
#' @param radius.ratio from the outside to the inside,
#' the size of the proportion of each layer of the circle plot
#' @param text.size the font size of each layer circle
#'
#' @return figure
#' @export
#' @importFrom circlize circos.par
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.trackPlotRegion
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom circlize highlight.sector
#' @importFrom circlize colorRamp2
#' @importFrom circlize circos.clear
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom utils data
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' xfbdf.compostion <- data.frame(
#'   herb = c(
#'     "mahuang", "kuxingren", "shengshigao",
#'     "shengyiren", "maocangzhu", "guanghuoxiang",
#'     "qinghaocao", "mabiancao", "ganlugen", "tinglizi",
#'     "huajuhong", "shenggancao", "huzhang"
#'   ),
#'   weight = c(6, 15, 30, 30, 10, 15, 12, 30, 30, 15, 15, 10, 20)
#' )
#' tcm_comp(xfbdf.compostion)
#' }
tcm_comp <- function(data,
                     start.degree = 180,
                     color = "Paired",
                     radius.ratio = c(0.4, 0.2, 0.2),
                     text.size = c(0.9, 0.9, 0.9)) {
  data$label <- paste0(data$weight, " g")
  # Set to take two digits after the decimal point
  data$ratio <- paste0(base::format(data$weight / sum(data$weight) * 100,
                                    digits = 2, nsmall = 2
  ), "%") %>%
    as.character()
  # color settings
  plot.color <- grDevices::colorRampPalette(
    if (nrow(data) <= 8) {
      RColorBrewer::brewer.pal(nrow(data), color)
    } else {
      RColorBrewer::brewer.pal(8, color)
    }
  )(nrow(data))
  names(plot.color) <- data$herb
  data$number <- rep(1, nrow(data))
  sector_lim <- c(
    rep(0, nrow(data)),
    data$weight / sum(data$weight)
  ) %>% matrix(ncol = 2)
  rownames(sector_lim) <- data$herb
  # start drawing
  circlize::circos.clear()
  # The first circle: the name of each drug
  circlize::circos.par(
    "start.degree" = start.degree,
    points.overflow.warning = FALSE
  ) # set the starting angle
  circlize::circos.initialize(
    factors = data$herb,
    xlim = sector_lim
  )
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    factors = data$herb,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        mean(xlim),
        mean(ylim),
        sector.index,
        cex = text.size[1],
        facing = "clockwise",
        niceFacing = TRUE
      )
    }, track.height = radius.ratio[1], bg.col = plot.color
  )
  # Second circle: dosage of each drug
  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 2,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
    }, track.height = radius.ratio[2]
  )
  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 2,
      sector.index = (data$herb)[i],
      col = paste0(plot.color[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "clockwise",
      cex = text.size[2],
      adj = c(0.5, 0.5),
      labels = (data$label)[i],
      sector.index = (data$herb)[i]
    )
  }
  # Third circle: the percentage of each drug dose in the total dose
  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 3,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
    }, track.height = radius.ratio[3]
  )
  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 3,
      sector.index = (data$herb)[i],
      col = paste0(plot.color[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "clockwise",
      cex = text.size[3],
      adj = c(0.5, 0.5),
      labels = (data$ratio)[i],
      sector.index = (data$herb)[i]
    )
  }
  circlize::circos.clear()
}
