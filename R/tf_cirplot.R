#' Display the relationship between transcription factors and target genes
#'
#' @param data data.frame
#' @param start.degree the angle at which the circle plot starts
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
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMNP")
#' tf_data <- tf_filter(xfbdf$target)
#' set.seed(1234)
#' data <- tf_data[, 1:2] %>%
#'   distinct() %>%
#'   sample_n(100)
#' tf_cirplot(data, color = "Spectral")
#' }
tf_cirplot <- function(data,
                       start.degree = 270,
                       text.size = 0.6,
                       color = "RdBu", ...) {
  # data processing
  if (is.data.frame(data)) {
    data <- data %>% dplyr::distinct()
    if (is.data.frame(data)) {
      if (length(unique(data[, 1])) > length(unique(data[, 2]))) {
        colnames(data) <- c("Target", "TF")
      } else {
        colnames(data) <- c("TF", "Target")
      }
    } else {
      print("The data must be a data frame with two columns.")
    }

    cols_number <- length(unique(data$TF))
    cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
    tf_col <- colorRampPalette(cols)(length(unique(data$TF)))
    names(tf_col) <- unique(data$TF)

    Target_col <- rep("lightgrey", length(unique(data$Target)))
    names(Target_col) <- unique(data$Target)
    grid.col <- c(tf_col, Target_col)

    # drawing
    circos.clear()
    circlize::circos.par(
      canvas.xlim = c(-1, 1),
      canvas.ylim = c(-1, 1),
      start.degree = start.degree
    )

    chordDiagram(data,
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
  } else {
    print("The data must be a data frame with two columns.")
  }
}
