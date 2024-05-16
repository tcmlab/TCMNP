#' Molecular docking scoring heatmap
#'
#' @param data Molecular docking results, data frame
#' @param shape "square" or "circle" 
#' @param color color: see "RColorBrewer::display.brewer.all()", brewer.pal(9,"Set1")
#' @param text.size text size
#' @param legend.mid legend median
#' @param legend.width legend width
#' @param legend.height legend height
#' @param ... additional parameters
#'
#' @return figure
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 unit
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom stats hclust
#' @importFrom stats dist
#' @importFrom aplot insert_left
#' @importFrom aplot insert_top
#' @importFrom ggtree ggtree
#' @importFrom ggtree layout_dendrogram
#'
#' @examples
#' data <- matrix(rnorm(81), 9, 9)
#' data[1:9, seq(1, 9, 2)] <- data[1:9, seq(1, 9, 2)] - 4
#' colnames(data) <- paste("molecule", 1:9, sep = "")
#' rownames(data) <- paste("target", 1:9, sep = "")
#' data <- round(data, digits = 2)
#' dock_plot(data,text.size = 4)
#' dock_plot(data, shape = "circle", text.size = 4,legend.height = 3)
dock_plot <- function(data,
                      shape = 'square',
                      color = c("#4393C3", "#F7F7F7", "#D6604D"),
                      text.size = 4,
                      legend.mid = -2,
                      legend.width = 0.5,
                      legend.height = 5, ...) {
  # remove rownames
  remove.rownames <- function(x) {
    stopifnot(is.data.frame(x))
    rownames(x) <- NULL
    return(x)
  }
  # Heat maps were generated in R using the ggplot2 package
  if (shape=='square') {
    p1 <- data %>%
      as.data.frame() %>%
      dplyr::mutate(rowname = rownames(data), .before = 1) %>%
      remove.rownames() %>%
      tidyr::pivot_longer(
        cols = 2:ncol(.),
        names_to = "molecule",
        values_to = "score"
      ) %>%
      ggplot(aes(x = molecule, y = rowname)) +
      geom_tile(aes(fill = score)) +
      scale_fill_gradient2(
        midpoint = legend.mid,
        low = color[1],
        mid = color[2],
        high = color[3]
      ) +
      scale_y_discrete(position = "right") +
      geom_text(aes(label = score), size = text.size) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          vjust = 0.5,
          hjust = 0.5,
          colour = "black",
          size = text.size*3
        ),
        axis.text.y = element_text(colour = "black",
                                   size = text.size*3),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      ) +
      guides(fill = guide_colorbar(
        direction = "vertical", reverse = FALSE,
        barwidth = unit(legend.width, "cm"),
        barheight = unit(legend.height, "cm")
      ))

    # cluster tree
    yclust <- stats::hclust(dist(data))
    xclust <- stats::hclust(dist(t(data)))
    p2 <- ggtree::ggtree(yclust)
    p3 <- ggtree::ggtree(xclust) + ggtree::layout_dendrogram()

    ## merge with aplot
    p4 <- p1 %>%
      aplot::insert_left(p2, width = 0.2) %>%
      aplot::insert_top(p3, height = 0.1)
    return(p4)
  } else if (shape == "circle") {
    # graphs were generated in R using the ggplot2 package
    p1 <- data %>%
      as.data.frame() %>%
      dplyr::mutate(rowname = rownames(data), .before = 1) %>%
      remove.rownames() %>%
      tidyr::pivot_longer(
        cols = 2:ncol(.),
        names_to = "molecule",
        values_to = "score"
      ) %>%
      ggplot(aes(x = molecule, y = rowname)) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()
      ) +
      geom_point(aes(size = abs(score), color = score),
        shape = 19, show.legend = TRUE
      ) +
      scale_colour_gradient2(
        midpoint = legend.mid,
        low = color[1],
        mid = color[2],
        high = color[3]
      ) +
      labs(x = NULL, y = NULL) +
      geom_text(aes(label = score),
        size = text.size,
        colour = "black"
      ) +
      scale_y_discrete(position = "right") +
      theme(
        axis.text.x = element_text(
          color = "black",
          angle = 45,
          vjust = 0.5,
          hjust = 0.5,
          size = text.size*3
        ),
        axis.text.y = element_text(color = "black",size = text.size*3),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(
          fill = NA,
          color = "grey80",
          size = 1,
          linetype = "solid"
        )
      ) +
      guides(color = guide_colorbar(
        reverse = FALSE, barwidth = unit(.5, "cm"),
        barheight = unit(5, "cm")
      ))

    # cluster tree
    yclust <- stats::hclust(dist(data))
    xclust <- stats::hclust(dist(t(data)))
    p2 <- ggtree::ggtree(yclust)
    p3 <- ggtree::ggtree(xclust) +
      ggtree::layout_dendrogram()

    # merge with aplot
    p4 <- p1 %>%
      aplot::insert_left(p2, width = 0.2) %>%
      aplot::insert_top(p3, height = 0.1)
    return(p4)
  }
}

