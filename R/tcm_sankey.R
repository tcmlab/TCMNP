#' Sankey diagram
#'
#' @param data_sankey data.frame
#' @param text.size text size
#' @param text.position text position 0, 0.5, 1
#' @param x.axis.text.size Location of x-axis labels in sankey plot
#' @param ... additional parameters
#'
#' @return figure
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 position_nudge
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom dplyr mutate
#' @importFrom ggsankey make_long
#' @importFrom ggsankey geom_sankey
#' @importFrom ggsankey geom_sankey_text
#' @importFrom ggsankey theme_sankey
#' @importFrom cols4all c4a
#'
#' @examples
#' data("xfbdf", package = "TCMNP")
#' data <- xfbdf %>% dplyr::sample_n(30, replace = FALSE)
#' tcm_sankey(data, text.position = 1)

tcm_sankey <- function(data_sankey,
                       text.size = 3,
                       text.position = 0,
                       x.axis.text.size = 12,
                       ...){
  # color settings
  df <- data_sankey %>%
        as.data.frame()%>%
        make_long(colnames(data_sankey))
  mycol <- cols4all::c4a("rainbow_wh_rd", length(unique(df$node)))
  mycol2 <- sample(mycol, length(mycol), replace = FALSE)
  # Sankey diagram
  p <- do.call(rbind, apply(data_sankey, 1, function(x) {
    data.frame(
      x = names(x), node = x,
      next_x = dplyr::lead(names(x)),
      next_node = dplyr::lead(x), row.names = NULL
    )
  })) %>%
    dplyr::mutate(
      x = factor(x, names(data_sankey)),
      next_x = factor(next_x, names(data_sankey))
    ) %>%
    ggplot(aes(
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
      size = text.size, color = "black",
      hjust = text.position,
      position = position_nudge(x = -0.05)
    ) +
    ggsankey::theme_sankey(base_size = 20) +
    scale_fill_manual(values = mycol2) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(
      size = x.axis.text.size,
      hjust = 0.5, vjust = 10,
      colour = "black"
    ))
  return(p)
}
