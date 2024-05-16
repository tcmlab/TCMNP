#' the degree of a network node
#'
#' @param data data.frame
#' the column name must contain three columns:
#' herb, molecule, and target
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
#'
#' @param top According to the order of P adjust value from small to large
#' the number of categories to show
#' @param text.size  text size
#' @param text.width text width
#' @param plot.set "horizontal" or "vertical"
#' @param title title
#' @param ... additional parameters
#'
#' @return bar plot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 coord_flip
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr desc
#' @importFrom stringr str_wrap
#'
#' @examples
#' degree_plot(xfbdf)
degree_plot <- function(data,
                        color = c("#A40545", "#7FCBA4", "#4B65AF"),
                        top = 10,
                        text.size = 12,
                        text.width = 20,
                        plot.set = "vertical",
                        title = "the degree of networks",
                        ...) {
  # data processing
  deg <- data %>% dplyr::mutate(., value = 1)
  herb_degree <- aggregate(deg$value, by = list(deg$herb), length) %>%
    as.data.frame() %>%
    dplyr::mutate(., type = "herb")
  molecule_degree <- aggregate(deg$value, by = list(deg$molecule), length) %>%
    as.data.frame() %>%
    dplyr::mutate(., type = "molecule")
  target_degree <- aggregate(deg$value, by = list(deg$target), length) %>%
    as.data.frame() %>%
    dplyr::mutate(., type = "target")
  deg2 <- rbind(herb_degree, molecule_degree, target_degree)
  colnames(deg2) <- c("id", "degree", "type")
  deg3 <- deg2 %>%
    dplyr::select("id", "type", "degree") %>%
    as.data.frame()
  data_sorted <- deg3 %>%
    dplyr::group_by(type) %>%
    top_n(top, wt = degree) %>%
    arrange(., type, desc(degree))
  data_sorted$id <- factor(data_sorted$id, levels = data_sorted$id)

  # Vertical version
  if (plot.set == "vertical") {
    p <- ggplot(data = data_sorted, aes(x = id, y = degree, fill = type)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width)) +
      scale_fill_manual(values = color) +
      theme_bw() +
      theme(axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = text.size + 1
      )) +
      ylab("degree") +
      xlab("") +
      ggtitle(title) +
      theme(
        axis.text = element_text(size = text.size, colour = "black"),
        plot.title = element_text(size = text.size + 2),
        legend.title = element_text(size = text.size),
        legend.text = element_text(size = text.size),
      ) +
      coord_flip() +
      scale_y_continuous(expand = c(0, 0))
    return(p)
  } else if (plot.set == "horizontal") {
    # hertical version
    p <- ggplot(data = data_sorted, aes(x = id, y = degree, fill = type)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = color) +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = text.width)) +
      ylab("degree") +
      theme(legend.title = element_blank(), legend.position = "top") +
      theme(
        legend.key.size = unit(text.size, "pt"),
        legend.text = element_text(size = text.size + 2)
      ) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = text.size + 1),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = text.size, colour = "black")
      ) +
      scale_y_continuous(expand = c(0, 0))
    return(p)
  } else {
    print("The plot.set must be 'horizontal' or 'vertical'.")
  }
}
