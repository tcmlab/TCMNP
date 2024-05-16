#' Bar graphs showed the results of tcm_prescription  analysis
#'
#' @param data tcm_prescription  results
#' @param top according to the order of p adjust value from small to large
#' the number of categories to show
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param text.size text size
#' @param text.width y-axis label length
#' @param text.bar The size of the text at the top of the bar graph
#' @param title  title
#' @param ... additional parameters
#'
#' @return barplot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom dplyr mutate
#' @importFrom dplyr slice
#' @importFrom stats reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#' data(disease_data, package = "TCMNP")
#' compound<- tcm_prescription(disease_data)[[2]]
#' tcm_pre_plot(compound)
#' }
tcm_pre_plot <- function(data,
                         top = 10,
                         color = "RdBu",
                         text.size = 10,
                         text.bar = 4,
                         text.width = 35,
                         title = NULL,
                         ...) {
  data<- data %>%
  as.data.frame()%>%
  arrange(desc(Count)) %>%
  dplyr::slice(1:top)
  data$Pvalue <- data$Pvalue %>% round(digits = 2)

  # ggplot2 plotting
  p <- ggplot(data) +
    geom_col(aes(x = Count, y = reorder(CompoundId, Count), fill = Pvalue),
             color = "black",
             width = 0.6
    ) +
    geom_text(aes(x = Count, y = reorder(CompoundId, Count), label = Count),
              hjust = -0.05,
              vjust = 0.5,
              size = text.bar
    ) +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = text.size),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = text.size,
        colour = "black"
      ),
      axis.text.y = element_text(
        angle = 0,
        size = text.size,
        face = "plain",
        colour = "black"
      )
    ) +
    scale_y_discrete(expand = c(0, -1)) +
    theme(legend.position = "right") +
    ggtitle(title) +
    ylab(NULL) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width))
  return(p)
}

