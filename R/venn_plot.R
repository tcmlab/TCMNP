#' Intersection of datasets
#'
#' @param data data.frame
#' column names must be two columns of data for gene and source
#' @param set.color color
#' Recommended color:
#' venn.color1 <- c("#b1d1e7", "#f0a4af", "#e5dfae", "#a8e1da", "#f8d068")
#' venn.color2 <- c("#d3e2b7", "#fbca50", "#c9d5d4", "#baa28a", "#4caecc")
#' venn.color3 <- c("#C1BCBF", "#B1D0A9", "#FFCCC3", "#FFD5AB", "#BCCBE5")
#' venn.color4 <- c("#e4c9b2", "#92C2DD", "#f49d98", "#fcd68f", "#629076")
#' venn.color5 <- c("#FF8748", "#5BAA56", "#B8BB5B", "#4186B7", "#8679BE")
#' venn.color6 <- c("#e6c09e", "#0d5888", "#cb8b3e", "#9cd6d6", "#dbcb09")
#' venn.color7 <- c("#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#75cbdc")
#' venn.color8 <- c("#604192", "#139177", "#ed9e08", "#f56f08", "#4caecc")
#' venn.color9 <- c("#d37a20", "#dbcb09", "#3a9cbc", "#dd7208", "#a30019")
#' @param bar.color top column color
#' @param type  1, 2
#' 1:take the intersection of 2-4 data sets
#' 2:take the intersection of >= 4 data sets
#' @param name.size name size
#' @param text.size text size
#' @param vratio The ratio of the histogram to the size of the matrix
#' @param ... additional parameters
#'
#' @return figure
#'
#' @importFrom dplyr filter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggvenn ggvenn
#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr pivot_wider
#' @importFrom UpSetR upset
#' @importFrom UpSetR intersects
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#'  \dontrun{
#' library(UpSetR)
#' data(venn_data, package = "TCMNP")
#' venn_plot(venn_data, type =1)
#' venn_plot(venn_data, type =2)
#' }
#'
venn_plot <- function(data,
                      set.color = c(
                        "#FF8748", "#5BAA56", "#4186B7",
                        "#8679BE", "#B8BB5B"
                      ),
                      name.size = 5,
                      text.size = 4,
                      bar.color = "Set1",
                      type = 1,
                      vratio = c(0.65, 0.35), ...) {
  # data processing
  data <- as.data.frame(data)
  if (is.data.frame(data)) {
    data <- data %>% dplyr::distinct()
    if (is.data.frame(data)) {
      if (length(unique(data[, 1])) > length(unique(data[, 2]))) {
        colnames(data) <- c("gene", "source")
      } else {
        colnames(data) <- c("source", "gene")
      }
    } else {
      print("The data must be a data frame with two columns.")
    }

    database <- unique(data$source)
    df <- list()
    df2 <- list()
    for (i in seq_along(database)) {
      df[[i]] <- data %>% dplyr::filter(source == database[i])
      df2[[i]] <- df[[i]]$gene
    }
    names(df2) <- database
    # draw a between-sample Venn diagram
    if (type == 1) {
      color_set <- set.color[seq_along(database)]
      df2 %>% ggvenn::ggvenn(
        show_percentage = TRUE,
        show_elements = FALSE,
        digits = 2,
        stroke_color = "white",
        stroke_alpha = 1,
        set_name_size = name.size,
        text_size = text.size,
        show_outside = c("auto"),
        fill_color = color_set,
        set_name_color = color_set
      )+
        scale_x_continuous(expand = expansion(mult = 0.2))
    } else if (type == 2) {
      # draw a third Venn diagram
      data2 <- data %>%
        dplyr::group_by(gene, source) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = source, values_from = n) %>%
        as.data.frame()
      data2[is.na(data2)] <- 0
      # create the color first, set as many colors as there are columns:
      data2 <- as.data.frame(data2)
      tmp2 <- unique(data2[, -1])
      tmp2 <- tmp2[rowSums(tmp2) != 0, ]
      # top column color
      colors2 <- colorRampPalette(brewer.pal(8, bar.color))(nrow(tmp2))
      query_list2 <- list()
      for (i in seq_len(nrow(tmp2))) {
        query_list2[[i]] <- list(
          query = UpSetR::intersects,
          params = list(colnames(tmp2)[which(tmp2[i, ] == 1)]),
          color = colors2[i],
          active = TRUE
        )
      }
      # left column color
      if (length(database) <= 5) {
        venn.color <- set.color[seq_along(database)]
      } else if (length(database) > 5) {
        venn.color <- grDevices::colorRampPalette(
          brewer.pal(length(set.color), set.color)(length(database))
        )
      }
      UpSetR::upset(data2,
                    nsets = ncol(tmp2),
                    text.scale = text.size/4,
                    # number of datasets
                    mb.ratio = vratio,
                    # the ratio of the histogram to the size of the matrix
                    #sets.bar.color = venn.color,
                    # modify the color of the left bar
                    order.by = "freq",
                    # column sorting
                    decreasing = TRUE,
                    queries = query_list2,
                    # modify the color of the bars and the color of the matrix scatter
                    shade.color = NA
                    # shadow color of rectangular scatter points
      )
    }
  }
}


