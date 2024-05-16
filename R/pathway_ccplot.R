#' circle chart showed the results of GO and KEGG analysis
#'
#' @param data R clusterprofiler package for KEGG and GO results
#' The data format must be S4 object with KEGG and GO results or
#' a data frame with KEGG and GO results or
#' a data frame with two columns.
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param label.name "ID" or "Description"
#' @param root root name
#' @param color.node node color
#' @param color.alpha color alpha
#' @param text.size text size
#' (1) first circle text size
#' (2) Second circle text size
#' @param ... additional parameters
#'
#' @return figure
#' @export
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate_at
#' @importFrom dplyr summarise
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_fan
#' @importFrom ggraph geom_edge_diagonal
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph scale_edge_width
#' @importFrom ggraph node_angle
#' @importFrom ggraph theme_graph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr separate_rows
#' @importFrom tidygraph tbl_graph
#' @importFrom tidyr unite
#' @importFrom utils tail
#' @examples
#' \dontrun{
#' data("KK", package = "TCMNP")
#' library(dplyr)
#' pathway_ccplot(KK, root = "KEGG")
#' }
#'
pathway_ccplot <- function(data,
                           top = 5,
                           label.name = "Description",
                           root = NULL,
                           color.node = "Paired",
                           color.alpha = 0.5,
                           text.size = c(3, 4),
                           ...) {
  # data processing
  if (isS4(data)) {
    data <- data@result %>% tidyr::drop_na()
  } else if (is.data.frame(data)) {
    data <- data %>% tidyr::drop_na()
  } else {
    print("The data format must be S4 object or data frame.")
  }

  if (all(c("ID", "Description", "geneID") %in% colnames(data))) {
    path <- separate_rows(data[1:top, ], geneID, sep = "/")
    if (label.name == "ID") {
      kegg.df <- path %>%
        dplyr::select(ID, geneID) %>%
        dplyr::mutate(frequen = as.numeric("1")) %>%
        dplyr::distinct()
    } else if (label.name == "Description") {
      kegg.df <- path %>%
        dplyr::select(Description, geneID) %>%
        dplyr::mutate(frequen = as.numeric("1")) %>%
        dplyr::distinct()
    } else {
      print("The label.name is 'ID' or 'Description'. ")
    }
  } else {
    path <- separate_rows(data[1:top, ], 2, sep = "/")
    colnames(path)[1:2] <- c(label.name, "geneID")
    kegg.df <- path %>%
      dplyr::mutate(frequen = as.numeric("1")) %>%
      dplyr::distinct()
  }

  # use the selected terms to build the data format
  se_index <- c(ifelse(label.name == "ID", "ID", "Description"), "geneID")
  value <- tail(colnames(kegg.df), 1)
  # build a data frame of nodes
  if (length(se_index) < 2) {
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(se_index), function(i) {
      dots <- se_index[1:i]
      kegg.df %>%
        dplyr::group_by(.dots = dots) %>%
        dplyr::summarise(
          node.size = sum(.data[[value]]),
          node.level = se_index[[i]],
          node.count = n(),
          .groups = 'drop'
        ) %>%
        dplyr::mutate(
          node.short_name = as.character(.data[[dots[[length(dots)]]]]),
          node.branch = as.character(.data[[dots[[1]]]])
        ) %>%
        tidyr::unite(node.name, dots, sep = "/")
    })
    newdata <- do.call("rbind", list) %>% as.data.frame()
    newdata$node.level <- factor(newdata$node.level, levels = se_index)
    if (is.null(root)) {
      nodes_kegg <- newdata
    } else {
      root_data <- data.frame(
        node.name = root,
        node.size = sum(kegg.df[[value]]),
        node.level = root,
        node.count = 1,
        node.short_name = root,
        node.branch = root,
        stringsAsFactors = F
      )
      newdata <- rbind(root_data, newdata)
      newdata$node.level <- factor(newdata$node.level, levels = c(root, se_index))
      nodes_kegg <- newdata
    }
  }
  # build a data frame of edges
  if (length(se_index) < 2) {
    stop("please specify at least two index column(s)")
  } else if (length(se_index) == 2) {
    newdata2 <- kegg.df %>%
      dplyr::mutate(from = .data[[se_index[[1]]]]) %>%
      tidyr::unite(., to, se_index, sep = "/") %>%
      dplyr::select(., from, to) %>%
      dplyr::mutate_at(., c("from", "to"), as.character)
  } else {
    list <- lapply(seq(2, length(se_index)), function(i) {
      dots <- se_index[1:i]
      kegg.df %>%
        tidyr::unite(from, dots[-length(dots)], sep = "/", remove = FALSE) %>%
        tidyr::unite(., to, dots, sep = "/") %>%
        dplyr::select(., from, to) %>%
        dplyr::mutate_at(., c("from", "to"), as.character)
    })
    newdata2 <- do.call("rbind", list) %>% as.data.frame(newdata2)
  }

  if (is.null(root)) {
    edges_kegg <- newdata2
  } else {
    root_data <- kegg.df %>%
      dplyr::group_by(.dots = se_index[[1]]) %>%
      dplyr::summarise(count = n(), .groups = 'drop') %>%
      dplyr::mutate(from = root, to = as.character(.data[[se_index[[1]]]])) %>%
      dplyr::select(., from, to)
    edges_kegg <- rbind(root_data, newdata2)
  }
  graph <- tidygraph::tbl_graph(nodes_kegg, edges_kegg)
  # start drawing
  ggraph(graph, layout = "dendrogram", circular = TRUE) +
    ggraph::geom_edge_diagonal(aes(color = node1.node.branch),
      alpha = color.alpha
    ) +
    ggraph::geom_node_point(aes(size = node.size, color = node.branch),
      alpha = color.alpha
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    scale_size(range = c(0.5, 30)) +
    geom_node_text(
      aes(
        x = 1.0175 * x,
        y = 1.0175 * y,
        label = node.short_name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = leaf, color = node.branch
      ),
      size = text.size[1], hjust = "outward"
    ) +
    scale_colour_manual(values = rep(
      RColorBrewer::brewer.pal(8, color.node),
      ifelse(label.name == "ID",
        length(unique(path$ID)),
        length(unique(path$"Description"))
      )
    )) +
    # add inner circle text label
    geom_node_text(
      show.legend = FALSE,
      aes(
        label = node.short_name,
        filter = !leaf,
        color = node.branch
      ),
      fontface = "bold",
      size = text.size[2]
    ) +
    xlim(-1.2, 1.2)+
    ylim(-1.2, 1.2)
}
