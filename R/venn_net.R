#' venn network
#'
#' @param data data frame
#' @param node.color  node color: see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size markup text size
#' @param label.degree The node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' venn.color1 <- c("#b1d1e7", "#f0a4af", "#e5dfae", "#a8e1da", "#f8d068")
#' venn.color2 <- c("#d3e2b7", "#fbca50", "#c9d5d4", "#baa28a", "#4caecc")
#' venn.color3 <- c("#C1BCBF", "#B1D0A9", "#FFCCC3", "#FFD5AB", "#BCCBE5")
#' venn.color4 <- c("#e4c9b2", "#92C2DD", "#f49d98", "#fcd68f", "#629076")
#' venn.color5 <- c("#FF8748", "#5BAA56", "#B8BB5B", "#4186B7", "#8679BE")
#' venn.color6 <- c("#e6c09e", "#0d5888", "#cb8b3e", "#9cd6d6", "#dbcb09")
#' venn.color7 <- c("#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#75cbdc")
#' venn.color8 <- c("#604192", "#139177", "#ed9e08", "#f56f08", "#4caecc")
#' venn.color9 <- c("#d37a20", "#dbcb09", "#3a9cbc", "#dd7208", "#a30019")
#' @param edge.alpha edge alpha
#' @param edge.type  "arc", "fan", "hive"
#' @param edge.width edge width
#' @param graph.layout Network Diagram Layout:
#' "kk", "nicely", "sphere", "circle",
#' "bipartite", "star","tree", "randomly",
#' "gem", "graphopt","lgl", "grid",
#' "mds", "sugiyama","fr"
#' @param label.repel label repel: "TRUE" or "FALSE"

#' @return venn network diagrams
#' @export
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 aes
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_arc
#' @importFrom ggraph geom_edge_fan
#' @importFrom ggraph geom_edge_hive
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph scale_edge_color_manual
#' @importFrom ggraph theme_graph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph V
#' @importFrom igraph E
#' @examples
#' \dontrun{
#' data(venn_data, package = "TCMNP")
#' set.seed(1234)
#' gene <- names(sort(table(venn_data$gene), decreasing = TRUE))[1:50]
#' data <- venn_data[venn.data$gene %in% gene, ]
#' data2 <- dplyr::sample_n(venn_data, 100) %>%
#'   rbind(data)
#' venn_net(data2, label.degree = 1,graph.layout = "fr")
#' venn_net(data2, edge.type = "hive", label.degree = 4)
#' }
venn_net <- function(data,
                     node.color = "Spectral",
                     node.size = c(1, 10),
                     label.size = 3,
                     label.degree = 1,
                     label.repel = TRUE,
                     edge.color = c("#5BAA56", "#f49d98", "#4186B7", "#8679BE"),
                     edge.alpha = 0.5,
                     edge.width = 0.5,
                     edge.type = "arc",
                     graph.layout = "kk") {
  # node data
  data <- data %>% dplyr::distinct()
  if (is.data.frame(data)) {
    if (length(unique(data[, 1])) > length(unique(data[, 2]))) {
      colnames(data) <- c("to", "from")
    } else {
      colnames(data) <- c("from", "to")
    }
  } else {
    print("The data must be a data frame with two columns.")
  }
  nodes <- data %>%
    {
      data.frame(gene = c(.$from, .$to))
    } %>%
    dplyr::distinct()

  # Create a network graph from data and nodes
  net <- igraph::graph_from_data_frame(
    d = data,
    vertices = nodes,
    directed = FALSE
  )
  # V and E are functions of the igraph package,
  # which are used to modify the nodes (nodes) and
  # connections (data) of the network graph respectively.
  igraph::V(net)$degree <- igraph::degree(net)
  igraph::V(net)$size <- igraph::degree(net)

  if (edge.type == "arc") {
    p <- ggraph(net, layout = graph.layout) +
      geom_edge_arc(
        aes(
          edge_width = edge.width,
          colour = as.character(data$from),
          edge_alpha = edge.alpha
        ),
        show.legend = FALSE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, node.color))
      ) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = label.repel
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans") +
      scale_edge_color_manual(values = edge.color)
    return(p)
  } else if (edge.type == "fan") {
    p <- ggraph(net, layout = graph.layout) +
      geom_edge_fan(
        aes(
          edge_width = edge.width,
          colour = as.character(data$from),
          edge_alpha = edge.alpha
        ),
        show.legend = FALSE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, node.color))
      ) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = label.repel
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans") +
      scale_edge_color_manual(edge.color)
    return(p)
  } else if (edge.type == "hive") {
    p <- ggraph(net, layout = graph.layout) +
      geom_edge_hive(
        aes(
          edge_width = edge.width,
          colour = as.character(data$from),
          edge_alpha = edge.alpha
        ),
        show.legend = FALSE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, node.color))
      ) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = label.repel
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans") +
      scale_edge_color_manual(values = edge.color)
    return(p)
  }
}
