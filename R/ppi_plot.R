#' Protein-protein interaction analysis
#'
#' @param data Human protein-protein interaction (PPI)
#' data were downloaded from the STRING database.
#' @param node.color  node color: see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size markup text size
#' @param label.degree The node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' @param edge.width edge width
#' @param rem.dis.inter remove single free unconnected nodes
#' @param graph.layout Network Diagram Layout:
#' "kk", "nicely", "sphere", "circle",
#' "bipartite", "star","tree", "randomly",
#' "gem", "graphopt","lgl", "grid",
#' "mds", "sugiyama"
#' @param label.repel label repel

#' @return network diagrams
#' @export
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_linewidth_continuous
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr count
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_fan
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph scale_edge_width
#' @importFrom ggraph theme_graph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph V
#' @importFrom igraph E

#' @examples
#' \dontrun{
#' data(ppi_data, package = "TCMNP")
#' ppi_plot(ppi_data, label.repel = FALSE)

#' }
ppi_plot <- function(data,
                     node.color = "RdBu",
                     node.size = c(1, 10),
                     label.size = 4,
                     label.degree = 5,
                     label.repel = TRUE,
                     edge.color = "lightgrey",
                     edge.width = c(0.2, 2),
                     rem.dis.inter = FALSE,
                     graph.layout = "kk") {
 # data processing
  if (is.data.frame(data)) {
    data <- data %>% dplyr::distinct()
    if (length(unique(data[, 1])) > length(unique(data[, 2]))) {
      colnames(data) <- c("to", "from", "weight")
    } else {
      colnames(data) <- c("from", "to", "weight")
    }
  } else {
    print("The data must be a data frame with three columns.")
  }
  links <- data
  # node data
  nodes <- links %>%
    {
      data.frame(gene = c(.$from, .$to))
    } %>%
    dplyr::distinct()

  # Create a network graph from links and nodes
  if (rem.dis.inter == FALSE) {
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight
    # Plot with ggraph
    ggraph(net, layout = graph.layout) +
      geom_edge_fan(
        aes(
          edge_width = score,
          colour = score
        ),
        color = edge.color,
        show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, node.color))
      ) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size,
        repel = label.repel
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_linewidth_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE) {
    # remove dissociative interactions
    # If the from of a link in the links data frame
    # only appears once, and the to appears only once, remove it
    links_2 <- links %>%
      dplyr::mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      dplyr::mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    # new node data
    nodes_2 <- links_2 %>%
      {
        data.frame(gene = c(.$from, .$to))
      } %>%
      dplyr::distinct()
    # create a network diagram
    net_2 <- igraph::graph_from_data_frame(
      d = links_2,
      vertices = nodes_2,
      directed = FALSE
    )
    # add the necessary parameters
    igraph::V(net_2)$degree <- igraph::degree(net_2)
    igraph::V(net_2)$size <- igraph::degree(net_2)
    igraph::E(net_2)$score <- igraph::E(net_2)$weight
    ggraph(net_2, layout = graph.layout) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(colours = rev(
        RColorBrewer::brewer.pal(8, node.color)
      )) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size, repel = label.repel
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_linewidth_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  }
}
