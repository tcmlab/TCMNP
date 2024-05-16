#' Network interaction diagram of herbs, ingredients, and targets
#'
#' @param network.data data frame
#' must contain herb, molecule, target three columns of data
#' @param node.color node color
#' see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size label size
#' @param label.degree
#' the node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' @param edge.width edge width
#' @param graph.layout etwork Diagram Layout:
#' @param graph.layout Network Diagram Layout:
#' "kk", "nicely", "circle", "sphere",
#' "bipartite", "star", "tree", "randomly",
#' "gem", "graphopt","lgl", "grid",
#' "mds", "sugiyama","fr"
#' @param rem.dis.inter remove single free unconnected nodes
#' @param label.repel label repel
#'
#' @return Network Diagram
#' @export
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 aes
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr count
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_fan
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph scale_edge_width
#' @importFrom ggraph geom_edge_link0
#' @importFrom ggraph theme_graph
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph graph_from_data_frame
#' @examples
#' \dontrun{
#' data("xfbdf", package = "TCMNP")
#' network.data <- xfbdf %>%
#'   dplyr::select(herb, molecule, target) %>%
#'   sample_n(100, replace = FALSE) %>%
#'   as.data.frame()
#' tcm_net(network.data, node.color= "Spectral",
#' label.degree = 1, rem.dis.inter = TRUE,graph.layout = "fr",
#' label.size = 3)
#' }
tcm_net <- function(network.data,
                    node.color = "RdBu",
                    node.size = c(2, 8),
                    label.size = 4,
                    label.degree = 3,
                    label.repel = TRUE,
                    edge.color = "lightgrey",
                    edge.width = c(0.2, 2),
                    graph.layout = "kk",
                    rem.dis.inter = FALSE) {
  network.data<-as.data.frame(network.data)
  links <- rbind(
    network.data %>%
      dplyr::select(herb, molecule) %>%
      dplyr::rename(from = herb, to = molecule) %>%
      dplyr::mutate(weight = 1),
    network.data %>%
      dplyr::select(molecule, target) %>%
      dplyr::rename(from = molecule, to = target) %>%
      dplyr::mutate(weight = 1)
  ) %>%
    dplyr::distinct()
  # Create a network graph from links and nodes
  if (rem.dis.inter == FALSE) {
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      dplyr::distinct()
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    # ggraph is a package based on ggplot2,
    # the syntax is similar to regular ggplot2
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight), edge_colour = edge.color) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree,
        size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      geom_node_text(aes(
        filter = degree >= label.degree,
        label = name
      ), size = label.size, repel = label.repel) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      ggraph::theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE) {
    links <- links %>%
      dplyr::mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      dplyr::mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      dplyr::distinct()
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    # Create a network graph from links and nodes
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight),
        edge_colour = edge.color
      ) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree, size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = label.repel
      ) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(
        name = "degree",
        range = node.size
      ) +
      ggraph::theme_graph(base_family = "sans")
  }
}
