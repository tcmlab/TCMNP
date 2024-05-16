#' @title Drug composition in prescription
#'
#' @param data data must consist of herb and weight columns
#' @param start.degree the angle at which the circle plot starts
#' @param first.color first circle color
#' @param second.color second circle color
#' @param third.color third circle color
#' @param fourth.color fouth circle color
#' @param fifth.color fifth circle color
#' see "RColorBrewer::display.brewer.all()"
#' @param radius.ratio from the outside to the inside,
#' the size of the proportion of each layer of the circle plot
#' @param text.size the font size of each layer circle
#'
#' @return figure
#' @export
#' @importFrom circlize circos.par
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.trackPlotRegion
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom circlize highlight.sector
#' @importFrom circlize colorRamp2
#' @importFrom circlize circos.clear
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom utils data
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' herb = c("ma huang", "ku xing ren", "hua ju hong",
#'            "cang zhu", "guang huo xiang", "yi yi ren",
#'            "hu zhang", "qing hao", "ma bian cao", "lu gen",
#'            "ting li zi", "shi gao", "gan cao")
#' xfbdf.comp = data.frame( herb = herb,
#'   weight = c(6, 15, 15, 10, 15, 30, 20, 12, 30, 30, 15, 30, 10),
#'   property = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Property,
#'   flavor = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Flavor,
#'   meridian = herb_pm[match(herb, herb_pm$Herb_name_pinyin), ]$Meridian
#' )
#' tcm_compound(xfbdf.comp)
#' }
tcm_compound <- function(data,
                          start.degree = 180,
                          first.color = "Spectral",
                          second.color = "RdBu",
                          third.color = c(
                            "#217AB3", "#A6CEE3", "#5EB9A8", "#33A02C",
                            "#FFFF99", "#FDBF6F", "#FB9A99", "#EF5B5B",
                            "#E31A1C", "lightgrey"
                          ),
                          fourth.color = "Spectral",
                          fifth.color = "Set2",
                          radius.ratio = c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                          text.size = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) {
  if (is.data.frame(data)) {
    data$label <- paste0(data$weight, " g")
    # Set to take two digits after the decimal point
    data$ratio <- paste0(base::format(data$weight / sum(data$weight) * 100,
      digits = 2, nsmall = 2
    ), "%") %>%
      as.character()
  } else {
    print("The data must be a data frame.")
  }
  # color settings
  plot.color <- grDevices::colorRampPalette(
    if (nrow(data) <= 8) {
      RColorBrewer::brewer.pal(nrow(data), first.color)
    } else {
      RColorBrewer::brewer.pal(8, first.color)
    }
  )(nrow(data))
  names(plot.color) <- data$herb
  data$number <- rep(1, nrow(data))
  sector_lim <- c(
    rep(0, nrow(data)),
    data$weight / sum(data$weight)
  ) %>% matrix(ncol = 2)
  rownames(sector_lim) <- data$herb
  # start drawing
  circlize::circos.clear()

  # The first circle: the name of each drug
  circlize::circos.par(
    "start.degree" = start.degree,
    points.overflow.warning = FALSE
  ) # set the starting angle
  circlize::circos.initialize(
    factors = data$herb,
    xlim = sector_lim
  )
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    factors = data$herb,
    track.height = radius.ratio[1],
    bg.col = plot.color,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        mean(xlim),
        mean(ylim),
        sector.index,
        cex = text.size[1],
        facing = "bending.inside",
        niceFacing = TRUE
      )
    }
  )
  # second circle: flavor of each drug
  flavor.number <- length(unique(data$flavor)) %>% as.numeric()
  flavor_color <- grDevices::colorRampPalette(
    if (flavor.number <= 8) {
      RColorBrewer::brewer.pal(flavor.number, second.color)
    } else {
      RColorBrewer::brewer.pal(8, second.color)
    }
  )(flavor.number)
  names(flavor_color) <- unique(data$flavor)

  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 2,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
    }, track.height = radius.ratio[2]
  )

  list1 <- list()
  for (i in seq_along(data$flavor)) {
    list1[[i]] <- str_split(data$flavor, ",")[[i]][1]
  }
  flavor.label1 <- unlist(list1)

  list2 <- list()
  for (i in seq_along(data$flavor)) {
    list2[[i]] <- str_split(data$flavor, ",")[[i]][2]
  }
  flavor.label2 <- unlist(list2)

  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 2,
      sector.index = (data$herb)[i],
      col = paste0((flavor_color[match(data$flavor, unique(data$flavor))])[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.outside",
      cex = text.size[2],
      adj = c(0.5, 2),
      labels = (flavor.label1)[i],
      sector.index = (data$herb)[i]
    )
  }

  for (i in seq_len(nrow(data))) {
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.outside",
      cex = text.size[2],
      adj = c(0.5, 0),
      labels = (flavor.label2)[i],
      sector.index = (data$herb)[i]
    )
  }

  # third circle: property of each drug
  property <- c("great cold", "cold", "mildly cold", "cool",
                "even", "mildly warm", "warm",
                "hot", "great hot", "others")
  property.color<-third.color
  names(property.color) <- property

  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 3,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
    }, track.height = radius.ratio[3]
  )
  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 3,
      sector.index = (data$herb)[i],
      col = paste0((property.color[match(data$property, property)])[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.inside",
      cex = text.size[3],
      adj = c(0.5, 0.5),
      labels = (data$property)[i],
      sector.index = (data$herb)[i]
    )
  }

  # fourth circle: meridian of each drug
  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 4,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
    }, track.height = radius.ratio[4]
  )

  list1 <- list()
  for (i in seq_along(data$meridian)) {
    list1[[i]] <- str_split(data$meridian, ",")[[i]][1]
  }
  meridian.label1 <- unlist(list1)

  list2 <- list()
  for (i in seq_along(data$meridian)) {
    list2[[i]] <- str_split(data$meridian, ",")[[i]][2]
  }
  meridian.label2 <- unlist(list2)

  list3 <- list()
  for (i in seq_along(data$meridian)) {
    list3[[i]] <- str_split(data$meridian, ",")[[i]][3]
  }
  meridian.label3 <- unlist(list3)


  meridian.number <- length(unique(meridian.label1)) %>% as.numeric()
  meridian.color <- grDevices::colorRampPalette(
    if (meridian.number <= 8) {
      RColorBrewer::brewer.pal(meridian.number, fourth.color)
    } else {
      RColorBrewer::brewer.pal(8, fourth.color)
    }
  )(meridian.number)
  names(meridian.color) <- unique(meridian.label1)

  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 4,
      sector.index = (data$herb)[i],
      col = paste0((meridian.color[match(meridian.label1, unique(meridian.label1))])[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.outside",
      cex = text.size[4],
      adj = c(0.5, 2.0),
      labels = (meridian.label1)[i],
      sector.index = (data$herb)[i]
    )
  }

  for (i in seq_len(nrow(data))) {
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.outside",
      cex = text.size[4],
      adj = c(0.5, 0.5),
      labels = (meridian.label2)[i],
      sector.index = (data$herb)[i]
    )
  }

  for (i in seq_len(nrow(data))) {
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.outside",
      cex = text.size[4],
      adj = c(0.5, -1),
      labels = (meridian.label3)[i],
      sector.index = (data$herb)[i]
    )
  }

  # fifth circle: dosage of each drug
  color.num <- length(unique(sort(as.numeric(data$weight)))) %>% as.numeric()
  dose.color <- grDevices::colorRampPalette(
    if (color.num <= 8) {
      RColorBrewer::brewer.pal(color.num, fifth.color)
    } else {
      RColorBrewer::brewer.pal(8, fifth.color)
    }
  )(color.num)
  dose.weight <- unique(sort(as.numeric(data$weight), decreasing = TRUE))
  names(dose.color) <- dose.weight

  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 5,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
    }, track.height = radius.ratio[5]
  )
  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 5,
      sector.index = (data$herb)[i],
      col = paste0((dose.color[match(data$weight, dose.weight)])[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "bending.inside",
      cex = text.size[5],
      adj = c(0.5, 0.5),
      labels = (data$label)[i],
      sector.index = (data$herb)[i]
    )
  }
  # sixth circle: the percentage of each drug dose in the total dose
  circlize::circos.trackPlotRegion(
    factors = data$herb,
    track.index = 6,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
    }, track.height = radius.ratio[6]
  )
  for (i in seq_len(nrow(data))) {
    highlight.sector(
      track.index = 6,
      sector.index = (data$herb)[i],
      col = paste0((dose.color[match(data$weight, dose.weight)])[i])
    )
    xlim <- get.cell.meta.data("xlim", sector.index = (data$herb)[i])
    ylim <- get.cell.meta.data("ylim", sector.index = (data$herb)[i])
    circos.text(
      mean(xlim),
      mean(ylim),
      niceFacing = TRUE,
      facing = "clockwise",
      cex = text.size[6],
      adj = c(0.5, 0.5),
      labels = (data$ratio)[i],
      sector.index = (data$herb)[i]
    )
  }
  circlize::circos.clear()
}
