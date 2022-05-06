#' Write code for loading libraries
#' @param lib vector of libraries
#' @rdname wr_lib
#' @export wr_lib
#'
wr_lib <- function(lib) {
  oup <- ""
  for (iLib in lib) {
    oup <- paste0(oup, "library(", iLib, ")\n")
  }
  glue::glue(paste0(oup,"\n"))
}

#' Write code for font
#' @param font Google font name
#' @rdname wr_font
#' @export wr_font
#'
wr_font <- function(font = "Lato") {
  paste0(
'
# load font for plot
sysfonts::font_add_google(name = "',font,'", family = "',font,'")
showtext::showtext_auto()

'
  )
}

#' Write code for loading objects for ui.R and server.R
#' @param prefix file prefix
#' @param tabs Vector of tab names to include
#' @rdname wr_load
#' @export wr_load
#'
wr_load <- function(prefix, tabs) {
x <- paste0('
if(!exists("{prefix}conf")) {prefix}conf = readRDS("{prefix}conf.rds")
if(!exists("{prefix}def")) {prefix}def  = readRDS("{prefix}def.rds")
if(!exists("{prefix}gene")) {prefix}gene = readRDS("{prefix}gene.rds")
if(!exists("{prefix}meta")) {prefix}meta = readRDS("{prefix}meta.rds")
')
if("mar" %in% tabs) x <- paste0(x,'if(!exists("{prefix}mar")) {prefix}mar = readRDS("{prefix}mar.rds")\n')
glue::glue(x,"\n")
}

#' Write code for fixed portion of server.R
#' @param font (Character) Google font
#' @rdname wr_sv_fix
#' @export wr_sv_fix
#'
wr_sv_fix <- function(font = "Lato") {
  glue::glue('

### Initialise variables and functions ----

# Colour palette
cList <- list(
  c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
  ),
  c(
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027"
  )[c(1, 1:9, 9)],
  c(
    "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
    "#2C728E", "#3B528B", "#472D7B", "#440154"
  )
)
names(cList) <- c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

# Panel sizes
pList <- c("400px", "600px", "800px")
names(pList) <- c("Small", "Medium", "Large")
pList2 <- c("500px", "700px", "900px")
names(pList2) <- c("Small", "Medium", "Large")
pList3 <- c("600px", "800px", "1000px")
names(pList3) <- c("Small", "Medium", "Large")
# baseplot font size
sList <- c(12, 14, 18, 22)
names(sList) <- c("Smaller", "Small", "Medium", "Large")
# ggrepel font size
lList <- c(5, 6, 7)
names(lList) <- c("Small", "Medium", "Large")

# Function to extract legend
g_legend <- function(a.gplot) {{
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}}

# progress indicator
show_progress <- function(...){{
  return(shinycssloaders::withSpinner(..., type = 7, color = "#95a5a6"))
}}

# Plot theme
# @description Custom ggplot theme
# @param base_size (Numeric) Base font size
# @param XYval (Logical) Show XY axes text?
# @param Xang (Numeric) X axis text angle
# @param XjusH (Numeric) X axis horizontal justification
# @param lpos (Character) Position of Legend
# @param font (Character) Google font
# @param col_text (Character) Text colour
# @param col_line (Character) Line colour
#
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5, lpos = "bottom", font = "{font}", col_text = "grey30", col_line = "grey60") {{
  oupTheme <- theme(
    text = element_text(size = base_size, family = font),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(colour = col_line),
    axis.ticks = element_line(colour = col_line),
    axis.title = element_text(colour = col_text),
    axis.text = element_text(size = base_size, colour = col_text),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    strip.background = element_rect(colour = "white"),
    legend.position = lpos,
    legend.key = element_rect(colour = NA, fill = NA)
  )
  if (!XYval) {{
    oupTheme <- oupTheme +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
  }}
  
  return(oupTheme)
}}

### Plotting functions ----

# @description DR scatterplot for gene expression
# @param dtab (Data.table) with columns X, Y, geneName and val
# @param bgCells (Logical) Background points
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scScatter <- function(dtab, bgCells = FALSE, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt){{
  
  if(any(!c("X", "Y", "geneName", "val") %in% colnames(dtab))) "Input missing one or more columns: X, Y, geneName, val."
  
  ggOut <- ggplot(dtab, aes(X, Y, color = val))
  rat <- (max(dtab$X) - min(dtab$X)) / (max(dtab$Y) - min(dtab$Y))
  ltitle <- dtab$geneName[1]
  
  if (bgCells) {{
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }}
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    scale_color_gradientn(ltitle, colours = cList[[inpcol]]) +
    #guides(color = guide_colorbar(barwidth = 20)) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt, lpos = "right")
  
  if (inpasp == "Square") {{
    ggOut <- ggOut + coord_fixed(ratio = rat)
  }} else if (inpasp == "Fixed") {{
    ggOut <- ggOut + coord_fixed()
  }}
  
  return(ggOut)
}}

# @description Plot gene expression on dimred for any number of input genes
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to expression h5
# @param inpGene (Numeric) Named gene expression vector
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inpncol (Integer) Number of rows of plots
# @details 
# Config table contains columns ID (Character, columns name in metadata), UI (Character, UI id), fID (Character, Levels for categorical data, | separated), fCL (Character, Colours for categorical data, | separated), fRow (Integer, number of rows), grp (Logical), dimred (Logical)
#
scFeature <- function(inpConf, inpMeta, inpdrX, inpdrY, inp, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inpncol = 0){{
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 36, "More than 36 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 0, "Please input at least 1 gene to plot!"))
  
  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {{
    tmp <- inpMeta[, c("sampleID",inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, inpConf[UI == inpsub1]$ID),with = FALSE]
    colnames(tmp) <- c("sampleID", "X", "Y", "sub")
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp))
  }}
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) ggData <- ggData[sub %in% inpsub2]
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }}
  
  if (inpord == "Max") {{
    ggData <- ggData[order(val)]
  }} else if (inpord == "Min") {{
    ggData <- ggData[order(-val)]
  }} else if (inpord == "Random") {{
    ggData <- ggData[sample(nrow(ggData))]
  }}
  
  ggDataSplit <- split(ggData, by=c("geneName"), flatten=FALSE)
  plist <- vector("list", length = length(ggDataSplit))
  for(i in seq_along(ggDataSplit)){{
    plist[[i]] <- scScatter(dtab = ggDataSplit[[i]], bgCells, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt)
  }}

  if(inpncol < 1) inpncol <- floor(sqrt(length(plist)))
  ggOut <- wrap_plots(plist) + plot_layout(ncol = inpncol)
  return(ggOut)
}}

# @description Plot cell information on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inplab 
#
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "val", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }}

  if (inpord == "Max") {{
    ggData <- ggData[order(val)]
  }} else if (inpord == "Min") {{
    ggData <- ggData[order(-val)]
  }} else if (inpord == "Random") {{
    ggData <- ggData[sample(nrow(ggData))]
  }}

  # Do factoring if required
  if (!is.na(inpConf[UI == inp1]$fCL)) {{
    ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\\\|")[[1]]
    names(ggCol) <- levels(ggData$val)
    ggLvl <- levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    ggData$val <- factor(ggData$val, levels = ggLvl)
    ggCol <- ggCol[ggLvl]
  }}

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) {{
  ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }}
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) +
    xlab(inpdrX) +
    ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)

  if (is.na(inpConf[UI == inp1]$fCL)) {{
    ggOut <- ggOut +
    scale_color_gradientn("", colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20)) 
  }} else {{
    sListX <- min(nchar(paste0(levels(ggData$val), collapse = "")), 200)
    sListX <- 0.75 * (sList - (1.5 * floor(sListX / 50)))
    ggOut <- ggOut + scale_color_manual("", values = ggCol) +
      guides(color = guide_legend(
        override.aes = list(size = 5),
        nrow = inpConf[UI == inp1]$fRow
      )) +
      theme(legend.text = element_text(size = sListX[inpfsz]))

    if (inplab) {{ ggData3 <- ggData[, .(X = mean(X), Y = mean(Y)), by = "val"]
      lListX <- min(nchar(paste0(ggData3$val, collapse = "")), 200)
      lListX <- lList - (0.25 * floor(lListX / 50))
      ggOut <- ggOut +
        geom_text_repel(
          data = ggData3, aes(X, Y, label = val),
          color = "grey10", bg.color = "grey95", bg.r = 0.15,
          size = lListX[inpfsz], seed = 42
        ) }}
    }}

  if (inpasp == "Square") {{ ggOut <- ggOut + coord_fixed(ratio = rat) }} else if (inpasp == "Fixed") {{ ggOut <- ggOut + coord_fixed() }}

  if(is.numeric(inpMeta[[inp1]])) ggOut <- ggOut + guides(color = guide_colorbar(barwidth = 25))
  return(ggOut)
}}

# @description Cell number, stats data
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Name of metadata column for coloring points
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Name/path to h5 object
# @param inpGene (Character) Gene to use
# @param inpsplt (Character) Whether to split table by Quartile or Decile
#
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inpsplt) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]

  colnames(ggData) <- c("group", "sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = T)))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{ ggData <- ggData[sub %in% inpsub2] }}

  # Split inp1 if necessary
  if (is.na(inpConf[UI == inp1]$fCL)) {{ if (inpsplt == "Quartile") {{ nBk <- 4 }}
    if (inpsplt == "Decile") {{ nBk <- 10 }}
    ggData$group <- cut(ggData$group, breaks = nBk) }}

  # Actual data.table
  ggData$express <- FALSE
  ggData[val2 > 0]$express <- TRUE
  ggData1 <- ggData[express == TRUE, .(nExpress = .N), by = "group"]
  ggData <- ggData[, .(nCells = .N), by = "group"]
  ggData <- ggData1[ggData, on = "group"]
  ggData <- ggData[, c("group", "nCells", "nExpress"), with = FALSE]
  ggData[is.na(nExpress)]$nExpress <- 0
  ggData$pctExpress <- 100 * ggData$nExpress / ggData$nCells
  ggData <- ggData[order(group)]
  colnames(ggData)[3] <- paste0(colnames(ggData)[3], "_", inp2)

  return(ggData)
}}

# @description Plot gene expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val < 0]$val <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{ bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2] }}

  if (inpord == "Max") {{
    ggData <- ggData[order(val)]
  }} else if (inpord == "Min") {{ 
    ggData <- ggData[order(-val)] 
  }} else if (inpord == "Random") {{ 
    ggData <- ggData[sample(nrow(ggData))] 
  }}

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) {{ 
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }}

  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") {{ 
    ggOut <- ggOut + coord_fixed(ratio = rat) 
  }} else if (inpasp == "Fixed") {{
    ggOut <- ggOut + coord_fixed()
  }}

  return(ggOut)
}}

# Plot gene coexpression on dimred
bilinear <- function(x, y, xy, Q11, Q21, Q12, Q22) {{
  oup <- (xy - x) * (xy - y) * Q11 + x * (xy - y) * Q21 + (xy - x) * y * Q12 + x * y * Q22
  oup <- oup / (xy * xy)
  return(oup)
}}

# @description Gene Co-expression on dimred with legend
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoexFull <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {{
  g1 <- scDRcoex(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                       inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                       inpasp, inptxt)
  g2 <- scDRcoexLeg(inp1, inp2, inpcol, inpfsz)
  g <- wrap_plots(g1, g2) + plot_layout(ncol = 2, nrow = 1, widths = c(10, 2))
  return(g)
}}

# @description Gene Co-expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }}

  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {{
    c10 <- c(255, 0, 0) 
  }} else if (cInp[1] == "Orange (Gene1)") {{
    c10 <- c(255, 140, 0)
  }} else {{
    c10 <- c(0, 255, 0)
  }}

  if (cInp[2] == "Green (Gene2)") {{ c01 <- c(0, 255, 0) }} else {{ c01 <- c(0, 0, 255) }}

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Map colours
  ggData$v1 <- round(nTot * ggData$val1 / max(ggData$val1))
  ggData$v2 <- round(nTot * ggData$val2 / max(ggData$val2))
  ggData$v0 <- ggData$v1 + ggData$v2
  ggData <- gg[ggData, on = c("v1", "v2")]

  if (inpord == "Max") {{
    ggData <- ggData[order(v0)]
  }} else if (inpord == "Min") {{
    ggData <- ggData[order(-v0)]
  }} else if (inpord == "Random") {{
    ggData <- ggData[sample(nrow(ggData))]
  }}

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y))
  if (bgCells) {{
    ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }}
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20, color = ggData$cMix) +
    xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[1]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") {{
    ggOut <- ggOut + coord_fixed(ratio = rat)
  }} else if (inpasp == "Fixed") {{
    ggOut <- ggOut + coord_fixed()
  }}
  
  return(ggOut)
}}
  
# Co-exp plot legend
# @param inp1
# @param inp2
# @param inpcol Colour
# @param inpfsz Font size
#
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz) {{
  
  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {{
    c10 <- c(255, 0, 0)
  }} else if (cInp[1] == "Orange (Gene1)") {{ 
    c10 <- c(255, 140, 0) 
  }} else {{ 
    c10 <- c(0, 255, 0) 
  }}

  if (cInp[2] == "Green (Gene2)") {{ 
    c01 <- c(0, 255, 0) 
  }} else {{ 
    c01 <- c(0, 0, 255) 
  }}

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Actual ggplot
  ggOut <- ggplot(gg, aes(v1, v2)) +
    geom_tile(fill = gg$cMix) +
    xlab(inp1) +
    ylab(inp2) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    sctheme(base_size = sList[inpfsz], XYval = TRUE)
  return(ggOut)
}}


scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene) {{
  
  if (is.null(inpsub1)) {{
  inpsub1 <- inpConf$UI[1]
  }}

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{ ggData <- ggData[sub %in% inpsub2] }}

  # Actual data.table
  ggData$express <- "none"
  ggData[val1 > 0]$express <- inp1
  ggData[val2 > 0]$express <- inp2
  ggData[val1 > 0 & val2 > 0]$express <- "both"
  ggData$express <- factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
  ggData <- ggData[, .(nCells = .N), by = "express"]
  ggData$percent <- 100 * ggData$nCells / sum(ggData$nCells)
  ggData <- ggData[order(express)]
  colnames(ggData)[1] <- "expression > 0"
  
  return(ggData)
}}

# Plot violin / boxplot / lineplot
# @description Violin plot gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) X axis cell info
# @param inp2 (Character) Y axis cell info / gene
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inptyp (Character) Plot type. "violin", "boxplot" or "lineplot".
# @param inppts (Logical) Should points be displayed?
# @param inpsiz (Numeric) Point size
# @param inpfsz (Character) Custom font size
# @param inpbarsz (Numeric) Bar size for lineplot
#
scVioBox <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inptyp, inppts, inpsiz, inpfsz, inpbarsz) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}
  if (is.null(inpbarsz)) inpbarsz <- 0.3

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]
  colnames(ggData) <- c("X", "sub")

  # Load in either cell meta or gene expr
  if (inp2 %in% inpConf$UI) {{
    ggData$val <- inpMeta[[inpConf[UI == inp2]$ID]]
  }} else {{
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
    ggData[val < 0]$val <- 0
    set.seed(42)
    tmpNoise <- rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val <- ggData$val + tmpNoise
    h5file$close_all()
  }}

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{
    ggData <- ggData[sub %in% inpsub2]
  }}

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\\\|")[[1]]
  names(ggCol) <- levels(ggData$X)
  ggLvl <- levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  ggData$X <- factor(ggData$X, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, val, fill = X))
  
  if (inptyp == "violin") {{
    ggOut <- ggOut + geom_violin(scale = "width")
  }} else if (inptyp == "boxplot")  {{
    ggOut <- ggOut + geom_boxplot()
  }} else if (inptyp == "lineplot") {{
    ggOut <- ggOut + geom_col(aes(col = X), position = position_dodge2(preserve = "single"), size = inpbarsz)
  }}
  
  if (inppts) ggOut <- ggOut + geom_jitter(size = inpsiz, shape = 20, alpha = 0.4)

  ggOut <- ggOut + xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_colour_manual("", values = ggCol) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
    
  return(ggOut)
  }}

# Plot proportion plot
# @description Proportion barplot
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inptyp (Character) Plot type. "violin" else boxplot.
# @param inpflp (Logical) Flip coordinates?
# @param inpfsz (Character) Custom font size
#
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inptyp, inpflp, inpfsz) {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "grp", "sub")

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{ ggData <- ggData[sub %in% inpsub2] }}
  ggData <- ggData[, .(nCells = .N), by = c("X", "grp")]
  ggData <- ggData[, {{ tot <- sum(nCells)
    .SD[, .(
      pctCells = 100 * sum(nCells) / tot,
      nCells = nCells
    ), by = "grp"] }}, by = "X"]

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp2]$fCL, "\\\\|")[[1]]
  names(ggCol) <- levels(ggData$grp)
  ggLvl <- levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)]
  ggData$grp <- factor(ggData$grp, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  if (inptyp == "Proportion") {{
    ggOut <- ggplot(ggData, aes(X, pctCells, fill = grp)) +
      geom_col() +
      ylab("Cell Proportion (%)") 
  }} else {{ 
    ggOut <- ggplot(ggData, aes(X, nCells, fill = grp)) +
      geom_col() +
      ylab("Number of Cells")
  }}

  if (inpflp) ggOut <- ggOut + coord_flip()

  ggOut <- ggOut + xlab(inp1) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "right")

  return(ggOut)
}}

# Get gene list
scGeneList <- function(inp, inpGene) {{
  geneList <- data.table(
  gene = inp,
  present = TRUE
  )
  geneList[!gene %in% names(inpGene)]$present <- FALSE
  return(geneList)
}}

# Plot gene expression bubbleplot / heatmap
# @description dotplot / heatmap gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp (Character) Gene names
# @param inpGrp
# @param inpPlt
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpScl
# @param inpRow
# @param inpCol
# @param inpcols
# @param inpfsz (Character) Custom font size
# @param col_line (Character) Line colour
#
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, inpcols, inpfsz, col_line = "grey60") {{
  
  if (is.null(inpsub1)) {{ inpsub1 <- inpConf$UI[1] }}

  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {{
    tmp <- inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) <- c("sampleID", "sub")
    tmp$grpBy <- inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp)) 
  }}
  h5file$close_all()
  
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {{
    ggData <- ggData[sub %in% inpsub2] 
  }}
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate
  ggData$val <- expm1(ggData$val)
  ggData <- ggData[, .(val = mean(val), prop = sum(val > 0) / length(sampleID)),
    by = c("geneName", "grpBy")
  ]
  ggData$val <- log1p(ggData$val)

  # Scale if required
  colRange <- range(ggData$val)
  if (inpScl) {{
    ggData[, val := scale(val), keyby = "geneName"]
    colRange <- c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  }}

  # hclust row/col if necessary
  ggMat <- dcast.data.table(ggData, geneName ~ grpBy, value.var = "val")
  tmp <- ggMat$geneName
  ggMat <- as.matrix(ggMat[, -1])
  rownames(ggMat) <- tmp
  if (inpRow) {{
    hcRow <- dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow <- ggplot() +
      coord_flip() +
      geom_segment(data = hcRow$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$grpBy)),
        labels = unique(ggData$grpBy), expand = c(0, 0)
      ) +
      scale_x_continuous(
        breaks = seq_along(hcRow$labels$label),
        labels = hcRow$labels$label, expand = c(0, 0.5)
      ) +
      sctheme(base_size = sList[inpfsz]) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "white", angle = 45, hjust = 1)
      )
    ggData$geneName <- factor(ggData$geneName, levels = hcRow$labels$label) 
  }} else {{ 
    ggData$geneName <- factor(ggData$geneName, levels = rev(geneList$gene)) 
  }}

  if (inpCol) {{
    hcCol <- dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol <- ggplot() +
      geom_segment(data = hcCol$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_x_continuous(
        breaks = seq_along(hcCol$labels$label),
        labels = hcCol$labels$label, expand = c(0.05, 0)
      ) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$geneName)),
        labels = unique(ggData$geneName), expand = c(0, 0)
      ) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white")
      )
    ggData$grpBy <- factor(ggData$grpBy, levels = hcCol$labels$label)
    }}

  # Actual plot according to plottype
  if (inpPlt == "Bubbleplot") {{
    # Bubbleplot
    ggOut <- ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion",
        range = c(0, 8),
        limits = c(0, 1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)
      ) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank(), legend.box = "vertical") 
  }} else {{
    # Heatmap
    ggOut <- ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank()) 
  }}

  # Final tidy
  ggLeg <- g_legend(ggOut)
  ggOut <- ggOut + theme(legend.position = "none")
  
  if (inpRow & inpCol) {{
    ggOut <- wrap_plots(ggCol, plot_spacer(), ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(1,7,2))
  }} else if (inpRow) {{
    ggOut <- wrap_plots(ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(7,2))
  }} else if (inpCol) {{
    ggOut <- wrap_plots(ggCol, ggOut, ggplotify::as.ggplot(ggLeg))+
               plot_layout(ncol = 1, heights = c(1,7,2))
  }} else {{
    ggOut <- wrap_plots(ggOut, ggplotify::as.ggplot(ggLeg))+
      plot_layout(ncol = 1, heights = c(7,2))
  }}
  
  return(ggOut)
}}


### Server code ----
shinyServer(function(input, output, session) {{

### For all tags and Server-side selectize
observe_helpers()

  '
  )
}

#' Write code for cell info vs gene expression
#'
wr_sv_civge <- function() {
paste0('  
### Tab civge cell info vs gene exp ----

{subst}  output${prefix}_civge_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civge_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_civge_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civge_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_civge_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civge_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_civge_oup1 <- renderPlot({{
  req(input${prefix}_civge_inp1)
  scDRcell({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp1, input${prefix}_civge_sub1, input${prefix}_civge_sub2, input${prefix}_civge_siz, input${prefix}_civge_col1, input${prefix}_civge_ord1, input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt, input${prefix}_civge_lab1)
}})

output${prefix}_civge_oup1.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_civge_oup1", height = pList[input${prefix}_civge_psz]))
}})

output${prefix}_civge_oup1.png <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX, "_", input${prefix}_civge_drY, "_", input${prefix}_civge_inp1, ".png")) }},
 content = function(file) {{
   ggsave(
   file, device = "png", dpi = input${prefix}_civge_oup1.res, bg = "white",
   plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp1,   input${prefix}_civge_sub1, input${prefix}_civge_sub2, input${prefix}_civge_siz, input${prefix}_civge_col1, input${prefix}_civge_ord1,  input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt, input${prefix}_civge_lab1)
   )
}})

output${prefix}_civge_oup1.pdf <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX,"_", input${prefix}_civge_drY,"_", input${prefix}_civge_inp1,".pdf")) }},
 content = function(file) {{
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white",
   plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp1,   input${prefix}_civge_sub1, input${prefix}_civge_sub2, input${prefix}_civge_siz, input${prefix}_civge_col1, input${prefix}_civge_ord1,  input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt, input${prefix}_civge_lab1)
   )
}})

output${prefix}_civge_oup1.svg <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX,"_", input${prefix}_civge_drY,"_", input${prefix}_civge_inp1,".svg")) }},
 content = function(file) {{
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp1,   input${prefix}_civge_sub1, input${prefix}_civge_sub2, input${prefix}_civge_siz, input${prefix}_civge_col1, input${prefix}_civge_ord1,  input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt, input${prefix}_civge_lab1)
   )
}})

output${prefix}_civge_.dt <- renderDataTable({{
 req(input${prefix}_civge_inp2)
 ggData = scDRnum({prefix}conf, {prefix}meta, input${prefix}_civge_inp1, input${prefix}_civge_inp2, input${prefix}_civge_sub1, input${prefix}_civge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
}})

output${prefix}_civge_oup2 <- renderPlot({{
 req(input${prefix}_civge_inp2)
 scDRgene({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp2, input${prefix}_civge_sub1, input${prefix}_civge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_civge_siz, input${prefix}_civge_col2, input${prefix}_civge_ord2, input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt)
}})

output${prefix}_civge_oup2.ui <- renderUI({{
 show_progress(imageOutput("{prefix}_civge_oup2", height = pList[input${prefix}_civge_psz]))
}})

output${prefix}_civge_oup2.png <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX,"_",input${prefix}_civge_drY,"_", input${prefix}_civge_inp2,".png")) }},
 content = function(file) {{
   ggsave(
   file, device = "png", dpi = input${prefix}_civge_oup2.res, bg = "white",
   plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp2, input${prefix}_civge_sub1, input${prefix}_civge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_civge_siz, input${prefix}_civge_col2, input${prefix}_civge_ord2, input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt)
   )
}})

output${prefix}_civge_oup2.pdf <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX,"_",input${prefix}_civge_drY,"_", input${prefix}_civge_inp2,".pdf")) }},
 content = function(file) {{
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white",
   plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp2,  input${prefix}_civge_sub1, input${prefix}_civge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_civge_siz, input${prefix}_civge_col2, input${prefix}_civge_ord2, input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt)
   )
}}) 

output${prefix}_civge_oup2.svg <- downloadHandler(
 filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civge_drX,"_",input${prefix}_civge_drY,"_", input${prefix}_civge_inp2,".svg")) }},
 content = function(file) {{
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_civge_drX, input${prefix}_civge_drY, input${prefix}_civge_inp2,  input${prefix}_civge_sub1, input${prefix}_civge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_civge_siz, input${prefix}_civge_col2, input${prefix}_civge_ord2, input${prefix}_civge_fsz, input${prefix}_civge_asp, input${prefix}_civge_txt)
   )
}}) # End of tab civge

')
}

#' Write code for cell info vs cell info
#'
wr_sv_civci <- function() {
paste0('
### Tab civci cell info vs cell info ----
  
{subst}  output${prefix}_civci_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civci_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_civci_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civci_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_civci_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_civci_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_civci_oup1 <- renderPlot({{
  req(input${prefix}_civci_inp1)
  scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp1, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col1, input${prefix}_civci_ord1, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab1)
}})

output${prefix}_civci_oup1.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_civci_oup1", height = pList[input${prefix}_civci_psz]))
}})

output${prefix}_civci_oup1.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX, "_", input${prefix}_civci_drY, "_", input${prefix}_civci_inp1, ".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", dpi = input${prefix}_civci_oup1.res, bg = "white",
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp1, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col1, input${prefix}_civci_ord1, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab1)
    )
}})

output${prefix}_civci_oup1.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX, "_", input${prefix}_civci_drY, "_", input${prefix}_civci_inp1, ".pdf")) }},
  content = function(file) {{ ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp1, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col1, input${prefix}_civci_ord1, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab1) )
}})

output${prefix}_civci_oup1.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX, "_", input${prefix}_civci_drY, "_", input${prefix}_civci_inp1, ".svg")) }},
  content = function(file) {{ ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp1, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col1, input${prefix}_civci_ord1, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab1) )
}})

output${prefix}_civci_oup2 <- renderPlot({{
  req(input${prefix}_civci_inp2)
  scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp2, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col2, input${prefix}_civci_ord2, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab2)
}})

output${prefix}_civci_oup2.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_civci_oup2", height = pList[input${prefix}_civci_psz]))
}})

output${prefix}_civci_oup2.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX,"_",input${prefix}_civci_drY,"_", input${prefix}_civci_inp2,".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_civci_oup2.res,
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp2, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col2, input${prefix}_civci_ord2, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab2)
    )
}}) 

output${prefix}_civci_oup2.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX,"_",input${prefix}_civci_drY,"_", input${prefix}_civci_inp2,".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp2, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col2, input${prefix}_civci_ord2, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab2) 
    )
}})

output${prefix}_civci_oup2.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_civci_drX,"_",input${prefix}_civci_drY,"_", input${prefix}_civci_inp2,".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell({prefix}conf, {prefix}meta, input${prefix}_civci_drX, input${prefix}_civci_drY, input${prefix}_civci_inp2, input${prefix}_civci_sub1, input${prefix}_civci_sub2, input${prefix}_civci_siz, input${prefix}_civci_col2, input${prefix}_civci_ord2, input${prefix}_civci_fsz, input${prefix}_civci_asp, input${prefix}_civci_txt, input${prefix}_civci_lab2) 
    )
}}) # End of tab civci

')
}

#' Write code for gene exp vs gene exp
#'
wr_sv_gevge <- function() {
paste0('
### Tab gevge gene exp vs gene exp ----

{subst}  output${prefix}_gevge_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gevge_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_gevge_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gevge_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_gevge_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gevge_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_gevge_oup1 <- renderPlot({{
  req(input${prefix}_gevge_inp1)
  scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp1, input${prefix}_gevge_sub1, input${prefix}_gevge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gevge_siz, input${prefix}_gevge_col1, input${prefix}_gevge_ord1, input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
}})

output${prefix}_gevge_oup1.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_gevge_oup1", height = pList[input${prefix}_gevge_psz]))
}})

output${prefix}_gevge_oup1.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp1,".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_gevge_oup1.res,
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp1, 
                    input${prefix}_gevge_sub1, input${prefix}_gevge_sub2,
                    "{prefix}gexpr.h5", {prefix}gene,
                    input${prefix}_gevge_siz, input${prefix}_gevge_col1, input${prefix}_gevge_ord1,
                    input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}})

output${prefix}_gevge_oup1.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp1, ".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp1, input${prefix}_gevge_sub1, input${prefix}_gevge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gevge_siz, input${prefix}_gevge_col1, input${prefix}_gevge_ord1, input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}})

output${prefix}_gevge_oup1.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp1, ".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp1, input${prefix}_gevge_sub1, input${prefix}_gevge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gevge_siz, input${prefix}_gevge_col1, input${prefix}_gevge_ord1, input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}})

output${prefix}_gevge_oup2 <- renderPlot({{
  req(input${prefix}_gevge_inp2)
  scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp2, input${prefix}_gevge_sub1, input${prefix}_gevge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gevge_siz, input${prefix}_gevge_col2, input${prefix}_gevge_ord2, input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
}})

output${prefix}_gevge_oup2.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_gevge_oup2", height = pList[input${prefix}_gevge_psz]))
}})

output${prefix}_gevge_oup2.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp2,".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_gevge_oup2.res,
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp2, input${prefix}_gevge_sub1, input${prefix}_gevge_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gevge_siz, input${prefix}_gevge_col2, input${prefix}_gevge_ord2, input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}}) 

output${prefix}_gevge_oup2.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp2,".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp2, 
                    input${prefix}_gevge_sub1, input${prefix}_gevge_sub2,
                    "{prefix}gexpr.h5", {prefix}gene,
                    input${prefix}_gevge_siz, input${prefix}_gevge_col2, input${prefix}_gevge_ord2,
                    input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}})

output${prefix}_gevge_oup2.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gevge_drX, "_", input${prefix}_gevge_drY, "_", input${prefix}_gevge_inp2,".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene({prefix}conf, {prefix}meta, input${prefix}_gevge_drX, input${prefix}_gevge_drY, input${prefix}_gevge_inp2, 
                    input${prefix}_gevge_sub1, input${prefix}_gevge_sub2,
                    "{prefix}gexpr.h5", {prefix}gene,
                    input${prefix}_gevge_siz, input${prefix}_gevge_col2, input${prefix}_gevge_ord2,
                    input${prefix}_gevge_fsz, input${prefix}_gevge_asp, input${prefix}_gevge_txt)
    )
}}) # End of tab gevge


')
}

#' Write code for gene co-expression
#'
wr_sv_gec <- function() {
paste0('
### Tab gec gene co-expression ----

{subst}  output${prefix}_gec_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gec_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_gec_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gec_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_gec_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_gec_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_gec_oup1 <- renderPlot({{
  scDRcoexFull({prefix}conf, {prefix}meta, input${prefix}_gec_drX, input${prefix}_gec_drY, input${prefix}_gec_inp1, input${prefix}_gec_inp2, input${prefix}_gec_sub1, input${prefix}_gec_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gec_siz, input${prefix}_gec_col1, input${prefix}_gec_ord1, input${prefix}_gec_fsz, input${prefix}_gec_asp, input${prefix}_gec_txt)
}})

output${prefix}_gec_oup1.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_gec_oup1", height = pList2[input${prefix}_gec_psz]))
}})

output${prefix}_gec_oup1.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gec_drX, "_", input${prefix}_gec_drY, "_", input${prefix}_gec_inp1, "_", input${prefix}_gec_inp2, ".png")) }},
  content = function(file) {{ ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_gec_oup1.res,
    plot = scDRcoexFull({prefix}conf, {prefix}meta, input${prefix}_gec_drX, input${prefix}_gec_drY, input${prefix}_gec_inp1, input${prefix}_gec_inp2, input${prefix}_gec_sub1, input${prefix}_gec_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gec_siz, input${prefix}_gec_col1, input${prefix}_gec_ord1, input${prefix}_gec_fsz, input${prefix}_gec_asp, input${prefix}_gec_txt) )
}})

output${prefix}_gec_oup1.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gec_drX, "_", input${prefix}_gec_drY, "_", input${prefix}_gec_inp1, "_", input${prefix}_gec_inp2, ".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scDRcoexFull({prefix}conf, {prefix}meta, input${prefix}_gec_drX, input${prefix}_gec_drY, input${prefix}_gec_inp1, input${prefix}_gec_inp2, input${prefix}_gec_sub1, input${prefix}_gec_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gec_siz, input${prefix}_gec_col1, input${prefix}_gec_ord1, input${prefix}_gec_fsz, input${prefix}_gec_asp, input${prefix}_gec_txt)
    )
}})

output${prefix}_gec_oup1.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gec_drX, "_", input${prefix}_gec_drY, "_", input${prefix}_gec_inp1, "_", input${prefix}_gec_inp2, ".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull({prefix}conf, {prefix}meta, input${prefix}_gec_drX, input${prefix}_gec_drY, input${prefix}_gec_inp1, input${prefix}_gec_inp2, input${prefix}_gec_sub1, input${prefix}_gec_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gec_siz, input${prefix}_gec_col1, input${prefix}_gec_ord1, input${prefix}_gec_fsz, input${prefix}_gec_asp, input${prefix}_gec_txt)
    )
}})

output${prefix}_gec_.dt <- renderDataTable({{
  ggData = scDRcoexNum({prefix}conf, {prefix}meta, input${prefix}_gec_inp1, input${prefix}_gec_inp2, input${prefix}_gec_sub1, input${prefix}_gec_sub2, "{prefix}gexpr.h5", {prefix}gene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}}) # End of tab gec

')
}

#' Write code for violinplot / boxplot
#'
wr_sv_vio <- function() {
paste0('
### Tab vio violinplot / boxplot / lineplot ----

{subst}  output${prefix}_vio_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_vio_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_vio_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_vio_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_vio_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_vio_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_vio_oup <- renderPlot({{
  gh5 <- ifelse(input${prefix}_vio_datatype == "normalised","{prefix}gexpr.h5","{prefix}gexpr2.h5")
  scVioBox({prefix}conf, {prefix}meta, input${prefix}_vio_inp1, input${prefix}_vio_inp2, input${prefix}_vio_sub1, input${prefix}_vio_sub2, gh5, {prefix}gene, input${prefix}_vio_typ, input${prefix}_vio_pts, input${prefix}_vio_siz, input${prefix}_vio_fsz, input${prefix}_vio_barsz)
}})

output${prefix}_vio_oup.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_vio_oup", height = pList2[input${prefix}_vio_psz]))
}})

output${prefix}_vio_oup.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_vio_typ, "_", input${prefix}_vio_inp1, "_", input${prefix}_vio_inp2,".png")) }},
  content = function(file) {{
    gh5 <- ifelse(input${prefix}_vio_datatype == "normalised","{prefix}gexpr.h5","{prefix}gexpr2.h5")
    ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_vio_oup.res,
    plot = scVioBox({prefix}conf, {prefix}meta, input${prefix}_vio_inp1, input${prefix}_vio_inp2, input${prefix}_vio_sub1, input${prefix}_vio_sub2, gh5, {prefix}gene, input${prefix}_vio_typ, input${prefix}_vio_pts, input${prefix}_vio_siz, input${prefix}_vio_fsz, input${prefix}_vio_barsz)
    )
}}) 

output${prefix}_vio_oup.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_vio_typ, "_", input${prefix}_vio_inp1, "_", input${prefix}_vio_inp2, ".pdf")) }},
  content = function(file) {{
    gh5 <- ifelse(input${prefix}_vio_datatype == "normalised","{prefix}gexpr.h5","{prefix}gexpr2.h5")
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scVioBox({prefix}conf, {prefix}meta, input${prefix}_vio_inp1, input${prefix}_vio_inp2, input${prefix}_vio_sub1, input${prefix}_vio_sub2, gh5, {prefix}gene, input${prefix}_vio_typ, input${prefix}_vio_pts, input${prefix}_vio_siz, input${prefix}_vio_fsz, input${prefix}_vio_barsz)
    )
}})

output${prefix}_vio_oup.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_vio_typ, "_", input${prefix}_vio_inp1, "_", input${prefix}_vio_inp2, ".svg")) }},
  content = function(file) {{
    gh5 <- ifelse(input${prefix}_vio_datatype == "normalised","{prefix}gexpr.h5","{prefix}gexpr2.h5")
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox({prefix}conf, {prefix}meta, input${prefix}_vio_inp1, input${prefix}_vio_inp2, input${prefix}_vio_sub1, input${prefix}_vio_sub2, gh5, {prefix}gene, input${prefix}_vio_typ, input${prefix}_vio_pts, input${prefix}_vio_siz, input${prefix}_vio_fsz, input${prefix}_vio_barsz)
    )
}}) # End of tab vio


')
}

#' Write code for pro
#'
wr_sv_pro <- function() {
paste0('
### Tab pro proportion plot ----

{subst}  output${prefix}_pro_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_pro_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_pro_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_pro_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_pro_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_pro_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_pro_oup <- renderPlot({{
  scProp({prefix}conf, {prefix}meta, input${prefix}_pro_inp1, input${prefix}_pro_inp2, input${prefix}_pro_sub1, input${prefix}_pro_sub2, input${prefix}_pro_typ, input${prefix}_pro_flp, input${prefix}_pro_fsz)
}})

output${prefix}_pro_oup.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_pro_oup", height = pList2[input${prefix}_pro_psz]))
}})

output${prefix}_pro_oup.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_pro_typ, "_", input${prefix}_pro_inp1, "_", input${prefix}_pro_inp2, ".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", bg = "white", dpi = input${prefix}_pro_oup.res,
    plot = scProp({prefix}conf, {prefix}meta, input${prefix}_pro_inp1, input${prefix}_pro_inp2, input${prefix}_pro_sub1, input${prefix}_pro_sub2, input${prefix}_pro_typ, input${prefix}_pro_flp, input${prefix}_pro_fsz)
    )
  }}) 
  
output${prefix}_pro_oup.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_pro_typ, "_", input${prefix}_pro_inp1, "_", input${prefix}_pro_inp2, ".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scProp({prefix}conf, {prefix}meta, input${prefix}_pro_inp1, input${prefix}_pro_inp2, input${prefix}_pro_sub1, input${prefix}_pro_sub2, input${prefix}_pro_typ, input${prefix}_pro_flp, input${prefix}_pro_fsz)
    )
  }})
  
output${prefix}_pro_oup.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_pro_typ, "_", input${prefix}_pro_inp1, "_", input${prefix}_pro_inp2, ".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp({prefix}conf, {prefix}meta, input${prefix}_pro_inp1, input${prefix}_pro_inp2, input${prefix}_pro_sub1, input${prefix}_pro_sub2, input${prefix}_pro_typ, input${prefix}_pro_flp, input${prefix}_pro_fsz)
    )
  }}) # End of tab pro

')
}

#' Write code for heatmap / dotplot
#'
wr_sv_hea <- function() {
paste0('
### Tab hea heatmap / dotplot ----

{subst}  output${prefix}_hea_sub1.ui <- renderUI({{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_hea_sub1]$fID, "\\\\|")[[1]]
{subst}    checkboxGroupInput("{prefix}_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
{subst}  }})
{subst}  observeEvent(input${prefix}_hea_sub1non, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_hea_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
{subst}  }})
{subst}  observeEvent(input${prefix}_hea_sub1all, {{
{subst}    sub = strsplit({prefix}conf[UI == input${prefix}_hea_sub1]$fID, "\\\\|")[[1]]
{subst}    updateCheckboxGroupInput(session, inputId = "{prefix}_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
{subst}  }})

output${prefix}_hea_oupTxt <- renderUI({{
  geneList = scGeneList(input${prefix}_hea_inp, {prefix}gene)
  if(nrow(geneList) > 50){{
    HTML("More than 50 input genes! Please reduce the gene list!")
  }}
}})

output${prefix}_hea_oup <- renderPlot({{
  scBubbHeat({prefix}conf, {prefix}meta, input${prefix}_hea_inp, input${prefix}_hea_grp, input${prefix}_hea_plt, input${prefix}_hea_sub1, input${prefix}_hea_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_hea_scl, input${prefix}_hea_row, input${prefix}_hea_col, input${prefix}_hea_cols, input${prefix}_hea_fsz)
}})

output${prefix}_hea_oup.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_hea_oup", height = pList3[input${prefix}_hea_psz]))
}})

output${prefix}_hea_oup.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_hea_plt,"_",input${prefix}_hea_grp,".png")) }},
  content = function(file) {{
    ggsave(
    file, device = "png", bg = "white", height = input${prefix}_hea_oup.height, width = input${prefix}_hea_oup.width, units = "cm", dpi = input${prefix}_hea_oup.res, plot = scBubbHeat({prefix}conf, {prefix}meta, input${prefix}_hea_inp, input${prefix}_hea_grp, input${prefix}_hea_plt, input${prefix}_hea_sub1, input${prefix}_hea_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_hea_scl, input${prefix}_hea_row, input${prefix}_hea_col, input${prefix}_hea_cols, input${prefix}_hea_fsz)
    )
}})

output${prefix}_hea_oup.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_hea_plt,"_",input${prefix}_hea_grp,".pdf")) }},
  content = function(file) {{
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", height = input${prefix}_hea_oup.height, width = input${prefix}_hea_oup.width, units = "cm", plot = scBubbHeat({prefix}conf, {prefix}meta, input${prefix}_hea_inp, input${prefix}_hea_grp, input${prefix}_hea_plt, input${prefix}_hea_sub1, input${prefix}_hea_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_hea_scl, input${prefix}_hea_row, input${prefix}_hea_col, input${prefix}_hea_cols, input${prefix}_hea_fsz)
    )
}})

output${prefix}_hea_oup.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_hea_plt,"_",input${prefix}_hea_grp,".svg")) }},
  content = function(file) {{
    ggsave(
    file, device = "svg", bg = "white", height = input${prefix}_hea_oup.height, width = input${prefix}_hea_oup.width, units = "cm", 
    plot = scBubbHeat({prefix}conf, {prefix}meta, input${prefix}_hea_inp, input${prefix}_hea_grp, input${prefix}_hea_plt, input${prefix}_hea_sub1, input${prefix}_hea_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_hea_scl, input${prefix}_hea_row, input${prefix}_hea_col, input${prefix}_hea_cols, input${prefix}_hea_fsz)
    )
}}) # End of tab hea      
       
')
}

#' Write code for gene expression multi
#'
wr_sv_gem <- function() {
  paste0('
### Tab gem gene expression multi ----

output${prefix}_gem_sub1.ui <- renderUI({{
  sub = strsplit({prefix}conf[UI == input${prefix}_gem_sub1]$fID, "\\\\|")[[1]]
  checkboxGroupInput("{prefix}_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
}})
observeEvent(input${prefix}_gem_sub1non, {{
  sub = strsplit({prefix}conf[UI == input${prefix}_gem_sub1]$fID, "\\\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "{prefix}_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
}})
observeEvent(input${prefix}_gem_sub1all, {{
  sub = strsplit({prefix}conf[UI == input${prefix}_gem_sub1]$fID, "\\\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "{prefix}_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
}})

output${prefix}_gem_oup1 <- renderPlot({{
  req(input${prefix}_gem_inp)
  
  scFeature({prefix}conf, {prefix}meta, input${prefix}_gem_drX, input${prefix}_gem_drY, input${prefix}_gem_inp, input${prefix}_gem_sub1, input${prefix}_gem_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gem_siz, input${prefix}_gem_col, input${prefix}_gem_ord, input${prefix}_gem_fsz, input${prefix}_gem_asp, input${prefix}_gem_txt, input${prefix}_gem_ncol)
}})

output${prefix}_gem_oup1.ui <- renderUI({{
  show_progress(imageOutput("{prefix}_gem_oup1", height = pList[input${prefix}_gem_psz]))
}})

output${prefix}_gem_oup1.png <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gem_drX, "_", input${prefix}_gem_drY, "_expression.png")) }},
  content = function(file) {{
    ggsave(
      file, device = "png", height = input${prefix}_gem_oup1.height, width = input${prefix}_gem_oup1.width, dpi = input${prefix}_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature({prefix}conf, {prefix}meta, input${prefix}_gem_drX, input${prefix}_gem_drY, input${prefix}_gem_inp, input${prefix}_gem_sub1, input${prefix}_gem_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gem_siz, input${prefix}_gem_col, input${prefix}_gem_ord, input${prefix}_gem_fsz, input${prefix}_gem_asp, input${prefix}_gem_txt, input${prefix}_gem_ncol)
    )
}}) 

output${prefix}_gem_oup1.pdf <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gem_drX, "_", input${prefix}_gem_drY, "_expression.pdf")) },
  content = function(file) {{ ggsave(
    file, device = "pdf", useDingbats = FALSE, height = input${prefix}_gem_oup1.height, width = input${prefix}_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature({prefix}conf, {prefix}meta, input${prefix}_gem_drX, input${prefix}_gem_drY, input${prefix}_gem_inp, input${prefix}_gem_sub1, input${prefix}_gem_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gem_siz, input${prefix}_gem_col, input${prefix}_gem_ord, input${prefix}_gem_fsz, input${prefix}_gem_asp, input${prefix}_gem_txt, input${prefix}_gem_ncol))
}})

output${prefix}_gem_oup1.svg <- downloadHandler(
  filename = function() {{ tolower(paste0("{prefix}", "_", input${prefix}_gem_drX, "_", input${prefix}_gem_drY, "_expression.svg")) },
  content = function(file) {{ ggsave(
    file, device = "svg", height = input${prefix}_gem_oup1.height, width = input${prefix}_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature({prefix}conf, {prefix}meta, input${prefix}_gem_drX, input${prefix}_gem_drY, input${prefix}_gem_inp, input${prefix}_gem_sub1, input${prefix}_gem_sub2, "{prefix}gexpr.h5", {prefix}gene, input${prefix}_gem_siz, input${prefix}_gem_col, input${prefix}_gem_ord, input${prefix}_gem_fsz, input${prefix}_gem_asp, input${prefix}_gem_txt, input${prefix}_gem_ncol))
}}) # End of tab gem

')
}


#' Write code for markers
#'
wr_sv_mar <- function() {
  paste0('
### Tab markers ----

output${prefix}_mar_table <- renderDataTable({{
  req(input${prefix}_mar_cls)
  datatable({prefix}mar[[input${prefix}_mar_cls]], rownames = FALSE, extensions = "Buttons", options = list(dom = "lftiprB", buttons = c("copy", "csv", "excel")))
}}) # End of tab mar

')
}

#' Write code for main block of server.R
#' @param prefix file prefix
#' @param subst Conditional
#' @param tabs Vector of tab names to include
#' @rdname wr_sv_main
#' @export wr_sv_main
#'
wr_sv_main <- function(prefix, subst = "", tabs = c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea")) {
  
glue::glue(
'optCrt="{{ option_create: function(data,escape) {{return(\'<div class=\\"create\\"><strong>\' + \'</strong></div>\');}} }}"
updateSelectizeInput(session, "{prefix}_civge_inp2", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_gevge_inp1", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_gevge_inp2", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_gec_inp1", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_gec_inp2", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_gem_inp", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_hea_inp", choices = names({prefix}gene), server = TRUE,
                     selected = {prefix}def$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "{prefix}_vio_inp2", server = TRUE,
                     choices = c({prefix}conf[is.na(fID)]$UI,names({prefix}gene)),
                     selected = {prefix}conf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length({prefix}conf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))',
paste(unlist(lapply(tabs, function(x) eval(parse(text = paste0("wr_sv_",x,"()"))))), collapse = "\n")
)
}

#' Write code for final portion of server.R
#'
#' @rdname wr_sv_end
#' @export wr_sv_end
#'
wr_sv_end <- function() {
  paste0(
'
})


'
  )
}

#' Write code for front portion of ui.R
#' @param title shiny app title
#' @param theme bootstrap theme
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#' @importFrom shinythemes shinytheme
#' @rdname wr_ui_single
#' @export wr_ui_single
#'
wr_ui_single <- function(title, theme = "flatly", ganalytics = NA) {
  if (!is.na(ganalytics)) {
    ga <- 'tags$head(includeHTML(("google-analytics.html"))),'
  } else {
    ga <- ""
  }
  
  glue::glue(
'

### UI code
shinyUI(
fluidPage(style="margin:0;padding:0;",
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
{ga}
theme = shinythemes::shinytheme("{theme}"),
navbarPage(
"{title}"
'
  )
}

#' Write code for civge
#'
wr_ui_civge <- function() {
paste0(',
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("{prefix}_civge_drX", "X-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[1]
            ),
            selectInput("{prefix}_civge_drY", "Y-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        {subst}column(
          {subst}4,
          {subst}div(
            {subst}class = "input-panel",
            {subst}checkboxInput("{prefix}_civge_togL", "Subset cells"),
            {subst}conditionalPanel(
              {subst}condition = "input.{prefix}_civge_togL == true",
              {subst}selectInput("{prefix}_civge_sub1", "Cell info to subset:",
                {subst}choices = {prefix}conf[grp == TRUE]$UI,
                {subst}selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_civge_sub1.ui"),
              {subst}actionButton("{prefix}_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            {subst})
          {subst})
        {subst}), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("{prefix}_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_civge_tog0 == true",
              sliderInput("{prefix}_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("{prefix}_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("{prefix}_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("{prefix}_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_civge_inp1", "Cell info:",
                  choices = {prefix}conf$UI,
                  selected = {prefix}def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_civge_tog1 == true",
                  radioButtons("{prefix}_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("{prefix}_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("{prefix}_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("{prefix}_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("{prefix}_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("{prefix}_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("{prefix}_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("{prefix}_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("{prefix}_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.{prefix}_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("{prefix}_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("{prefix}_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_civge_tog2 == true",
                  radioButtons("{prefix}_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("{prefix}_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("{prefix}_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("{prefix}_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("{prefix}_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("{prefix}_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("{prefix}_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge
')
}

#' Write code for civci
#'
#'
wr_ui_civci <- function() {
paste0(',
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("{prefix}_civci_drX", "X-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[1]
            ),
            selectInput("{prefix}_civci_drY", "Y-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        {subst}column(4,
          {subst}div(
            {subst}class = "input-panel",
            {subst}checkboxInput("{prefix}_civci_togL", "Subset cells"),
            {subst}conditionalPanel(
              {subst}condition = "input.{prefix}_civci_togL == true",
              {subst}selectInput("{prefix}_civci_sub1", "Cell info to subset:",
               {subst} choices = {prefix}conf[grp == TRUE]$UI,
               {subst} selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_civci_sub1.ui"),
              {subst}actionButton("{prefix}_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            {subst})
          {subst})
        {subst}), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("{prefix}_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_civci_tog0 == true",
              sliderInput("{prefix}_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("{prefix}_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("{prefix}_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("{prefix}_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_civci_inp1", "Cell info:",
                  choices = {prefix}conf$UI,
                  selected = {prefix}def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_civci_tog1 == true",
                  radioButtons("{prefix}_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("{prefix}_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("{prefix}_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("{prefix}_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("{prefix}_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("{prefix}_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("{prefix}_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("{prefix}_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_civci_inp2", "Cell info:",
                  choices = {prefix}conf$UI,
                  selected = {prefix}def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_civci_tog2 == true",
                  radioButtons("{prefix}_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("{prefix}_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("{prefix}_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("{prefix}_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("{prefix}_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("{prefix}_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("{prefix}_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("{prefix}_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci
')
}

#' Write code for gevge
#'
wr_ui_gevge <- function() {
paste0(',
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("{prefix}_gevge_drX", "X-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[1]
            ),
            selectInput("{prefix}_gevge_drY", "Y-axis:",
              choices = {prefix}conf[dimred == TRUE]$UI,
              selected = {prefix}def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        {subst}column(4,
          {subst}div(
            {subst}class = "input-panel",
            {subst}checkboxInput("{prefix}_gevge_togL", "Subset cells"),
            {subst}conditionalPanel(
             {subst} condition = "input.{prefix}_gevge_togL == true",
              {subst}selectInput("{prefix}_gevge_sub1", "Cell info to subset:",
               {subst} choices = {prefix}conf[grp == TRUE]$UI,
                {subst}selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_gevge_sub1.ui"),
              {subst}actionButton("{prefix}_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            {subst})
          {subst})
        {subst}), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("{prefix}_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_gevge_tog0 == true",
              sliderInput("{prefix}_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("{prefix}_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("{prefix}_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("{prefix}_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_gevge_tog1 == true",
                  radioButtons("{prefix}_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("{prefix}_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("{prefix}_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("{prefix}_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("{prefix}_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("{prefix}_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("{prefix}_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("{prefix}_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("{prefix}_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.{prefix}_gevge_tog2 == true",
                  radioButtons("{prefix}_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("{prefix}_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("{prefix}_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("{prefix}_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("{prefix}_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("{prefix}_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("{prefix}_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge
')
}

#' Write code for gec
#'
wr_ui_gec <- function() {
paste0(',
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("{prefix}_gec_drX", "X-axis:",
                               choices = {prefix}conf[dimred == TRUE]$UI,
                               selected = {prefix}def$dimred[1]
                   ),
                   selectInput("{prefix}_gec_drY", "Y-axis:",
                               choices = {prefix}conf[dimred == TRUE]$UI,
                               selected = {prefix}def$dimred[2]
                   ),
                   selectInput("{prefix}_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("{prefix}_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   {subst}checkboxInput("{prefix}_gec_togL", "Subset cells"),
                   {subst}conditionalPanel(
                    {subst}condition = "input.{prefix}_gec_togL == true",
                    {subst}selectInput("{prefix}_gec_sub1", "Cell info to subset:", choices = {prefix}conf[grp == TRUE]$UI, selected = {prefix}def$grp1),
                     {subst}uiOutput("{prefix}_gec_sub1.ui"),
                     {subst}actionButton("{prefix}_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     {subst}actionButton("{prefix}_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   {subst}),
                   checkboxInput("{prefix}_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.{prefix}_gec_tog0 == true",
                     radioButtons("{prefix}_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("{prefix}_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("{prefix}_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("{prefix}_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("{prefix}_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("{prefix}_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("{prefix}_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("{prefix}_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("{prefix}_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("{prefix}_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("{prefix}_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("{prefix}_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("{prefix}_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec
')
}

#' Write code for vio
#'
wr_ui_vio<- function() {
paste0(',
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("{prefix}_vio_inp1", "Cell info (X-axis):",
              choices = {prefix}conf[grp == TRUE]$UI,
              selected = {prefix}def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("{prefix}_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("{prefix}_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("{prefix}_vio_pts", "Show data points", value = FALSE),
            {subst}checkboxInput("{prefix}_vio_togL", "Subset cells"),
            {subst}conditionalPanel(
              {subst}condition = "input.{prefix}_vio_togL == true",
              {subst}selectInput("{prefix}_vio_sub1", "Cell info to subset:",
               {subst} choices = {prefix}conf[grp == TRUE]$UI,
               {subst} selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_vio_sub1.ui"),
              {subst}actionButton("{prefix}_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            {subst}),
            checkboxInput("{prefix}_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_vio_tog == true",
              sliderInput("{prefix}_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("{prefix}_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("{prefix}_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.{prefix}_vio_typ == \'lineplot\'",
              sliderInput("{prefix}_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("{prefix}_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("{prefix}_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("{prefix}_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("{prefix}_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("{prefix}_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("{prefix}_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio
')
}

#' Write code for pro
#'
wr_ui_pro <- function() {
paste0(',
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("{prefix}_pro_inp1", "Cell info to plot (X-axis):",
              choices = {prefix}conf[grp == TRUE]$UI,
              selected = {prefix}def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("{prefix}_pro_inp2", "Cell info to group / colour by:",
              choices = {prefix}conf[grp == TRUE]$UI,
              selected = {prefix}def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("{prefix}_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("{prefix}_pro_flp", "Flip X/Y", value = FALSE),
            {subst}checkboxInput("{prefix}_pro_togL", "Subset cells"),
            {subst}conditionalPanel(
             {subst} condition = "input.{prefix}_pro_togL == true",
             {subst} selectInput("{prefix}_pro_sub1", "Cell info to subset:",
              {subst}  choices = {prefix}conf[grp == TRUE]$UI,
               {subst} selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_pro_sub1.ui"),
              {subst}actionButton("{prefix}_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("{prefix}_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_pro_tog == true",
              radioButtons("{prefix}_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("{prefix}_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("{prefix}_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("{prefix}_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("{prefix}_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("{prefix}_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro
')
}

#' Write code for hea
#'
wr_ui_hea <- function() {
paste0(',
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("{prefix}_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("{prefix}_hea_grp", "Group by:",
              choices = {prefix}conf[grp == TRUE]$UI,
              selected = {prefix}conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("{prefix}_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("{prefix}_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("{prefix}_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("{prefix}_hea_col", "Cluster columns (samples)", value = FALSE),
            {subst}checkboxInput("{prefix}_hea_togL", "Subset cells"),
            {subst}conditionalPanel(
              {subst}condition = "input.{prefix}_hea_togL == true",
              {subst}selectInput("{prefix}_hea_sub1", "Cell info to subset:",
                {subst}choices = {prefix}conf[grp == TRUE]$UI,
                {subst}selected = {prefix}def$grp1
              {subst}),
              {subst}uiOutput("{prefix}_hea_sub1.ui"),
              {subst}actionButton("{prefix}_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              {subst}actionButton("{prefix}_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            {subst}),
            checkboxInput("{prefix}_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_hea_tog == true",
              radioButtons("{prefix}_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("{prefix}_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("{prefix}_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("{prefix}_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("{prefix}_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("{prefix}_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("{prefix}_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("{prefix}_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("{prefix}_hea_oupTxt")),
          uiOutput("{prefix}_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea
')
}

#' Write code for gem
#'
wr_ui_gem <- function() {
paste0(',
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("{prefix}_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("{prefix}_gem_drX", "X-axis:",
                        choices = {prefix}conf[dimred == TRUE]$UI,
                        selected = {prefix}def$dimred[1]
            ),
            selectInput("{prefix}_gem_drY", "Y-axis:",
                        choices = {prefix}conf[dimred == TRUE]$UI,
                        selected = {prefix}def$dimred[2]
            ),
            checkboxInput("{prefix}_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.{prefix}_gem_togL == true",
              selectInput("{prefix}_gem_sub1", "Cell info to subset:",
                          choices = {prefix}conf[grp == TRUE]$UI,
                          selected = {prefix}def$grp1
              ),
              uiOutput("{prefix}_gem_sub1.ui"),
              actionButton("{prefix}_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("{prefix}_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("{prefix}_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.{prefix}_gem_tog0 == true",
              sliderInput("{prefix}_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("{prefix}_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("{prefix}_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("{prefix}_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("{prefix}_gem_txt", "Show axis text", value = FALSE),
              radioButtons("{prefix}_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("{prefix}_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("{prefix}_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("{prefix}_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("{prefix}_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("{prefix}_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("{prefix}_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("{prefix}_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("{prefix}_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("{prefix}_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem
')
}

#' Write code for mar
#'
wr_ui_mar <- function() {
  paste0(',
# tab mar ----
tabPanel(
  "Markers",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Markers"),
          p("Explore markers for different clustering.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(3,
               fluidRow(
                 column(
                   12,
                   div(
                     class = "input-panel input-panel-section",
                     selectInput("{prefix}_mar_cls","Select clustering:", choices = names({prefix}mar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("{prefix}_mar_table")
      )
      ), # end of row 3
      hr()
    )
  )
) # End of tab mar
')
}


#' Write code for about page
#'
wr_ui_about <- function() {
paste0('\n,
# about ----
tabPanel(
  "About",
  fluidRow(
    class = "container page",
    column(
      12,
      includeMarkdown("about.md")
    )
  )
)'
)
}

#' Write code for main block of ui.R
#' @param prefix file prefix
#' @param subst Conditional
#' @param ptsiz Point size
#' @param tabs Vector of tab names to include
#' @param about Should about page be added as a tab?
#' @rdname wr_ui_main
#' @export wr_ui_main
#'
wr_ui_main <- function(prefix, subst = "", ptsiz = "1.25", tabs = c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea"), about = TRUE) {
  glue::glue(
    paste(unlist(lapply(tabs, function(x) eval(parse(text = paste0("wr_ui_",x,"()"))))), collapse = "\n"),
    "\n"
  )
}

#' Write code for final portion of ui.R
#'
#' @rdname wr_ui_end
#' @export wr_ui_end
#'
wr_ui_end <- function() {
  glue::glue(
'
)
)
)\n
'
  )
}

#' Write code for google-analytics.html
#' @param gaID Google analytics tracking ID (e.g. "UA-123456789-0")
#'
#' @rdname wr_ui_ga
#' @export wr_ui_ga
#'
wr_ui_ga <- function(gaID) {
  glue::glue(
    '<!-- Global site tag (gtag.js) - Google Analytics -->
      <script async src="https://www.googletagmanager.com/gtag/js?id={gaID}"></script>
      <script>\n",
      window.dataLayer = window.dataLayer || [];
    function gtag(){{dataLayer.push(arguments);}}
    gtag("js", new Date());
    
    gtag("config", "{gaID}");
    </script>\n'
  )
}

#' Write code for about. about.md
#' @rdname wr_about
#' @export wr_about
#'
wr_about <- function(){
  paste('
### About

App description and author info.

')  
}

#' Write code for custom css. www/styles.css
#' @rdname wr_css
#' @export wr_css
#'
wr_css <- function(){
  paste('
/* easyshiny */
/* custom css styles */

.page {
    margin: auto;
}

.tab-section {
  margin-bottom: 1.5em;
}

.input-panel {
  background-color: #f3f6f4;
  padding: 15px;
  border-radius: 6px;
}

.input-panel-section {
  margin-bottom: 1em;
}

a:hover,
a:focus,
a:active,
.navbar-default .navbar-nav>li>a:hover,
.navbar-default .navbar-nav>li>a:focus,
.navbar-default .navbar-nav>li>a:active,
.navbar-default .navbar-brand:hover,
.navbar-default .navbar-brand:focus,
.navbar-default .navbar-brand:active {
  color: #bcbcbc;
}

.btn {
  margin: 5px 0px;
}

.dt-buttons {
  padding-top: 0.25em;
  
}

.dataTables_info,
.dataTables_paginate,
.dt-buttons {
  margin: 0.5em 0em;
}

')  
}
