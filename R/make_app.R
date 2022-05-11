#' Make a shiny app
#'
#' Make a shiny app based on the shinycell config data.table and single-cell
#' data object.
#'
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment
#'   data or input file path for h5ad / loom files
#' @param scConf shinycell config data.table
#' @param assay1 assay in single-cell data object to use for plotting
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay,
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts",
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'     \item{loom files}: "matrix" or any assay in "layers",
#'       default is "matrix"
#'   }
#' @param slot1 slot in single-cell assay to plot. This is only used
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param assay2 secondary assay to use (Optional)
#' @param slot2 secondary slot to use (Optional). This is only used
#'   for Seurat objects (v3+).
#' @param gene.mapping specifies whether to convert human / mouse Ensembl gene
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols.
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users
#'   can supply a named vector where \code{names(gene.mapping)} correspond
#'   to the actual gene identifiers in the gene expression matrix and
#'   \code{gene.mapping} correspond to new identifiers to map to
#' @param shiny.title title for shiny app
#' @param shiny.dir specify directory to create the shiny app in. Default is
#'   to create a new directory named "shinyApp"
#' @param shiny.prefix specify file prefix 
#' @param enableSubset specify whether to enable "Toggle to subset cells"
#'   functionality in the shiny app. Default is to enable this functionality
#' @param defPtSiz specify default point size for single cells. For example, a
#'   smaller size can be used if you have many cells in your dataset
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying the default genes to
#'   show in bubbleplot / heatmap
#' @param default.dimred character vector specifying the two default dimension
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param tabs Vector of tab names to include
#' @param theme Bootstrap theme
#' @param font Google font for plots
#'
#' @return directory containing shiny app
#'
#' @author John F. Ouyang
#'
#' @export
#' 
make_app <- function(obj, scConf, assay1 = NA, slot1 = c("data", "scale.data", "counts"), assay2 = NA, slot2 = NA, gene.mapping = FALSE, shiny.title = "scRNA-seq shiny app", shiny.dir = "shinyApp/", shiny.prefix = "sc1", enableSubset = TRUE, defPtSiz = 1.25, ganalytics = NA, default.gene1 = NA, default.gene2 = NA, default.multigene = NA, default.dimred = NA, tabs = c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea"), theme = "flatly", font = "Lato") {

  tbs <- c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea", "mar")
  if(length(tabs) < 1) stop("At least 1 tab must be specified.")
  if(any(!tabs %in% tbs)) stop(paste("One of more tabs are incorrect. Tab options are:",paste(tbs,collapse=", "),"."))
  
  # Checks are performed in respective functions
  # Wrapper for two main functions
  make_file(
    obj = obj, scConf = scConf,
    assay1 = assay1[1], slot1 = slot1[1],
    assay2 = assay2[1], slot2 = slot2[1],
    gene.mapping = gene.mapping,
    shiny.prefix = shiny.prefix, shiny.dir = shiny.dir,
    default.gene1, default.gene2, default.multigene, default.dimred,
    tabs = tabs
  )

  make_code(
    shiny.title = shiny.title,
    shiny.prefix = shiny.prefix, shiny.dir = shiny.dir,
    theme = theme, enableSubset = enableSubset, defPtSiz = defPtSiz,
    ganalytics = ganalytics, tabs = tabs, font = font
  )
}
