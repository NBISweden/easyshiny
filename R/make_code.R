#' Generate code files required for shiny app (one dataset)
#'
#' Generate code files required for shiny app containing only one dataset. In
#' particular, two R scripts will be generated, namely \code{server.R} and
#' \code{ui.R}. If users want to include multiple dataset in one shiny app,
#' please use \code{make_code_multi()} instead. Note that both
#' \code{make_file} and \code{make_code} functions are ran when
#' running the wrapper function \code{make_app}.
#'
#' @param shiny.title title for shiny app
#' @param shiny.prefix specify file prefix
#' @param shiny.dir specify directory to create the shiny app in
#' @param enableSubset specify whether to enable "Toggle to subset cells"
#'   functionality in the shiny app. Default is to enable this functionality
#' @param defPtSiz specify default point size for single cells. For example, a
#'   smaller size can be used if you have many cells in your dataset
#' @param theme Bootsrap theme
#' @param tabs Vector of tab names to include
#' @param about Should about page be added as a tab?
#' @param font Google font for plots. Defaults to "Lato".
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr glue shiny
#' @importFrom utils packageVersion
#' @importFrom ggplotify as.ggplot
#'
#' @export
#' 
make_code <- function(shiny.title = "App", shiny.prefix = "sc1", shiny.dir = "app", enableSubset = TRUE, defPtSiz = 1.25, theme = "flatly", tabs = c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea"), about = TRUE, font = "Lato", ganalytics = NA) {
  
  tbs <- c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea", "mar")
  if(length(tabs) < 1) stop("At least 1 tab must be specified.")
  if(any(!tabs %in% tbs)) stop(paste("One of more tabs are incorrect. Tab options are:",paste(tbs,collapse=", "),"."))
  
  
  subst <- "#"
  if (enableSubset) {
    subst <- ""
  }
  defPtSiz <- as.character(defPtSiz)
  slibs <- c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr", "ggplot2", "ggplotify", "ggrepel", "hdf5r", "ggdendro", "grid", "shinycssloaders", "patchwork")
  ulibs <- c("shiny", "shinyhelper", "shinythemes", "showtext", "data.table", "Matrix", "DT", "magrittr")

  ### Write code for server.R
  fname <- paste0(shiny.dir, "/server.R")
  readr::write_file(wr_lib(slibs), file = fname)
  readr::write_file(wr_font(font = font), append = TRUE, file = fname)
  readr::write_file(wr_load(shiny.prefix, tabs = tabs), append = TRUE, file = fname)
  readr::write_file(wr_sv_fix(font = font), append = TRUE, file = fname)
  readr::write_file(wr_sv_main(shiny.prefix, subst, tabs = tabs), append = TRUE, file = fname)
  readr::write_file(wr_sv_end(), append = TRUE, file = fname)


  ### Write code for ui.R
  fname <- paste0(shiny.dir, "/ui.R")
  readr::write_file(wr_lib(ulibs), file = fname)
  readr::write_file(wr_load(shiny.prefix, tabs = tabs), append = TRUE, file = fname)
  readr::write_file(wr_ui_single(shiny.title, theme = theme, ganalytics = ganalytics), append = TRUE, file = fname)
  readr::write_file(wr_ui_main(shiny.prefix, subst, defPtSiz, tabs = tabs, about = about), append = TRUE, file = fname)
  ## readr::write_file('\n)\n', append = TRUE, file = fname)
  readr::write_file(wr_ui_end(), append = TRUE, file = fname)

  ### Write code for google-analytics.html
  if (!is.na(ganalytics)) {
    fname <- paste0(shiny.dir, "/google-analytics.html")
    readr::write_file(wr_ui_ga(ganalytics), file = fname)
  }

  ### Write extra css
  if (!dir.exists(file.path(shiny.dir, "www"))) dir.create(path = file.path(shiny.dir, "www"))
  path_css <- file.path(shiny.dir, "www", "styles.css")
  if (!file.exists(path_css)) file.create(path_css)
  readr::write_file(wr_css(), file = path_css)

  ### Write about
  if (about) {
    path_about <- file.path(shiny.dir, "about.md")
    if (!file.exists(path_about)) file.create(path_about)
    readr::write_file(wr_about(), file = path_about)
  }
}
