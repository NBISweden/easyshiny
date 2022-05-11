#' Generate code files required for shiny app (multi datasets)
#'
#' Generate code files required for shiny app containing multiple datasets. In
#' particular, two R scripts will be generated, namely \code{server.R} and
#' \code{ui.R}. Note that \code{make_file} has to be ran prior to
#' generate the necessary data files for each of the dataset to be included.
#' The prefix used in \code{make_file} have to be then supplied in this
#' function.
#'
#' @param shiny.title specify the overall title for shiny app
#' @param shiny.prefix specify file prefix for each dataset. Must match the
#'   prefix used in \code{make_file}
#' @param shiny.headers specify the tab header names for each dataset. Length
#'   must match that of \code{shiny.prefix}
#' @param shiny.dir specify directory to create the shiny app in
#' @param enableSubset specify whether to enable "Toggle to subset cells"
#'   functionality in the shiny app. Default is to enable this functionality
#' @param defPtSiz specify default point size for single cells. For example, a
#'   smaller size can be used if you have many cells in your dataset. A single
#'   value can be specified to set the point size for all datasets. Otherwise,
#'   users have to specify one value for each dataset
#' @param theme Bootsrap theme
#' @param tabs Vector of tab names to include
#' @param about Should about page be added as a tab?
#' @param font Google font for plots
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
make_code_multi <- function(shiny.title, shiny.prefix, shiny.headers, shiny.dir, enableSubset = TRUE, defPtSiz = 1.25, theme = "flatly", tabs = c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea"), about = TRUE, font = "Lato", ganalytics = NA) {
  
  ### Checks
  tbs <- c("civge", "civci", "gevge", "gem", "gec", "vio", "pro", "hea", "mar")
  if(length(tabs) < 1) stop("At least 1 tab must be specified.")
  if(any(!tabs %in% tbs)) stop(paste("One of more tabs are incorrect. Tab options are:",paste(tbs,collapse=", "),"."))
  
  if (length(shiny.prefix) != length(shiny.headers)) {
    stop("length of shiny.prefix and shiny.headers does not match!")
  }
  subst <- "#"
  if (enableSubset) {
    subst <- ""
  }
  if (length(shiny.prefix) != length(defPtSiz)) {
    defPtSiz <- rep(defPtSiz[1], length(shiny.prefix))
  }
  defPtSiz <- as.character(defPtSiz)
  slibs <- c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr", "ggplot2", "ggplotify", "ggrepel", "hdf5r", "ggdendro", "grid", "shinycssloaders","patchwork")
  ulibs <- c("shiny", "shinyhelper", "shinythemes", "showtext", "data.table", "Matrix", "DT", "magrittr")


  ### Write code for server.R
  fname <- paste0(shiny.dir, "/server.R")
  readr::write_file(wr_lib(slibs), file = fname)
  readr::write_file(wr_font(font = font), append = TRUE, file = fname)
  for (i in shiny.prefix) {
    readr::write_file(wr_load(i, tabs = tabs), append = TRUE, file = fname)
  }
  readr::write_file(wr_sv_fix(font = font), append = TRUE, file = fname)
  for (i in shiny.prefix) {
    readr::write_file(wr_sv_main(i, subst, font = font, tabs = tabs), append = TRUE, file = fname)
  }
  readr::write_file(wr_sv_end(), append = TRUE, file = fname)


  ### Write code for ui.R
  fname <- paste0(shiny.dir, "/ui.R")
  readr::write_file(wr_lib(ulibs), file = fname)
  for (i in shiny.prefix) {
    readr::write_file(wr_load(i, tabs = tabs), append = TRUE, file = fname)
  }
  readr::write_file(wr_ui_single(shiny.title, theme = theme, ganalytics = ganalytics), append = TRUE, file = fname)
  for (i in seq_along(shiny.prefix)) {
    hhh <- shiny.headers[i]
    readr::write_file(glue::glue('\n\n,navbarMenu("{hhh}"'), append = TRUE, file = fname)
    readr::write_file(wr_ui_main(shiny.prefix[i], subst, defPtSiz[i], tabs = tabs, about = FALSE), append = TRUE, file = fname)
    readr::write_file("\n)", append = TRUE, file = fname)
  }
  if (about) {
    readr::write_file(glue::glue('\n\n,navbarMenu("About"'), append = TRUE, file = fname)
    readr::write_file(wr_ui_about(), append = TRUE, file = fname)
    readr::write_file("\n)", append = TRUE, file = fname)
  }
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
