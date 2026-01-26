#' Locate MAGMA binary
#'
#' Finds the path to the MAGMA executable. MAGMA can be located by:
#' \enumerate{
#'   \item Setting the R option \code{options(magma.path = "/path/to/magma")}
#'   \item Having 'magma' in your system PATH
#' }
#'
#' @return Full path to the magma executable
#'
#' @examples
#' \dontrun{
#' # Get MAGMA path
#' mp <- magma_path()
#' print(mp)
#'
#' # Set custom MAGMA path
#' options(magma.path = "/usr/local/bin/magma")
#' }
#'
#' @seealso \code{\link{magma_gene}}
#' @export
magma_path <- function() {
  # 1) user can set an option in .Rprofile:
  # options(magma.path = "/path/to/magma")
  opt <- getOption("magma.path")
  if (!is.null(opt) && file.exists(opt)) {
    return(opt)
  }

  # 2) otherwise look in PATH
  path <- Sys.which("magma")
  if (nzchar(path)) {
    return(path)
  }

  stop(
    "MAGMA binary not found. ",
    "Install MAGMA for your local device) and either:\n",
    "  * put 'magma' in your PATH, or\n",
    "  * set options(magma.path = '/full/path/to/magma')",
    call. = FALSE
  )
}
