#' Load specific version of MMVmalaria package
#'
#' Installs the requested version from the **local** MMVmalaria
#' git repository and loads the package.
#' Please make sure that your local repository
#' is synchronized with the github repository to be sure that
#' all version numbers are locally available.
#'
#' The function requires that the user has created MMVmalaria
#' options via
#'
#' ```
#' options(MMVmalaria = list(
#'   user.name = "Daniel Kaschek, IntiQuan",
#'   user.email = "daniel.kaschek@intiquan.com",
#'   local.repo = "C:/PROJTOOLS/MMVmalaria"
#' ))
#' ```
#'
#' The general procedure to create new versions of the MMVmalaria package is:
#' 1. Modify the DESCRIPTION file and set the new version number
#' 2. Commit all changes
#'     ```
#'     git add .
#'     git commit -m "commit message"
#'     ```
#' 3. Add tag to your commit with the new version number, e.g.
#'     ```
#'     git tag v0.1.2
#'     ```
#' 4. Synchronize with github repository
#'     ```
#'     git pull
#'     git push
#'     git push origin --tags
#'     ```
#'
#' @param version version number as charager, e.g., "v0.1.2".
#' @export
#' @md
#' @family Package Administration
loadversion_MMVmalaria <- function(version) {

  myop <- getOption("MMVmalaria")
  repo <- git2r::repository(myop$local.repo)
  git2r::config(repo, user.name = myop$user.name, user.email = myop$user.email)

  devtools::install_git(myop$local.repo, subdir = "MMVmalaria", branch = version)
}


#' Build a dependency tree of function sources and targets from a package
#'
#' @param package package name
#' @return Data frame of all functions found in the package. The row "n(from)"
#' indicates the number of other functions where the actual function is used. The
#' row "n(to)" indicates the number of other functions used within the actual function.
#'
#' The output carries several attributes which are needed in \code{\link{get_targets}} and
#' \code{\link{get_sources}}.
#'
#' @export
#'
#' @examples
#' # Get all functions from MMVmalaria package
#' functions <- listfunctions_MMVmalaria()
#'
#' # Copy a list of all functions to the clipboard to paste into Excel (or Google Docs)
#' write.table(
#'   data.frame(Function = names(functions)),
#'   file = "clipboard",
#'   sep = "\t",
#'   row.names = FALSE
#' )
listfunctions_MMVmalaria <- function(package = "MMVmalaria", onlyExported = FALSE) {

  mynamespace <- getNamespace(package)
  allfns <- as.character(lsf.str(mynamespace))
  if (onlyExported) allfns <- intersect(allfns, getNamespaceExports(package))

  nomfun <- data.frame(id = seq_along(allfns), label = allfns)
  fromto <- do.call(rbind, lapply(seq_along(allfns), function(from) {
    from.fn <- allfns[from]
    from.fn.body <- as.character(body(getFromNamespace(from.fn, mynamespace)))
    to <- which(sapply(allfns, function(x) {
      direct <- any(grepl(paste0(x, "("), from.fn.body, fixed = TRUE) & !grepl(paste0("\"", x, "("), from.fn.body, fixed = TRUE) & !grepl(paste0("_", x, "("), from.fn.body, fixed = TRUE))
      incall <- any(grepl(paste0("do.call(", x), from.fn.body, fixed = TRUE))
      inmethod <- strsplit(x, split = ".", fixed = TRUE)[[1]][1] == from.fn & x != from.fn
      direct | incall | inmethod
    }))
    if (length(to) == 0) return()
    data.frame(from = from, to = to)
  }))


  n_from <- table(fromto[["from"]])
  n_to <- table(fromto[["to"]])

  out <- matrix(0, nrow = 2, ncol = length(nomfun[["id"]]))
  out[1, ] <- n_from[as.character(nomfun[["id"]])]
  out[2, ] <- n_to[as.character(nomfun[["id"]])]

  out <- as.data.frame(out)
  colnames(out) <- nomfun[["label"]]
  rownames(out) <- c("n(from)", "n(to)")

  attr(out, "allfns") <- allfns
  attr(out, "fromto") <- fromto
  attr(out, "nomfun") <- nomfun

  # Output:
  return(out)
}

