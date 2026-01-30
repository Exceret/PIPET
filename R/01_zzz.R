# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
}

#' @importFrom data.table `:=` `%chin%`
NULL

#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()
