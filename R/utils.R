# ==========================================================================
# SVMC: small internal utilities
# Null-coalescing operator, character coalesce, URL existence check.
# ==========================================================================

`%||%` <- function(a,b) ifelse(is.null(a), b, a)


coalesce_chr <- function(...) {
  x <- list(...)
  out <- x[[1]]
  for (i in seq_along(x)) out <- dplyr::coalesce(out, x[[i]])
  out
}


file_exists_at_url <- 
  function(url) {
    con <- tryCatch({
      suppressWarnings({
        con <- url(url, "r")  # Try to open the URL
        close(con)            # Close the connection if it was successful
        TRUE                  # Return TRUE if the file exists
      })
    }, error = function(e) {
      FALSE  # Return FALSE if there's an error (file doesn't exist)
    })
    return(con)
  }
