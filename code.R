#!/usr/bin/env Rscript
library(HMM)

read_1 <- function(filepath) {
  lines <- readLines(filepath)
  lines <- lines[nzchar(lines)]
  seqs <- list(normal = list(), mutation = list())
  current_label <- NULL
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (grepl("normal", tolower(line))) {
        current_label <- "normal"
      } else if (grepl("mutated|mutation|mutant", tolower(line))) {
        current_label <- "mutation"
      }
    } else if (!is.null(current_label)) {
      seqs[[current_label]] <- append(seqs[[current_label]], list(strsplit(line, "")[[1]]))
    }
  }
  return(seqs)
}
