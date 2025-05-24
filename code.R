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

ematrix <- function(seqs, states, symbols) {
  ecounts <- matrix(1, nrow = length(states), ncol = length(symbols))  # Laplace smoothing
  rownames(ecounts) <- states
  colnames(ecounts) <- symbols
  for (seq in seqs) {
    for (i in seq_along(seq)) {
      if (i <= length(states)) {
        aa <- seq[i]
        if (aa %in% symbols) {
          ecounts[states[i], aa] <- ecounts[states[i], aa] + 1
        }
      }
    }
  }
  emission_probs <- ecounts / rowSums(ecounts)
  return(emission_probs)
}

prob_path <- function(hmm, obs_seq, viterbi_path) {
  log_prob <- log(hmm$startProbs[viterbi_path[1]]) +
              log(hmm$emissionProbs[viterbi_path[1], obs_seq[1]])
  for (t in 2:length(obs_seq)) {
    prev_state <- viterbi_path[t - 1]
    curr_state <- viterbi_path[t]
    symbol <- obs_seq[t]
    log_prob <- log_prob + log(hmm$transProbs[prev_state, curr_state]) +
                             log(hmm$emissionProbs[curr_state, symbol])
  }
  return(log_prob)
}
