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

args <- commandArgs(trailingOnly=TRUE)
train_file <- args[1]
test_file <- args[2]
train_seqs <- read_1(train_file)
normal_seqs <- train_seqs$normal
mutation_seqs <- train_seqs$mutation
max_len <- max(sapply(c(normal_seqs, mutation_seqs), length))
num_match_states <- min(50, max_len)
num_insert_states <- num_match_states
states <- c(paste0("M", 1:num_match_states), paste0("I", 1:num_insert_states))
symbols <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
emission_normal <- ematrix(normal_seqs, states, symbols)
emission_mutation <- ematrix(mutation_seqs, states, symbols)
trans_mat <- matrix(0, nrow = length(states), ncol = length(states))
rownames(trans_mat) <- states
colnames(trans_mat) <- states
start_probs <- rep(0, length(states))
start_probs[1] <- 1.0
hmm_normal <- initHMM(states, symbols, start_probs, trans_mat, emission_normal)
hmm_mutation <- initHMM(states, symbols, start_probs, trans_mat, emission_mutation)
test_seq <- readLines(test_file)
test_seq <- test_seq[nzchar(test_seq)]
test_seq <- strsplit(test_seq[1], "")[[1]]
test_seq <- test_seq[test_seq %in% symbols]
viterbi_normal <- viterbi(hmm_normal, test_seq)
viterbi_mutation <- viterbi(hmm_mutation, test_seq)
logprob_normal <- prob_path(hmm_normal, test_seq, viterbi_normal)
logprob_mutation <- prob_path(hmm_mutation, test_seq, viterbi_mutation)
if (logprob_mutation < logprob_normal) {
  cat("The sequence is mutated\n")
} else {
  cat("The sequence is normal\n")
}
