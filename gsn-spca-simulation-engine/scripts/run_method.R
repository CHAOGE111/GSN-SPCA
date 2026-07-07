args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

opt <- parse_args(args)
required <- c("method", "work-dir", "k", "k-group", "we", "t", "niter", "err", "num-init", "edge-mode")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
}

method <- opt[["method"]]
work_dir <- opt[["work-dir"]]
edge_mode <- opt[["edge-mode"]]
hardcoded_edge_n <- if ("hardcoded-edge-n" %in% names(opt) && nzchar(opt[["hardcoded-edge-n"]])) {
  as.integer(opt[["hardcoded-edge-n"]])
} else {
  NA_integer_
}

params <- list(
  k = as.integer(opt[["k"]]),
  k_group = as.numeric(opt[["k-group"]]),
  we = as.numeric(opt[["we"]]),
  t = as.numeric(opt[["t"]]),
  niter = as.integer(opt[["niter"]]),
  err = as.numeric(opt[["err"]]),
  num_init = as.integer(opt[["num-init"]])
)

setwd(work_dir)
set.seed(1)

read_edges <- function(path, mode = "nrow", hardcoded_n = NA_integer_) {
  result <- read.table(path, quote = "\"", comment.char = "")
  if (mode == "hardcoded") {
    if (is.na(hardcoded_n)) {
      stop("edge_mode=hardcoded requires hardcoded_edge_n")
    }
    if (nrow(result) < hardcoded_n) {
      stop(sprintf("Edge file has %d rows, less than hardcoded_edge_n=%d", nrow(result), hardcoded_n))
    }
    edge_n <- hardcoded_n
  } else {
    edge_n <- nrow(result)
  }
  edges <- vector("list", edge_n)
  for (i in seq_len(edge_n)) {
    edges[[i]] <- c(result[i, 1], result[i, 2])
  }
  edges
}

read_groups <- function(path) {
  lines <- readLines(path, warn = FALSE)
  groups <- vector("list", length(lines))
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], "\t")[[1]]
    parts <- parts[nzchar(parts)]
    groups[[i]] <- as.integer(parts)
  }
  groups
}

write_outputs <- function(out) {
  write.table(out$U, file = "U_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(out$V, file = "V_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(out$D, file = "D_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
}

process_awge_power <- function(gene_file = "gene_new_o.txt", result_file = "result-1_p2.txt") {
  gene_count <- length(readLines(gene_file))
  pairs <- do.call(
    rbind,
    lapply(
      strsplit(readLines(result_file), "\\s+"),
      function(x) if (length(x) >= 2) as.numeric(x[1:2]) else NULL
    )
  )
  pairs_norm <- t(apply(pairs, 1, function(x) if (x[1] != x[2] && x[1] > x[2]) x[2:1] else x))
  unique_pairs <- unique(pairs_norm)
  freq <- table(c(unique_pairs[, 1], unique_pairs[unique_pairs[, 1] != unique_pairs[, 2], 2]))
  counts <- rep(0, gene_count)
  counts[as.numeric(names(freq))] <- as.numeric(freq)
  norm_val <- ifelse(counts > 0, 1 + (counts / max(counts)), 1)
  writeLines(as.character(norm_val), "power.txt")
}

process_dm_power <- function(gene_file = "gene_new_o.txt", output_file = "power.txt") {
  raw_data <- readLines(gene_file, warn = FALSE)
  raw_data <- trimws(raw_data)
  raw_data <- raw_data[nzchar(raw_data)]
  n_genes <- length(raw_data)
  results <- numeric(n_genes)
  pb <- txtProgressBar(min = 0, max = n_genes, style = 3)
  for (i in seq_len(n_genes)) {
    data <- suppressWarnings(as.numeric(unlist(strsplit(raw_data[i], "\\s+"))))
    data <- data[is.finite(data)]
    if (length(data) < 2) {
      results[i] <- 1.5
      setTxtProgressBar(pb, i)
      next
    }
    split_point <- max(1, floor(length(data) / 2))
    group1 <- data[1:split_point]
    group2 <- data[(split_point + 1):length(data)]
    p_value <- tryCatch({
      if (sd(group1) == 0 && sd(group2) == 0) { 1.0 }
      else {
        test <- t.test(group1, group2, var.equal = FALSE)
        max(min(test$p.value, 1.0), 0.0)
      }
    }, error = function(e) 1.0)
    results[i] <- max(min(round(2 - p_value, 4), 2.0), 1.0)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  writeLines(format(results, nsmall = 4, scientific = FALSE), output_file)
  cat(sprintf("DM power generated: %d genes\n", n_genes))
}

gene_new_o <- read.delim("gene_new_o.txt", header = FALSE)
gene <- t(data.matrix(gene_new_o))

if (method == "ESPCA") {
  source("fun_SPCA.R")
  source("fun_ESPCA.R")
  edges <- read_edges("result-1_p2.txt", edge_mode, hardcoded_edge_n)
  out <- ESPCA(
    X = gene,
    k = params$k,
    overlap.group = edges,
    k.group = params$k_group,
    we = params$we,
    t = params$t,
    niter = params$niter,
    err = params$err,
    Num.init = params$num_init
  )
  write_outputs(out)
} else if (method == "DM-ESPCA") {
  process_dm_power()
  source("fun_SPCA.R")
  source("fun_DM_ESPCA.R")
  edges <- read_edges("result-1_p2.txt", edge_mode, hardcoded_edge_n)
  out <- DMESPCA(
    X = gene,
    k = params$k,
    overlap.group = edges,
    k.group = params$k_group,
    we = params$we,
    t = params$t,
    niter = params$niter,
    err = params$err,
    Num.init = params$num_init
  )
  write_outputs(out)
} else if (method == "AWGE-ESPCA") {
  process_awge_power()
  source("fun_AWGE.R")
  edges <- read_edges("result-1_p2.txt", edge_mode, hardcoded_edge_n)
  out <- AWGE(
    X = gene,
    k = params$k,
    overlap.group = edges,
    k.group = params$k_group,
    we = params$we,
    t = params$t,
    niter = params$niter,
    err = params$err,
    Num.init = params$num_init
  )
  write_outputs(out)
} else if (method == "GSN-SPCA") {
  if (file.exists("result_max_values.txt")) {
    file.remove("result_max_values.txt")
  }
  source("fun_GSN-SPCA.R", local = TRUE)
  edges <- read_groups("result-1_p2.txt")
  out <- ESPCA(
    X = gene,
    k = params$k,
    overlap.group = edges,
    k.group = params$k_group,
    we = params$we,
    t = params$t,
    niter = params$niter,
    err = params$err,
    Num.init = params$num_init
  )
  write_outputs(out)
} else {
  stop(sprintf("Unknown method: %s", method))
}

cat(sprintf("completed %s in %s\n", method, work_dir))

