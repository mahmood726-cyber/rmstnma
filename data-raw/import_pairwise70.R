# Import Cochrane pairwise/NMA datasets from an external repo (read-only)
# without modifying that repo. Heuristically standardize to arm/contrast and
# write into inst/extdata here.

args <- commandArgs(trailingOnly = TRUE)
src_dir <- if (length(args) >= 1) args[[1]] else "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"

message("Source dir: ", src_dir)
if (!dir.exists(src_dir)) stop("Source directory not found: ", src_dir)

dir.create(file.path('inst','extdata'), showWarnings = FALSE, recursive = TRUE)

normalize_names <- function(x) {
  nm <- tolower(gsub("[^a-zA-Z0-9]+", "_", names(x)))
  # de-duplicate underscores
  nm <- gsub("_+", "_", nm)
  nm <- sub("^_", "", nm)
  names(x) <- nm
  x
}

pick <- function(df, candidates) {
  nms <- names(df)
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) NULL else df[[hit[1]]]
}

standardize_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- normalize_names(df)
  nms <- names(df)

  # Try arm-level binary (experimental/control)
  exp_n <- pick(df, c('experimental_n','exp_n','treated_n','intervention_n'))
  con_n <- pick(df, c('control_n','cntl_n','comparator_n','placebo_n'))
  exp_e <- pick(df, c('experimental_events','exp_events','treated_events','intervention_events','events_e','experimental_cases'))
  con_e <- pick(df, c('control_events','cntl_events','comparator_events','placebo_events','events_c','control_cases'))
  exp_trt <- pick(df, c('experimental_treatment','experimental_intervention','exp_treatment','treated','intervention','experimental','intervention_name','experimental_arm','exp_arm'))
  con_trt <- pick(df, c('control_treatment','control_intervention','comparator','placebo','control','comparator_name','control_arm','cntl_arm'))
  study <- pick(df, c('study','study_id','studylab','trial','author_year','author'))
  if (!is.null(exp_n) && !is.null(con_n) && !is.null(exp_trt) && !is.null(con_trt)) {
    S <- if (!is.null(study)) as.character(study) else as.character(seq_len(nrow(df)))
    arm <- rbind(
      data.frame(study = S, treatment = as.character(exp_trt), r = if (!is.null(exp_e)) as.numeric(exp_e) else NA_real_, n = as.numeric(exp_n), stringsAsFactors = FALSE),
      data.frame(study = S, treatment = as.character(con_trt), r = if (!is.null(con_e)) as.numeric(con_e) else NA_real_, n = as.numeric(con_n), stringsAsFactors = FALSE)
    )
    return(list(type = 'arm', data = arm))
  }

  # Try arm-level continuous (means/sd/n)
  exp_mean <- pick(df, c('experimental_mean','exp_mean','treated_mean','intervention_mean','mean_e'))
  con_mean <- pick(df, c('control_mean','cntl_mean','comparator_mean','placebo_mean','mean_c'))
  exp_sd <- pick(df, c('experimental_sd','exp_sd','treated_sd','intervention_sd','sd_e','se_e'))
  con_sd <- pick(df, c('control_sd','cntl_sd','comparator_sd','placebo_sd','sd_c','se_c'))
  if (!is.null(exp_mean) && !is.null(con_mean) && !is.null(exp_n) && !is.null(con_n) && !is.null(exp_trt) && !is.null(con_trt)) {
    if (is.null(study)) study <- seq_len(nrow(df))
    arm <- rbind(
      data.frame(study = as.character(study), treatment = as.character(exp_trt), y = as.numeric(exp_mean), sd = as.numeric(exp_sd), n = as.numeric(exp_n), stringsAsFactors = FALSE),
      data.frame(study = as.character(study), treatment = as.character(con_trt), y = as.numeric(con_mean), sd = as.numeric(con_sd), n = as.numeric(con_n), stringsAsFactors = FALSE)
    )
    return(list(type = 'arm', data = arm))
  }

  # Try contrast-level if TE + seTE available or CI to compute se
  te <- pick(df, c('te','effect','lnor','logor','y','mean'))
  sete <- pick(df, c('sete','se','stderr'))
  lcl <- pick(df, c('ci_low','lcl','lower_ci','ci_start'))
  ucl <- pick(df, c('ci_high','ucl','upper_ci','ci_end'))
  t1 <- pick(df, c('treat1','t1','treatment1','experimental_treatment'))
  t2 <- pick(df, c('treat2','t2','treatment2','control_treatment'))
  if (!is.null(te) && (!is.null(sete) || (!is.null(lcl) && !is.null(ucl)))) {
    if (is.null(sete)) sete <- (as.numeric(ucl) - as.numeric(lcl)) / (2*1.96)
    # Default treatment labels if none present
    if (is.null(t1)) t1 <- rep('Experimental', length.out = nrow(df))
    if (is.null(t2)) t2 <- rep('Control', length.out = nrow(df))
    con <- data.frame(study = if (!is.null(study)) as.character(study) else NA_character_,
                      treat1 = as.character(t1), treat2 = as.character(t2),
                      TE = as.numeric(te), seTE = as.numeric(sete), stringsAsFactors = FALSE)
    return(list(type = 'contrast', data = con))
  }

  NULL
}

write_one <- function(id, obj) {
  norm <- standardize_df(obj)
  if (is.null(norm)) return(FALSE)
  if (norm$type == 'arm') {
    utils::write.csv(norm$data, file.path('inst','extdata', paste0(id, '_studies.csv')), row.names = FALSE)
  } else if (norm$type == 'contrast') {
    utils::write.csv(norm$data, file.path('inst','extdata', paste0(id, '_contrast.csv')), row.names = FALSE)
  }
  # Nodes
  if ('treatment' %in% names(norm$data)) {
    tr <- sort(unique(na.omit(as.character(norm$data$treatment))))
  } else if (all(c('treat1','treat2') %in% names(norm$data))) {
    tr <- sort(unique(na.omit(c(as.character(norm$data$treat1), as.character(norm$data$treat2)))))
  } else tr <- character()
  if (length(tr)) utils::write.csv(data.frame(id = tr, name = tr, stringsAsFactors = FALSE),
                                   file.path('inst','extdata', paste0(id, '_nodes.csv')), row.names = FALSE)
  TRUE
}

files <- list.files(src_dir, pattern = '\\.rda$', full.names = TRUE)
ok <- 0L; fail <- 0L
for (f in files) {
  e <- new.env(parent = emptyenv())
  try(load(f, envir = e), silent = TRUE)
  nms <- ls(envir = e)
  got <- FALSE
  for (nm in nms) {
    x <- get(nm, envir = e)
    if (is.data.frame(x)) {
      id <- paste0('cochrane__', tools::file_path_sans_ext(basename(f)))
      if (write_one(id, x)) { ok <- ok + 1L; got <- TRUE; break }
    }
  }
  if (!got) fail <- fail + 1L
}
cat('Pairwise70 import complete. Files processed:', length(files), ' ok:', ok, ' fail:', fail, '\n')
