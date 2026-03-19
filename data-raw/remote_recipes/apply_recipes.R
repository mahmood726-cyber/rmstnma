# Remote recipes: light transforms to coerce common shapes into
# standard arm/contrast tables before normalization.

# Apply a set of heuristics/recipes. Return a transformed object (data.frame or list)
# or NULL if no recipe matched.
apply_remote_recipe <- function(id, obj) {
  # Only handle data.frames here; lists will go through normal heuristics later
  if (!is.data.frame(obj)) return(NULL)
  df <- obj
  nms <- tolower(names(df))

  # 0) Simple rename synonyms to standard arm columns
  syn <- nms
  syn[nms == 'studlab'] <- 'study'
  syn[nms %in% c('treat','tx','arm')] <- 'treatment'
  syn[nms %in% c('responders','events','event','r')] <- 'r'
  syn[nms %in% c('samplesize','sample','total','n_total','n1','n2','n3')] <- 'n'
  names(df) <- syn; nms <- tolower(names(df))

  # 0b) If TE/seTE variants exist, normalize to TE/seTE
  if ('te.w' %in% nms && !('te' %in% nms)) names(df)[match('te.w', nms)] <- 'te'
  if ('sete.w' %in% nms && !('sete' %in% nms)) names(df)[match('sete.w', nms)] <- 'sete'
  nms <- tolower(names(df))

  # 1) Pairwise binary wide -> long (event1/n1/event2/n2 + treat1/treat2)
  if (all(c('treat1','treat2') %in% nms) && all(c('event1','n1','event2','n2') %in% nms)) {
    dff <- df; names(dff) <- nms
    S <- if ('study' %in% nms) as.character(dff$study) else as.character(seq_len(nrow(dff)))
    long <- rbind(
      data.frame(study = S, treatment = as.character(dff$treat1), r = as.numeric(dff$event1), n = as.numeric(dff$n1), stringsAsFactors = FALSE),
      data.frame(study = S, treatment = as.character(dff$treat2), r = as.numeric(dff$event2), n = as.numeric(dff$n2), stringsAsFactors = FALSE)
    )
    return(long)
  }

  # 1b) Generalized k-arm binary wide -> long using numbered suffixes r1/n1, r2/n2, ... + t1/t2/...
  rcols <- grep('^(event|responders|r)[0-9]+$', nms, value = TRUE)
  ncols <- grep('^n[0-9]+$', nms, value = TRUE)
  tcols <- grep('^(treat|trt|tx|t)[0-9]+$', nms, value = TRUE)
  if (length(rcols) >= 2 && length(ncols) >= 2 && length(tcols) >= 2) {
    dff <- df; names(dff) <- nms
    idx <- sort(as.integer(gsub('^.*?([0-9]+)$','\u0001', rcols)))
    build_long <- function(irow){
      out <- list(); S <- if ('study' %in% nms) as.character(dff$study[irow]) else as.character(irow)
      for (k in idx) {
        rk <- paste0(sub('[0-9]+$','', rcols[1]), k)
        nk <- paste0('n', k)
        tk <- paste0(sub('[0-9]+$','', tcols[1]), k)
        if (rk %in% names(dff) && nk %in% names(dff) && tk %in% names(dff)) {
          out[[length(out)+1]] <- data.frame(study = S,
                                             treatment = as.character(dff[[tk]][irow]),
                                             r = as.numeric(dff[[rk]][irow]),
                                             n = as.numeric(dff[[nk]][irow]),
                                             stringsAsFactors = FALSE)
        }
      }
      if (length(out)) do.call(rbind, out) else NULL
    }
    rows <- lapply(seq_len(nrow(dff)), build_long)
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) return(do.call(rbind, rows))
  }

  # 2) Pairwise continuous wide -> contrast (y1,y2,se1,se2[,y3,se3]; t1,t2[,t3])
  if (all(c('y1','y2','se1','se2','t1','t2') %in% nms)) {
    dff <- df; names(dff) <- nms
    build_rows <- function(row) {
      out <- list()
      if (!is.na(row$y1) && !is.na(row$y2) && !is.na(row$se1) && !is.na(row$se2)) {
        out[[length(out)+1]] <- data.frame(study = NA_character_,
                                           treat1 = as.character(row$t1),
                                           treat2 = as.character(row$t2),
                                           TE = as.numeric(row$y2) - as.numeric(row$y1),
                                           seTE = sqrt(as.numeric(row$se1)^2 + as.numeric(row$se2)^2),
                                           stringsAsFactors = FALSE)
      }
      if (all(c('y3','se3','t3') %in% names(row))) {
        if (!is.na(row$y3) && !is.na(row$se3)) {
          out[[length(out)+1]] <- data.frame(study = NA_character_,
                                             treat1 = as.character(row$t1),
                                             treat2 = as.character(row$t3),
                                             TE = as.numeric(row$y3) - as.numeric(row$y1),
                                             seTE = sqrt(as.numeric(row$se1)^2 + as.numeric(row$se3)^2),
                                             stringsAsFactors = FALSE)
        }
      }
      if (length(out)) do.call(rbind, out) else NULL
    }
    rows <- lapply(seq_len(nrow(dff)), function(i) build_rows(dff[i, , drop = FALSE]))
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) return(do.call(rbind, rows))
  }

  # 3) BUGSnet-like binary with Study/Treatment/n and an event column
  if (all(c('study','treatment','n') %in% nms) && any(c('diabetes','event','events','r') %in% nms)) {
    evcol <- intersect(c('diabetes','event','events','r'), nms)[1]
    dff <- df; names(dff) <- nms
    long <- data.frame(study = as.character(dff$study),
                       treatment = as.character(dff$treatment),
                       r = as.numeric(dff[[evcol]]),
                       n = as.numeric(dff$n),
                       stringsAsFactors = FALSE)
    return(long)
  }

  # 4) Simple arm-level if study/treatment/r/n present via synonyms
  if (all(c('study','treatment','r','n') %in% nms)) {
    dff <- df; names(dff) <- nms
    return(data.frame(study = as.character(dff$study),
                      treatment = as.character(dff$treatment),
                      r = as.numeric(dff$r), n = as.numeric(dff$n),
                      stringsAsFactors = FALSE))
  }

  # 5) Contrast if treat1/treat2 and TE/seTE variants exist
  if (any(c('treat1','t1') %in% nms) && any(c('treat2','t2') %in% nms) && any(c('te','lnor','effect') %in% nms) && any(c('sete','se','stderr','selnor') %in% nms)) {
    dff <- df; names(dff) <- nms
    t1 <- dff[[ intersect(c('treat1','t1'), names(dff))[1] ]]
    t2 <- dff[[ intersect(c('treat2','t2'), names(dff))[1] ]]
    te <- dff[[ intersect(c('te','lnor','effect'), names(dff))[1] ]]
    se <- dff[[ intersect(c('sete','se','stderr','selnor'), names(dff))[1] ]]
    con <- data.frame(study = if ('study' %in% names(dff)) as.character(dff$study) else NA_character_,
                      treat1 = as.character(t1), treat2 = as.character(t2),
                      TE = as.numeric(te), seTE = as.numeric(se), stringsAsFactors = FALSE)
    return(con)
  }

  # 6) BUGSnet contrast-like: differences + se.diffs + coded treatment pairs per study
  if (all(c('differences','se.diffs') %in% nms) && 'treatment' %in% nms && 'study' %in% nms) {
    dff <- df; names(dff) <- nms
    # Attempt to pair rows within each study: even rows vs odd rows
    dff$.row <- seq_len(nrow(dff))
    con <- do.call(rbind, lapply(split(dff, dff$study), function(dd){
      if (nrow(dd) < 2) return(NULL)
      dd <- dd[order(dd$.row), ]
      i1 <- seq(1, nrow(dd)-1, by = 2); i2 <- i1 + 1
      data.frame(study = as.character(dd$study[i1]),
                 treat1 = as.character(dd$treatment[i1]),
                 treat2 = as.character(dd$treatment[i2]),
                 TE = as.numeric(dd$differences[i2]),
                 seTE = as.numeric(dd$`se.diffs`[i2]),
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(con) && nrow(con)) return(con)
  }

  NULL
}
