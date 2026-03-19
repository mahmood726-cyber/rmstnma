# Embed known datasets from installed packages using package-specific heuristics

library(utils)

dir.create('inst/extdata', recursive = TRUE, showWarnings = FALSE)

nodes_from <- function(tr){
  tr <- sort(unique(na.omit(as.character(tr))))
  data.frame(id = tr, name = tr, stringsAsFactors = FALSE)
}

write_tables <- function(id, arm = NULL, contrast = NULL, nodes = NULL){
  base <- file.path('inst','extdata', id)
  if (!is.null(arm))      write.csv(arm,      paste0(base, '_studies.csv'), row.names = FALSE)
  if (!is.null(contrast)) write.csv(contrast, paste0(base, '_contrast.csv'), row.names = FALSE)
  if (!is.null(nodes))    write.csv(nodes,    paste0(base, '_nodes.csv'), row.names = FALSE)
}

# Generic mappers
to_arm_from_pairwise_df <- function(df){
  nms <- tolower(names(df)); names(df) <- nms
  t1 <- (c('treat1','t1','treatment1','tx1')[c('treat1','t1','treatment1','tx1')%in%nms])[1]
  t2 <- (c('treat2','t2','treatment2','tx2')[c('treat2','t2','treatment2','tx2')%in%nms])[1]
  t3 <- (c('treat3','t3','treatment3','tx3')[c('treat3','t3','treatment3','tx3')%in%nms])[1]
  stud <- (c('studlab','study','trial','author','id')[c('studlab','study','trial','author','id')%in%nms])[1]
  # binary
  ev1 <- (c('event1','responders1','r1')[c('event1','responders1','r1')%in%nms])[1]
  ev2 <- (c('event2','responders2','r2')[c('event2','responders2','r2')%in%nms])[1]
  ev3 <- (c('event3','responders3','r3')[c('event3','responders3','r3')%in%nms])[1]
  n1  <- (c('n1','total1','sample1','N1')[c('n1','total1','sample1','N1')%in%nms])[1]
  n2  <- (c('n2','total2','sample2','N2')[c('n2','total2','sample2','N2')%in%nms])[1]
  n3  <- (c('n3','total3','sample3','N3')[c('n3','total3','sample3','N3')%in%nms])[1]
  # continuous
  y1  <- (c('y1','mean1')[c('y1','mean1')%in%nms])[1]
  y2  <- (c('y2','mean2')[c('y2','mean2')%in%nms])[1]
  y3  <- (c('y3','mean3')[c('y3','mean3')%in%nms])[1]
  sd1 <- (c('sd1','stdev1')[c('sd1','stdev1')%in%nms])[1]
  sd2 <- (c('sd2','stdev2')[c('sd2','stdev2')%in%nms])[1]
  sd3 <- (c('sd3','stdev3')[c('sd3','stdev3')%in%nms])[1]

  if (is.na(t1) || is.na(t2)) return(NULL)

  S <- if (!is.na(stud)) as.character(df[[stud]]) else as.character(seq_len(nrow(df)))
  rows <- list()
  # decide mode per row: prefer binary if event cols present; else continuous
  bin_ok <- !is.na(ev1) && !is.na(n1)
  cont_ok <- !is.na(y1) && (!is.na(sd1) || !is.na(n1))
  if (!bin_ok && !cont_ok) return(NULL)

  if (bin_ok) {
    rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t1]]), r=as.numeric(df[[ev1]]), n=as.numeric(df[[n1]]), stringsAsFactors=FALSE)
    rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t2]]), r=as.numeric(df[[ev2]]), n=as.numeric(df[[n2]]), stringsAsFactors=FALSE)
    if (!is.na(t3) && !is.na(ev3) && !is.na(n3)) rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t3]]), r=as.numeric(df[[ev3]]), n=as.numeric(df[[n3]]), stringsAsFactors=FALSE)
  } else {
    # continuous
    rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t1]]), y=as.numeric(df[[y1]]), sd=as.numeric(df[[sd1]]), n=as.numeric(df[[n1]]), stringsAsFactors=FALSE)
    rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t2]]), y=as.numeric(df[[y2]]), sd=as.numeric(df[[sd2]]), n=as.numeric(df[[n2]]), stringsAsFactors=FALSE)
    if (!is.na(t3) && !is.na(y3)) rows[[length(rows)+1]] <- data.frame(study=S, treatment=as.character(df[[t3]]), y=as.numeric(df[[y3]]), sd=as.numeric(df[[sd3]]), n=as.numeric(df[[n3]]), stringsAsFactors=FALSE)
  }
  arm <- do.call(rbind, rows)
  nodes <- nodes_from(arm$treatment)
  list(arm=arm, nodes=nodes)
}

to_contrast_from_effects_df <- function(df){
  nms <- tolower(names(df)); names(df) <- nms
  t1 <- (c('treat1','t1','treatment1','tx1')[c('treat1','t1','treatment1','tx1')%in%nms])[1]
  t2 <- (c('treat2','t2','treatment2','tx2')[c('treat2','t2','treatment2','tx2')%in%nms])[1]
  te <- (c('te','lnor','logor','effect','y')[c('te','lnor','logor','effect','y')%in%nms])[1]
  se <- (c('sete','selnor','se','stderr')[c('sete','selnor','se','stderr')%in%nms])[1]
  stud <- (c('studlab','study','trial','author','id')[c('studlab','study','trial','author','id')%in%nms])[1]
  if (is.na(t1) || is.na(t2) || is.na(te) || is.na(se)) return(NULL)
  df2 <- data.frame(
    study = if (!is.na(stud)) as.character(df[[stud]]) else NA_character_,
    treat1 = as.character(df[[t1]]),
    treat2 = as.character(df[[t2]]),
    TE = as.numeric(df[[te]]),
    seTE = as.numeric(df[[se]]),
    stringsAsFactors = FALSE
  )
  nodes <- nodes_from(c(df2$treat1, df2$treat2))
  list(contrast=df2, nodes=nodes)
}

# Package-specific embedders
embed_netmeta <- function(id, obj){
  if (is.data.frame(obj)){
    out <- to_arm_from_pairwise_df(obj)
    if (is.null(out)) out <- to_contrast_from_effects_df(obj)
    if (!is.null(out)) { write_tables(id, arm=out$arm, contrast=out$contrast, nodes=out$nodes); return(TRUE) }
  }
  FALSE
}

embed_multinma <- function(id, obj){
  if (is.list(obj)){
    if (!is.null(obj$agd_arm)){
      df <- as.data.frame(obj$agd_arm)
      names(df) <- tolower(names(df))
      # expect study, treatment and either r/n or mean/sd/n
      arm <- df
      nodes <- nodes_from(df$treatment)
      write_tables(id, arm=arm, nodes=nodes)
      return(TRUE)
    }
    if (!is.null(obj$agd_contrast)){
      df <- as.data.frame(obj$agd_contrast)
      names(df) <- tolower(names(df))
      out <- to_contrast_from_effects_df(df)
      if (!is.null(out)) { write_tables(id, contrast=out$contrast, nodes=out$nodes); return(TRUE) }
    }
    # fallback: first df inside list
    dfs <- Filter(is.data.frame, obj)
    for (d in dfs){
      out <- to_arm_from_pairwise_df(d); if (!is.null(out)) { write_tables(id, arm=out$arm, nodes=out$nodes); return(TRUE) }
      out2 <- to_contrast_from_effects_df(d); if (!is.null(out2)) { write_tables(id, contrast=out2$contrast, nodes=out2$nodes); return(TRUE) }
    }
  }
  FALSE
}

embed_generic_list <- function(id, obj){
  if (is.list(obj)){
    dfs <- Filter(is.data.frame, obj)
    for (d in dfs){
      out <- to_arm_from_pairwise_df(d); if (!is.null(out)) { write_tables(id, arm=out$arm, nodes=out$nodes); return(TRUE) }
      out2 <- to_contrast_from_effects_df(d); if (!is.null(out2)) { write_tables(id, contrast=out2$contrast, nodes=out2$nodes); return(TRUE) }
    }
  }
  FALSE
}

get_from_pkg <- function(pkg, obj){
  e <- new.env(parent=emptyenv())
  suppressWarnings(utils::data(list=obj, package=pkg, envir=e))
  if (!exists(obj, envir=e, inherits=FALSE)) {
    obj2 <- sub("\\s*\\(.*$", "", obj)
    suppressWarnings(utils::data(list=obj2, package=pkg, envir=e))
    if (exists(obj2, envir=e, inherits=FALSE)) return(get(obj2, envir=e, inherits=FALSE))
    return(NULL)
  }
  get(obj, envir=e, inherits=FALSE)
}

reg <- read.csv('inst/registry/registry.csv', stringsAsFactors = FALSE)

targets <- subset(reg, package %in% c('netmeta','multinma','MBNMAdose','pcnetmeta','bnma'))

embedded <- 0L
for (i in seq_len(nrow(targets))){
  id <- targets$id[i]; pkg <- targets$package[i]; obj <- targets$object[i]
  if (file.exists(file.path('inst','extdata', paste0(id, '_studies.csv'))) ||
      file.exists(file.path('inst','extdata', paste0(id, '_contrast.csv')))) next
  x <- get_from_pkg(pkg, obj)
  if (is.null(x)) next
  ok <- FALSE
  if (pkg == 'netmeta') ok <- embed_netmeta(id, x)
  else if (pkg == 'multinma') ok <- embed_multinma(id, x)
  else ok <- embed_generic_list(id, x)
  if (ok) embedded <- embedded + 1L
}
cat('Embedded via recipes: ', embedded, '\n', sep='')
