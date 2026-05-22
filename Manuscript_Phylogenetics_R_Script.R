#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  required_pkgs <- c("Biostrings","data.table","stringr","ape","phangorn","R.utils")
  for (p in required_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(Biostrings); library(data.table); library(stringr); library(ape); library(phangorn); library(R.utils)
})

# ---------------------------
# USER CONFIG (EDIT)
# ---------------------------
base_dir <- "/Users/jamesking/My_Files/University_of_Oxford/MBiol_Biology_Year_4/RESEARCH/Phylogenetic_Analysis"
outdir   <- file.path(base_dir, "results_hif_pathway_analysis_pipeline")
pfam_dir <- file.path(base_dir, "pfam_db")
threads  <- 4L
incE_default <- 1e-2
keep_old_results <- FALSE

family_params <- list(
  HIF1a = list(incE = 1e-2, min_len = 300L, max_len = 1100L),
  PHD2  = list(incE = 1e-2, min_len = 300L, max_len = 600L),
  VHL   = list(incE = 1e-1, min_len = 100L, max_len = 250L)
)
# ---------------------------

msg <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")
warnmsg <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)
stopmsg <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
safe_normalize <- function(p) tryCatch(normalizePath(p, winslash = "/", mustWork = TRUE), error = function(e) p)

run <- function(cmd, args = character(), to_file = NULL, echo = TRUE, stop_on_fail = FALSE) {
  if (echo) cat("[RUN]", cmd, paste(shQuote(as.character(args)), collapse = " "), "\n")
  exe <- Sys.which(cmd)
  if (!nzchar(exe)) {
    txt <- sprintf("[ERROR] Command not found in PATH: %s", cmd)
    if (stop_on_fail) stop(txt) else { warning(txt); return(invisible(NA)) }
  }
  if (!is.null(to_file)) {
    dir.create(dirname(to_file), recursive = TRUE, showWarnings = FALSE)
    errf <- paste0(to_file, ".stderr.log")
    status <- tryCatch(system2(exe, args = as.character(args), stdout = to_file, stderr = errf), error = function(e) e)
    if (inherits(status, "error")) {
      msg <- sprintf("[ERROR] Failed to run '%s': %s", cmd, conditionMessage(status))
      if (stop_on_fail) stop(msg) else { warning(msg); return(invisible(NA)) }
    }
    if (is.numeric(status) && status == 0) return(invisible(0L))
    msg <- sprintf("[WARN] Command '%s' exited with status %s. See %s and %s", cmd, paste(status, collapse=","), to_file, errf)
    if (stop_on_fail) stop(msg) else { warning(msg); return(invisible(status)) }
  } else {
    out <- tryCatch(system2(exe, args = as.character(args), stdout = TRUE, stderr = TRUE), error = function(e) e)
    if (inherits(out, "error")) {
      msg <- sprintf("[ERROR] Failed to run '%s': %s", cmd, conditionMessage(out))
      if (stop_on_fail) stop(msg) else { warning(msg); return(invisible(NA)) }
    }
    return(out)
  }
}

anchors_candidates <- list.files(base_dir, pattern = "Anchor.*(fa|fasta|faa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
anchors_path <- if (length(anchors_candidates)) anchors_candidates[[1]] else stopmsg("Anchor FASTA not found under base_dir. Place an anchors FASTA under base_dir.")
print(anchors_candidates)
print(anchors_path)

ncbi_proteomes_dir <- file.path(base_dir, "NCBI_Proteomes")
if (!dir.exists(ncbi_proteomes_dir)) stopmsg("NCBI_Proteomes directory not found under base_dir: %s", ncbi_proteomes_dir)

proteomes_paths <- list.files(ncbi_proteomes_dir, pattern="\\.protein\\.(fa|fasta|faa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(proteomes_paths)) {
  proteomes_paths <- list.files(ncbi_proteomes_dir, pattern="(protein|proteome).*\\.(fa|fasta|faa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}
if (!length(proteomes_paths)) stopmsg("No proteome FASTA files found under NCBI_Proteomes (%s)", ncbi_proteomes_dir)

anchors_path <- safe_normalize(anchors_path)
proteomes_paths <- vapply(proteomes_paths, safe_normalize, FUN.VALUE = character(1))

ncbi_proteomes_dir <- normalizePath(ncbi_proteomes_dir, winslash = "/", mustWork = TRUE)

in_ncbi <- vapply(proteomes_paths, function(p) {
  ok <- FALSE
  tryCatch({
    pn <- normalizePath(p, winslash = "/", mustWork = TRUE)
    ok <- startsWith(pn, ncbi_proteomes_dir)
  }, error = function(e) { ok <- FALSE })
  ok
}, logical(1))

exclude_pattern <- "(?i)(wormbase|wbps|worm_base|WBPS)"
not_excluded <- !grepl(exclude_pattern, proteomes_paths, perl = TRUE)
proteomes_paths <- proteomes_paths[in_ncbi & not_excluded]
if (!length(proteomes_paths)) stopmsg("No proteome FASTA files found in NCBI_Proteomes after filtering. Confirm files and patterns.")

allowed_species <- basename(proteomes_paths)
out_subdirs <- c("hmmsearch", "candidates", "pfam")
for (d in out_subdirs) {
  dpath <- file.path(outdir, d)
  if (!dir.exists(dpath)) next
  files <- list.files(dpath, recursive = TRUE, full.names = TRUE)
  if (!length(files)) next
  files_sp <- vapply(basename(files), function(fn) sub("\\.([^.]+)$", "", fn), character(1))
  to_remove <- files[! files_sp %in% allowed_species]
  if (length(to_remove)) {
    for (f in to_remove) {
      tryCatch({
        file.remove(f)
        p <- dirname(f)
        if (length(list.files(p)) == 0) unlink(p, recursive = TRUE, force = TRUE)
      }, error = function(e) warning("Failed to remove stale output file: ", f))
    }
  }
}

msg("Anchors: %s", anchors_path)
msg("Proteomes found: %d", length(proteomes_paths))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dirs <- file.path(outdir, c("anchors","hmms","alignments","hmmsearch","pfam","candidates","phylogeny","tables","report","logs"))
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

required_bins <- c("mafft","hmmbuild","hmmsearch","hmmscan","hmmpress","iqtree")
missing_bins <- required_bins[!nzchar(Sys.which(required_bins))]
if (length(missing_bins)) stopmsg("Missing required binaries: %s. Install them and ensure they are discoverable.", paste(missing_bins, collapse=", "))
msg("[OK] external binaries present")

get_version_robust <- function(cmd) {
  path <- Sys.which(cmd)
  if (!nzchar(path)) return(sprintf("%s: NOT FOUND", cmd))
  candidates <- c("--version", "-V", "-v", "-h", "-help")
  for (carg in candidates) {
    out <- tryCatch({ suppressWarnings(system2(path, args = carg, stdout = TRUE, stderr = TRUE)) }, error = function(e) NA_character_)
    if (is.character(out) && length(out) > 0) {
      txt <- paste(out, collapse = "\n")
      lines <- unlist(strsplit(txt, "\n"))
      pick <- grep("version|Version|v[0-9]+\\.[0-9]+|HMMER|IQ-TREE|MAFFT", lines, ignore.case = TRUE, value = TRUE)
      if (length(pick) > 0) { return(sprintf("%s: %s", cmd, trimws(pick[1]))) } else {
        first_nonempty <- lines[which(nzchar(trimws(lines)))[1]]
        if (!is.na(first_nonempty) && nzchar(first_nonempty)) return(sprintf("%s: %s", cmd, trimws(first_nonempty)))
      }
    }
  }
  out0 <- tryCatch({ suppressWarnings(system2(path, args = character(), stdout = TRUE, stderr = TRUE)) }, error = function(e) NA_character_)
  if (is.character(out0) && length(out0) > 0) {
    txt0 <- paste(out0, collapse = "\n")
    first_line <- unlist(strsplit(txt0, "\n"))[1]
    return(sprintf("%s: %s", cmd, trimws(first_line)))
  }
  sprintf("%s: present at %s (no textual version found)", cmd, path)
}

vinfo <- sapply(required_bins, get_version_robust, USE.NAMES = TRUE)
writeLines(c(paste("R:", R.version.string), vinfo), file.path(outdir, "report", "software_versions.txt"))
msg("[OK] software versions logged")

get_pfam_db_strict <- function(dest_dir = pfam_dir,
                               url = "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz") {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dest_dir, "Pfam-A.hmm")
  gz <- paste0(f, ".gz")
  if (file.exists(f)) {
    idxs <- paste0(f, c(".h3f",".h3i",".h3m",".h3p"))
    if (all(file.exists(idxs))) return(normalizePath(f))
    msg("Pfam present but indexes missing; running hmmpress")
    run("hmmpress", c(f), stop_on_fail = TRUE)
    if (all(file.exists(idxs))) return(normalizePath(f))
    stopmsg("hmmpress failed to create index files for Pfam-A.hmm")
  }
  msg("Downloading Pfam-A.hmm.gz to %s", gz)
  tryCatch(utils::download.file(url, gz, mode = "wb", quiet = FALSE), error = function(e) warning("download.file failed: ", conditionMessage(e)))
  if (!file.exists(gz)) {
    if (nzchar(Sys.which("curl"))) system(sprintf("curl -L -o %s %s", shQuote(gz), shQuote(url)))
    else if (nzchar(Sys.which("wget"))) system(sprintf("wget -O %s %s", shQuote(gz), shQuote(url)))
    else stopmsg("Cannot download Pfam-A.hmm.gz; place it in %s and retry", dest_dir)
  }
  if (!file.exists(gz)) stopmsg("Pfam download failed; place Pfam-A.hmm.gz in %s", dest_dir)
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils", repos = "https://cloud.r-project.org")
  R.utils::gunzip(gz, overwrite = TRUE)
  if (!file.exists(f)) stopmsg("gunzip did not produce Pfam-A.hmm")
  run("hmmpress", c(f), stop_on_fail = TRUE)
  idxs <- paste0(f, c(".h3f",".h3i",".h3m",".h3p"))
  if (!all(file.exists(idxs))) stopmsg("hmmpress did not create index files for Pfam-A.hmm")
  normalizePath(f)
}

pfam_db <- get_pfam_db_strict()
msg("[OK] Pfam DB at %s", pfam_db)

read_tblout <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  if (!length(lines)) return(NULL)
  parts <- strsplit(trimws(lines), "\\s+")
  targets <- vapply(parts, function(p) if (length(p)>=1) p[1] else NA_character_, character(1))
  data.frame(target = targets, stringsAsFactors = FALSE)
}

# read_domtbl_simple extracts domain-level i-Evalue and bitscore (i-score)
# Column indices assume HMMER3 domtblout order where i-Evalue at column 13 and i-score (bits) at column 14.
read_domtbl_simple <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  if (!length(lines)) return(NULL)
  parts <- strsplit(trimws(lines), "\\s+")
  extract <- function(p, i) if (length(p) >= i) p[[i]] else NA_character_
  df <- data.frame(
    target = vapply(parts, extract, character(1), i = 1),
    pfam_name = vapply(parts, extract, character(1), i = 4),
    pfam_acc = vapply(parts, extract, character(1), i = 5),
    # dom_i_evalue (c-Evalue) and i_evalue taken from typical domtblout columns
    dom_c_evalue = suppressWarnings(as.numeric(vapply(parts, extract, character(1), i = 12))),
    dom_i_evalue = suppressWarnings(as.numeric(vapply(parts, extract, character(1), i = 13))),
    bitscore = suppressWarnings(as.numeric(vapply(parts, extract, character(1), i = 14))),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$target) & nzchar(df$target), , drop = FALSE]
  unique(df)
}

# Filter domtbl rows to require paper thresholds:
# i-E-value <= 1e-5 AND bitscore >= 10
filter_domtbl_by_threshold <- function(df, evalue_thresh = 1e-5, bitscore_thresh = 10) {
  if (is.null(df) || !nrow(df)) return(df[FALSE, , drop = FALSE])
  # (require numeric values)
  keep <- !is.na(df$dom_i_evalue) & !is.na(df$bitscore) & (df$dom_i_evalue <= evalue_thresh) & (df$bitscore >= bitscore_thresh)
  df[keep, , drop = FALSE]
}

read_fasta_to_dt <- function(fa) {
  if (!file.exists(fa)) return(data.table(id=character(0), seq=character(0)))
  s <- readAAStringSet(fa)
  ids <- vapply(names(s), function(h) strsplit(h, "\\s+")[[1]][1], character(1))
  data.table(id = ids, seq = as.character(s))
}

extract_by_ids <- function(fa, ids, outfa) {
  if (!file.exists(fa) || !length(ids)) return(0L)
  s <- readAAStringSet(fa)
  hdr_tokens <- vapply(names(s), function(h) strsplit(h, "\\s+")[[1]][1], character(1))
  keep <- which(hdr_tokens %in% ids)
  if (length(keep)) { writeXStringSet(s[keep], outfa); return(length(keep)) }
  keep2 <- unique(unlist(lapply(ids, function(id) grep(id, names(s), fixed = TRUE))))
  if (length(keep2)) { writeXStringSet(s[keep2], outfa); return(length(keep2)) }
  0L
}

odd_scan <- function(seq) {
  m <- gregexpr("AP", seq, ignore.case = TRUE)[[1]]
  if (length(m) == 0 || m[1] == -1) return(list(n = 0L, pos = integer(0)))
  n <- nchar(seq); pos <- m
  keep <- (pos > 1) & ((pos + 1) < n)
  pos <- pos[keep]; list(n = length(pos), pos = pos)
}

hif_regex  <- "(?i)(HIF[\\s_-]*1(\\b|[^0-9])|Hypoxia[\\s-]*inducible[\\s-]*factor[^\\n]*1[^\\n]*(alpha|α)|HIF1A|HIF-1alpha|HIF1-alpha)"
phd2_regex <- "(?i)(PHD2|EGLN1|EGL-?9\\b|Prolyl[- ]hydroxylase[- ]domain[- ]containing\\s*protein\\s*2|EGLN-1)"
vhl_regex  <- "(?i)(VHL|von\\s*Hippel[- ]Lindau|von[- ]Hippel[- ]Lindau|vonHippelLindau)"

anchors <- tryCatch(readAAStringSet(anchors_path), error = function(e) stopmsg("Cannot read anchors FASTA: %s", conditionMessage(e)))
if (!length(anchors)) stopmsg("Anchor FASTA appears empty")
ids <- names(anchors)
hif_idx  <- grepl(hif_regex, ids, perl = TRUE)
phd_idx  <- grepl(phd2_regex, ids, perl = TRUE)
vhl_idx  <- grepl(vhl_regex, ids, perl = TRUE)

anchors_out_dir <- normalizePath(file.path(outdir, "anchors"), winslash = "/", mustWork = FALSE)
if (!dir.exists(anchors_out_dir)) dir.create(anchors_out_dir, recursive = TRUE, showWarnings = FALSE)
if (any(hif_idx)) writeXStringSet(anchors[hif_idx], file.path(anchors_out_dir, "HIF1a_anchors.fasta") ) else warnmsg("No HIF1a anchors matched")
if (any(phd_idx)) writeXStringSet(anchors[phd_idx], file.path(anchors_out_dir, "PHD2_anchors.fasta")) else warnmsg("No PHD2 anchors matched")
if (any(vhl_idx)) writeXStringSet(anchors[vhl_idx], file.path(anchors_out_dir, "VHL_anchors.fasta")) else warnmsg("No VHL anchors matched")

families <- c("HIF1a","PHD2","VHL")
for (fam in families) {
  anc <- file.path(outdir, "anchors", paste0(fam, "_anchors.fasta"))
  aln <- file.path(outdir, "alignments", paste0(fam, "_aligned.fasta"))
  hmm <- file.path(outdir, "hmms", paste0(fam, ".hmm"))
  if (!file.exists(anc) || file.info(anc)$size == 0) { warnmsg("%s: no anchors; skipping HMM build", fam); next }
  anc_dt <- read_fasta_to_dt(anc)
  if (nrow(anc_dt) < 2) { warnmsg("%s: fewer than 2 anchor sequences; skipping HMM build", fam); next }
  msg("Aligning %s anchors with MAFFT...", fam)
  run("mafft", c("--quiet","--thread", as.character(threads),"--maxiterate","1000","--localpair", anc), to_file = aln, stop_on_fail = TRUE)
  if (!file.exists(aln) || file.info(aln)$size == 0) { warnmsg("%s: alignment missing; skipping hmmbuild", fam); next }
  msg("Building HMM for %s with hmmbuild...", fam)
  run("hmmbuild", c(hmm, aln), stop_on_fail = TRUE)
  if (!(file.exists(hmm) && file.info(hmm)$size > 0)) warnmsg("%s: hmmbuild produced no HMM", fam)
}

for (fam in families) {
  hmm <- file.path(outdir, "hmms", paste0(fam, ".hmm"))
  if (!file.exists(hmm) || file.info(hmm)$size == 0) { warnmsg("%s: missing HMM; skipping hmmsearch", fam); next }
  fam_dir <- file.path(outdir, "hmmsearch", fam); if (!dir.exists(fam_dir)) dir.create(fam_dir, recursive = TRUE, showWarnings = FALSE)
  incE_fam <- if (!is.null(family_params[[fam]]$incE)) family_params[[fam]]$incE else incE_default
  for (prot in proteomes_paths) {
    if (!file.exists(prot)) { warnmsg("Proteome missing, skipping: %s", prot); next }
    sp <- basename(prot)
    tbl   <- file.path(fam_dir, paste0(sp, ".tblout"))
    domtb <- file.path(fam_dir, paste0(sp, ".domtblout"))
    args <- c("--incE", format(incE_fam, scientific = TRUE), "--cpu", as.character(threads),
              "--tblout", tbl, "--domtblout", domtb, hmm, prot)
    run("hmmsearch", args)
  }
}

candidate_registry <- list()
for (fam in families) {
  fam_dir <- file.path(outdir, "hmmsearch", fam)
  cand_dir <- file.path(outdir, "candidates", fam); if (!dir.exists(cand_dir)) dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)
  hits <- list()
  for (prot in proteomes_paths) {
    sp <- basename(prot)
    tbl <- file.path(fam_dir, paste0(sp, ".tblout"))
    dom <- file.path(fam_dir, paste0(sp, ".domtblout"))
    tbl_df <- read_tblout(tbl)
    dom_df <- read_domtbl_simple(dom)
    ids <- unique(na.omit(c(if (!is.null(tbl_df) && nrow(tbl_df)) tbl_df$target else character(0),
                            if (!is.null(dom_df) && nrow(dom_df)) dom_df$target else character(0))))
    ids <- ids[nzchar(ids)]
    outfa <- file.path(cand_dir, paste0(sp, ".candidates.fasta"))
    if (length(ids) && file.exists(prot)) {
      n <- extract_by_ids(prot, ids, outfa)
      hits[[sp]] <- data.frame(species = sp, family = fam, n_hits = n, stringsAsFactors = FALSE)
    } else {
      if (file.exists(outfa)) file.remove(outfa)
      hits[[sp]] <- data.frame(species = sp, family = fam, n_hits = 0L, stringsAsFactors = FALSE)
    }
  }
  candidate_registry[[fam]] <- rbindlist(hits, fill = TRUE)
}
hmm_summary <- rbindlist(candidate_registry, fill = TRUE)
fwrite(hmm_summary, file.path(outdir, "tables", "HMM_hits_summary.csv"))
msg("Wrote HMM hits summary to tables/HMM_hits_summary.csv")

pfam_wanted <- list(HIF1a = c("HLH","PAS"), PHD2 = c("2OG-Fe(II)Oxy"), VHL = c("VHL"))
pfam_ids <- list(HLH="PF00010", PAS="PF00989", `2OG-Fe(II)Oxy`="PF13640", VHL="PF01847")
pfam_exclude <- c("PF02373","PF02375","PF00096","PF00567","PF01388","PF13538","PF18781","PF02008","PF07645","PF00008")

is_potential_PHD2_one <- function(accs, seq_len, min_len=300L, max_len=550L) {
  accs <- if (length(accs)) accs else character(0)
  has_core <- "PF13532" %in% accs
  has_bad  <- any(accs %in% pfam_exclude)
  size_ok  <- is.na(seq_len) || (seq_len >= min_len && seq_len <= max_len)
  list(keep = (has_core && !has_bad && size_ok),
       reason = if (!has_core) "no_PF13532" else if (has_bad) paste0("excluded_domain:", intersect(accs, pfam_exclude)[1]) else if (!size_ok) "size_out_of_range" else "ok")
}

pfam_annot_results <- list()
for (fam in families) {
  cand_dir <- file.path(outdir, "candidates", fam)
  pf_dir <- file.path(outdir, "pfam", fam); if (!dir.exists(pf_dir)) dir.create(pf_dir, recursive = TRUE, showWarnings = FALSE)
  cand_fastas <- list.files(cand_dir, pattern = "\\.fasta$", full.names = TRUE)
  if (!length(cand_fastas)) { pfam_annot_results[[fam]] <- list(); next }
  for (fa in cand_fastas) {
    sp <- gsub("\\.candidates\\.fasta$", "", basename(fa))
    domtbl <- file.path(pf_dir, paste0(sp, ".pfam.domtblout"))
    # run hmmscan of Pfam DB against the candidate FASTA
    run("hmmscan", c("--cpu", as.character(threads), "--domtblout", domtbl, pfam_db, fa))
  }
  entries <- list()
  for (fa in cand_fastas) {
    sp <- gsub("\\.candidates\\.fasta$", "", basename(fa))
    domtbl <- file.path(pf_dir, paste0(sp, ".pfam.domtblout"))
    dt <- if (file.exists(domtbl)) read_domtbl_simple(domtbl) else NULL
    # **Apply paper thresholds here**: keep only domtbl rows with i-Evalue <= 1e-5 AND bitscore >= 10
    dt_filtered <- filter_domtbl_by_threshold(dt, evalue_thresh = 1e-5, bitscore_thresh = 10)
    fa_dt <- read_fasta_to_dt(fa)
    motifs <- lapply(fa_dt$seq, odd_scan)
    dt_motif <- data.frame(target = fa_dt$id, ODD_AP_count = vapply(motifs, function(x) x$n, integer(1)), stringsAsFactors = FALSE)
    filter_df <- NULL
    if (identical(fam, "PHD2")) {
      acc_by_target <- if (!is.null(dt_filtered) && nrow(dt_filtered)) split(dt_filtered$pfam_acc, dt_filtered$target) else list()
      seq_len_map <- setNames(nchar(fa_dt$seq), fa_dt$id)
      pd_rows <- vector("list", length = nrow(fa_dt))
      if (nrow(fa_dt) > 0) {
        for (i in seq_len(nrow(fa_dt))) {
          id <- fa_dt$id[i]
          accs <- acc_by_target[[id]]
          x <- is_potential_PHD2_one(accs = accs, seq_len = seq_len_map[[id]])
          keep_logical <- identical(x$keep, TRUE) || (is.character(x$keep) && tolower(x$keep) %in% c("true", "t", "1")) || (is.numeric(x$keep) && x$keep == 1)
          pd_rows[[i]] <- data.frame(
            target = id,
            is_PHD2_like = if (keep_logical) "Yes" else "No",
            PHD2_filter_reason = as.character(x$reason),
            stringsAsFactors = FALSE
          )
        }
      }
      filter_df <- if (length(pd_rows)) do.call(rbind, pd_rows) else data.frame(target=character(0), is_PHD2_like=character(0), PHD2_filter_reason=character(0), stringsAsFactors = FALSE)
    }
    entries[[sp]] <- list(pfam = if (is.null(dt_filtered)) data.frame(target=character(0), pfam_acc=character(0), pfam_name=character(0), dom_i_evalue=double(0), bitscore=double(0)) else dt_filtered,
                          motifs = dt_motif, filter = filter_df)
  }
  pfam_annot_results[[fam]] <- entries
}
msg("Pfam annotation done")

pfam_wanted <- pfam_wanted
pfam_ids <- pfam_ids

domain_present <- function(df, fam) {
  if (is.null(df) || !nrow(df)) return(FALSE)
  wanted <- unlist(pfam_wanted[[fam]])
  acc_map <- unlist(pfam_ids[wanted], use.names = FALSE)
  any(df$pfam_acc %in% acc_map, na.rm = TRUE)
}

species_names <- basename(proteomes_paths)
summary_rows <- list()
for (sp in species_names) {
  row <- list(species = sp)
  for (fam in families) {
    tblfile <- file.path(outdir, "hmmsearch", fam, paste0(sp, ".tblout"))
    has_hit <- FALSE
    if (file.exists(tblfile)) {
      tdf <- read_tblout(tblfile)
      has_hit <- !is.null(tdf) && nrow(tdf) > 0
    }
    pf_entries <- if (!is.null(pfam_annot_results[[fam]])) pfam_annot_results[[fam]][[sp]] else NULL
    pf_df <- if (!is.null(pf_entries)) pf_entries$pfam else NULL
    pf_ok <- domain_present(pf_df, fam)
    motif_count <- NA_integer_
    if (!is.null(pf_entries) && "motifs" %in% names(pf_entries)) {
      mdf <- pf_entries$motifs
      motif_count <- if (!is.null(mdf) && nrow(mdf) > 0) max(mdf$ODD_AP_count, na.rm = TRUE) else 0L
    }
    row[[paste0(fam, "_HMMhit")]] <- if (has_hit) "Yes" else "No"
    row[[paste0(fam, "_PfamOK")]]  <- if (pf_ok) "Yes" else "No"
    if (fam == "HIF1a") row[["ODD_AP_count"]] <- motif_count
  }
  summary_rows[[sp]] <- as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}
summary_dt <- rbindlist(summary_rows, fill = TRUE)
out_csv <- file.path(outdir, "tables", "Orthologs_Table.csv")
fwrite(summary_dt, out_csv)
msg("[OK] Wrote orthologs table to %s", out_csv)

msg("Pipeline finished. Inspect outputs under: %s", outdir)
msg("Key outputs: anchors/*.fasta, hmms/*.hmm, HMM_hits_summary.csv, Orthologs_Table.csv, candidates/*/*.candidates.fasta, pfam/*/*.pfam.domtblout, logs/*")
cat("\n--- Quick checklist ---\n")
cat("1) Pfam DB (Pfam-A.hmm + .h3*) at:", pfam_db, "\n")
cat("2) HMM hits summary:", file.path(outdir, "tables", "HMM_hits_summary.csv"), "\n")
cat("3) Orthologs table:", out_csv, "\n")
cat("4) Candidate FASTAs (one file per family/species):", file.path(outdir, "candidates"), "\n")