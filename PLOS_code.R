# For PLOS

# co-occurence plot 
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(forcats)
  library(ggplot2); library(patchwork); library(scales); library(cowplot)
})

# ===================== 0) Input & harmonization =========================
df <- dat_withrescued %>% mutate(across(everything(), as.character))

# Normalize SL name and coalesce SL sources if needed
if (!"SL" %in% names(df) && "sl" %in% names(df)) df <- dplyr::rename(df, SL = sl)
if ("sl_2" %in% names(df)) df <- df %>% mutate(SL = dplyr::coalesce(SL, sl_2))
df <- df %>% mutate(SL = str_squish(SL))

# Optional: enforce a consistent site order in legends/bars
df <- df %>%
  mutate(
    site = factor(site, levels = c(site_order, setdiff(sort(unique(site)), site_order)))
  )

# ===================== 1) Column discovery ==============================
has <- function(x) names(df)[str_detect(names(df), regex(x, ignore_case = TRUE))]

ESBL_cols <- unique(c(
  has("^bla_esbl_acquired$"),
  has("^bla_esbl_inh_r_acquired$"),
  has("^rescue_.*_bla_esbl_acquired$"),
  has("^rescue_.*_bla_esbl_inh_r_acquired$"),
  has("^bla_acquired$"),
  has("^rescue_esbl_list$")
))

CP_cols <- unique(c(
  has("^bla_carb_acquired$"),
  has("^rescue_.*_bla_carb_acquired$"),
  has("^rescue_carb_family$")
))

# ---- AmpC cols block you requested ----
AmpC_cols <- unique(c(
  has("^bla_acquired$"),
  has("^beta_lactam_genes_final$"),
  has("^rescue_beta_lactam_genes$"),
  has("^rescue_all_genes_genes$"),
  has("^klebsiella_pneumo_complex__amr__bla_acquired$"),
  has("^rescue_.*_bla_acquired$"),
  has("^ampc_plasmid$")
))

ESBL_cols <- intersect(ESBL_cols, names(df))
CP_cols   <- intersect(CP_cols,   names(df))
AmpC_cols <- intersect(AmpC_cols, names(df))

message("ESBL_cols used: ", paste(ESBL_cols, collapse = ", "))
message("CP_cols used: ",   paste(CP_cols,   collapse = ", "))
message("AmpC_cols used: ", paste(AmpC_cols, collapse = ", "))

# ===================== 2) Tokenizer for multi-valued gene cells =========
split_pat <- "[,;/|\\s]+"  # commas / semicolons / slashes / pipes / whitespace

explode_cols <- function(dat, cols) {
  if (length(cols) == 0) {
    return(tibble(row_id=integer(), site=character(), SL=character(), col=character(), val=character(), tok=character()))
  }
  
  dat %>%
    transmute(row_id = row_number(), site, SL, across(all_of(cols))) %>%
    pivot_longer(-c(row_id, site, SL), names_to = "col", values_to = "val") %>%
    filter(!is.na(val), str_squish(val) != "", val != "-") %>%
    separate_rows(val, sep = split_pat) %>%
    mutate(
      tok = toupper(str_squish(val)) %>%
        # drop version marks and adornments (handles .v1, .V2, trailing markers)
        gsub("\\.V\\d+.*$", "", ., perl = TRUE) %>%
        gsub("\\.v\\d+.*$", "", ., perl = TRUE) %>%
        gsub("[\\^\\*\\?%]+$", "", ., perl = TRUE) %>%
        gsub("\\(.*?\\)", "", ., perl = TRUE) %>%
        gsub("\\[.*?\\]", "", ., perl = TRUE) %>%
        gsub("[\"“”]+", "", ., perl = TRUE) %>%
        gsub("[^A-Z0-9\\-_/]", "", ., perl = TRUE)
    ) %>%
    filter(tok != "")
}

# ===================== 3) Label mappers =================================
to_esbl_label <- function(tok){
  case_when(
    str_detect(tok, "^CTX[-_ ]?M[-_ ]?\\d+$") ~ str_replace(tok, "^CTX[-_ ]?M[-_ ]?", "CTX-M-"),
    str_detect(tok, "^SHV[-_ ]?\\d+$")        ~ str_replace(tok, "^SHV[-_ ]?",      "SHV-"),
    str_detect(tok, "^TEM[-_ ]?\\d+$")        ~ str_replace(tok, "^TEM[-_ ]?",      "TEM-"),
    str_detect(tok, "^PER[-_ ]?\\d+$")        ~ str_replace(tok, "^PER[-_ ]?",      "PER-"),
    str_detect(tok, "^VEB[-_ ]?\\d+$")        ~ str_replace(tok, "^VEB[-_ ]?",      "VEB-"),
    str_detect(tok, "^GES[-_ ]?\\d+$")        ~ str_replace(tok, "^GES[-_ ]?",      "GES-"),
    TRUE ~ tok
  )
}
is_esbl_family <- function(tok) {
  str_detect(tok, regex("^(CTX[-_ ]?M|SHV|TEM|PER|VEB|GES)\\b", ignore_case = TRUE))
}

# ===================== 4) Mutation sourcing + normalization (SHV/TEM) ========
take_mut_col <- function(df, candidates) {
  avail <- candidates[candidates %in% names(df)]
  if (!length(avail)) return(rep(NA_character_, nrow(df)))
  for (nm in avail) {
    v <- df[[nm]]
    if (any(!is.na(v) & trimws(v) != "")) return(v)
  }
  df[[avail[1]]]
}

aa3_to_1 <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E", GLN="Q", GLY="G",
  HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P",
  SER="S", THR="T", TRP="W", TYR="Y", VAL="V"
)
collapse_triplets <- function(x) {
  if (all(is.na(x))) return(x)
  y <- toupper(x)
  for (k in names(aa3_to_1)) {
    y <- gsub(paste0("\\b", k, "\\b"), aa3_to_1[[k]], y)
    y <- gsub(k, aa3_to_1[[k]], y, fixed = TRUE)
  }
  y
}
norm_mut_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- collapse_triplets(x)
  str_squish(toupper(x))
}
mut_to <- function(txt, pos, aa){
  x <- norm_mut_text(txt)
  pat <- paste0("\\b", pos, "\\s*[-_:/]?\\s*", aa, "\\b|\\b[A-Z]{1,3}\\s*", pos, "\\s*[-_:/]?\\s*", aa, "\\b")
  str_detect(x, regex(pat, ignore_case = TRUE))
}
mut_any <- function(txt, pos){
  x <- norm_mut_text(txt)
  pat <- paste0("\\b", pos, "\\s*[A-Z]\\b|\\b[A-Z]{1,3}\\s*", pos, "\\s*[A-Z]\\b")
  str_detect(x, regex(pat, ignore_case = TRUE))
}

shv_txt_raw <- take_mut_col(df, c("shv_mut_text","shv_mutations","shvmut_text","shvmut_base","shv_mut","shv_mutation_text"))
tem_txt_raw <- take_mut_col(df, c("tem_mut_text","tem_mutations","temmut_text","tem_mut","tem_mutation_text"))
if (all(is.na(tem_txt_raw) | trimws(tem_txt_raw) == "")) tem_txt_raw <- shv_txt_raw

# SHV by mutation rule
shv_238X <- mut_any(shv_txt_raw, 238)
shv_179X <- mut_any(shv_txt_raw, 179)
shv_169X <- mut_any(shv_txt_raw, 169)
shv_148X <- mut_any(shv_txt_raw, 148)
shv_240K <- mut_to(shv_txt_raw, 240, "K")
shv_35Q  <- mut_to(shv_txt_raw,  35, "Q")

df$ESBL_SHV_by_mut <- (!is.na(shv_txt_raw) & trimws(shv_txt_raw) != "") &
  (shv_238X | shv_179X | shv_169X | shv_148X | (shv_240K & shv_35Q))

# TEM by mutation rule
df$ESBL_TEM_by_mut <- (!is.na(tem_txt_raw) & trimws(tem_txt_raw) != "") & (
  mut_to(tem_txt_raw, 238, "S") |
    mut_to(tem_txt_raw, 164, "S") |
    mut_to(tem_txt_raw, 104, "K") |
    mut_to(tem_txt_raw, 240, "K") |
    mut_to(tem_txt_raw, 237, "T") |
    mut_to(tem_txt_raw, 182, "T")
)

# ===================== 5) Build long data frames (ESBL, AmpC, CP) ============
# ---- ESBL tokens (drop SHV-12 by name; drop TEM-by-name; keep TEM via mutations only) ----
esbl_long <- explode_cols(df, ESBL_cols) %>%
  filter(is_esbl_family(tok)) %>%
  mutate(Group = to_esbl_label(tok)) %>%
  filter(Group != "SHV-12", !str_detect(Group, "^TEM-")) %>%
  distinct(row_id, site, SL, Group)

mut_long <- df %>%
  mutate(row_id = row_number()) %>%
  transmute(
    row_id, site, SL,
    `SHV (mutation)` = ESBL_SHV_by_mut %in% TRUE,
    `TEM (mutation)` = ESBL_TEM_by_mut %in% TRUE
  ) %>%
  pivot_longer(-c(row_id, site, SL), names_to = "Group", values_to = "flag") %>%
  filter(flag) %>%
  select(-flag)

esbl_long <- bind_rows(esbl_long, mut_long) %>% distinct(row_id, site, SL, Group)

esbl_dat <- bind_rows(
  esbl_long,
  df %>% mutate(row_id = row_number()) %>%
    filter(!(row_id %in% esbl_long$row_id)) %>%
    transmute(row_id, site, SL, Group = "None")
) %>%
  filter(!is.na(SL), SL != "") %>%
  distinct(row_id, site, SL, Group)

write_xlsx(esbl_long, "esbl_long_cooccurance.xlsx")

# ---- AmpC (DHA/CMY/LAP) WITHOUT Mixed; co-occurrence allowed ----
ampc_long <- explode_cols(df, AmpC_cols) %>%
  mutate(tok2 = gsub("^BLA[-_]*", "", tok, perl = TRUE)) %>%
  mutate(
    Group = case_when(
      # anchored so LAPJ01000014 is NOT called LAP
      str_detect(tok2, "^DHA([-_ ]?\\d+)?$") ~ "DHA",
      str_detect(tok2, "^CMY([-_ ]?\\d+)?$") ~ "CMY",
      str_detect(tok2, "^LAP([-_ ]?\\d+)?$") ~ "LAP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  distinct(row_id, site, SL, Group)

ampc_dat <- bind_rows(
  ampc_long,
  df %>% mutate(row_id = row_number()) %>%
    filter(!(row_id %in% ampc_long$row_id)) %>%
    transmute(row_id, site, SL, Group = "None")
) %>%
  filter(!is.na(SL), SL != "") %>%
  distinct(row_id, site, SL, Group)

write_xlsx(ampc_long, "ampc_long_co_occurance.xlsx")

# ---- Carbapenemase (token-based) ----
cp_long <- explode_cols(df, CP_cols) %>%
  mutate(tok2 = gsub("^BLA[-_]*", "", tok, perl = TRUE)) %>%
  mutate(
    Group = case_when(
      str_detect(tok2, "^NDM[-_ ]?1\\b")       ~ "NDM-1",
      str_detect(tok2, "^NDM[-_ ]?5\\b")       ~ "NDM-5",
      str_detect(tok2, "^NDM(\\b|[-_ ]?\\d+)") ~ "NDM (other)",
      
      str_detect(tok2, "^OXA[-_ ]?204\\b")     ~ "OXA-204",
      str_detect(tok2, "^OXA[-_ ]?48\\b")      ~ "OXA-48",
      str_detect(tok2, "^OXA(\\b|[-_ ]?\\d+)") ~ "OXA (other)",
      
      str_detect(tok2, "^KPC(\\b|[-_ ]?\\d+)") ~ "KPC",
      str_detect(tok2, "^VIM(\\b|[-_ ]?\\d+)") ~ "VIM",
      str_detect(tok2, "^IMP(\\b|[-_ ]?\\d+)") ~ "IMP",
      str_detect(tok2, "^GES(\\b|[-_ ]?\\d+)") ~ "GES",
      
      TRUE ~ "Other/Unclear"
    )
  ) %>%
  distinct(row_id, site, SL, Group)

cp_dat <- bind_rows(
  cp_long,
  df %>% mutate(row_id = row_number()) %>%
    filter(!(row_id %in% cp_long$row_id)) %>%
    transmute(row_id, site, SL, Group = "None")
) %>%
  filter(!is.na(SL), SL != "") %>%
  distinct(row_id, site, SL, Group)

write_xlsx(cp_long, "cp_long_co_occurance.xlsx")

# ===================== 6) Palette & plotter ===============================
plos_no_black <- c("#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

make_named_site_palette <- function(dat) {
  sites <- levels(factor(dat$site))
  pal <- rep(plos_no_black, length.out = length(sites))
  names(pal) <- sites
  pal
}

plot_dot_bars <- function(dat, y_title, main_title, min_n_SL = 2, group_order = NULL){
  
  # filter SLs by UNIQUE isolates (not by determinant rows)
  if (min_n_SL > 0) {
    keep <- dat %>%
      distinct(row_id, SL) %>%
      count(SL, name = "n_iso") %>%
      filter(n_iso >= min_n_SL) %>%
      pull(SL)
    dat <- dat %>% filter(SL %in% keep)
  }
  
  # group ordering
  if (!is.null(group_order)) {
    dat <- dat %>% mutate(Group = factor(Group, levels = group_order))
    group_levels <- group_order
  } else {
    dat <- dat %>% mutate(Group = fct_infreq(factor(Group)))
    dat <- dat %>% mutate(Group = fct_relevel(Group, "None", after = Inf))
    group_levels <- levels(dat$Group)
  }
  
  dat <- dat %>%
    mutate(
      SL   = fct_infreq(factor(SL)),
      site = factor(site)
    )
  
  pal <- make_named_site_palette(dat)
  
  # Center (jittered points)
  p_center <- ggplot(dat, aes(x = SL, y = Group, color = site)) +
    geom_jitter(width = 0.25, height = 0.25, alpha = 0.85, size = 1.9) +
    scale_color_manual(values = pal, name = "Hospital") +
    scale_y_discrete(limits = group_levels, drop = FALSE) +
    labs(title = main_title, x = "Sequence lineage (SL)", y = y_title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  # Top (counts by SL × site) using UNIQUE isolates
  top_df <- dat %>% distinct(row_id, SL, site) %>% count(SL, site, name = "n")
  p_top <- ggplot(top_df, aes(x = SL, y = n, fill = site)) +
    geom_col(position = "stack") +
    geom_text(aes(label = n),
              position = position_stack(vjust = 0.5), size = 2.6, color = "white") +
    scale_fill_manual(values = pal) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "none"
    )
  
  # Right (counts by Group × site) using UNIQUE isolates per group
  side_df <- dat %>% distinct(row_id, Group, site) %>%
    count(Group, site, name = "n") %>%
    mutate(Group = factor(Group, levels = group_levels))
  p_side <- ggplot(side_df, aes(x = n, y = Group, fill = site)) +
    geom_col(position = "stack") +
    geom_text(aes(label = ifelse(n >= 2, n, "")),
              position = position_stack(vjust = 0.5), size = 2.6, color = "white") +
    scale_fill_manual(values = pal) +
    scale_y_discrete(limits = group_levels, drop = FALSE) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "none"
    )
  
  (p_top + patchwork::plot_spacer()) / (p_center + p_side) +
    patchwork::plot_layout(heights = c(1, 4), widths = c(4, 1))
}

# ===================== 7) Build plots =====================================
p_ESBL <- plot_dot_bars(
  esbl_dat,
  y_title    = "ESBL determinant",
  main_title = "ESBL determinants by SL and hospital",
  min_n_SL   = 2
)

p_AmpC <- plot_dot_bars(
  ampc_dat,
  y_title     = "AmpC determinant",
  main_title  = "AmpC determinants (DHA/CMY/LAP; no mixed) by SL and hospital",
  min_n_SL    = 1,  # <-- key: show singleton SLs so LAP doesn't disappear
  group_order = c("DHA","CMY","LAP","None")
)

p_CP <- plot_dot_bars(
  cp_dat,
  y_title    = "Carbapenemase determinant",
  main_title = "Carbapenemase determinants by SL and hospital",
  min_n_SL   = 2
)

fig_plos <- cowplot::plot_grid(
  p_ESBL, p_AmpC, p_CP,
  ncol = 1,
  rel_heights = c(1, 1, 1),
  labels = c("A", "B", "C"),
  label_fontface = "bold",
  label_size = 12,
  label_x = 0.01, label_y = 0.99,
  hjust = 0, vjust = 1
)

ggsave("Fig_SL_ESBL_AmpC_CP_PLOS_triple.pdf", plot = fig_plos, width = 15, height = 18, units = "in")
ggsave("Fig_SL_ESBL_AmpC_CP_PLOS_triple.png", plot = fig_plos, width = 15, height = 18, units = "in", dpi = 300)

# Optional single-panel exports
ggsave("Fig_ESBL_SL_dot_bars_clean.pdf", p_ESBL, width = 20, height = 10)
ggsave("Fig_AmpC_SL_dot_bars_clean.pdf", p_AmpC, width = 20, height = 10)
ggsave("Fig_CP_SL_dot_bars_clean.pdf",   p_CP,   width = 20, height = 10)



## UpSet Plot #####
#----- new AMRgen plot
df_markers <- tibble(
  id      = seq_len(nrow(dat0)),
  pheno   = pheno,          # keeps S/R colouring
  disk    = disk_caz,       # assay values (for ordering + boxplots)
  # (optional bookkeeping columns; not used as rows in upset)
  R       = CAZ_R,
  NWT     = 0L,
  
  # ESBL determinants
  CTX_M_15  = as.integer(has_CTXM_15),
  CTX_M_14  = as.integer(has_CTXM_14),
  CTX_M_156 = as.integer(has_CTXM_156),
  CTX_M_3   = as.integer(has_CTXM_3),
  SHV_mut   = as.integer(SHV_mut_effective),
  TEM_mut   = as.integer(dat0$tem_is_esbl),
  
  # AmpC & Porins
  AmpC_plasmid = as.integer(AmpC_plasmid),
  OmpK36_mut   = as.integer(Porin_OmpK36),
  OmpK35_mut   = as.integer(Porin_OmpK35)
  
  # NOTE: deliberately NOT adding CAZ_R as a marker row
)

# Order for UpSet (NO CAZ_R here)
marker_levels <- c(
  "CTX_M_15","CTX_M_14","CTX_M_156","CTX_M_3",
  "SHV_mut","TEM_mut",
  "AmpC_plasmid","OmpK36_mut","OmpK35_mut"
)

df_upset <- df_markers %>%
  select(id, pheno, disk, R, NWT, all_of(marker_levels)) %>%
  pivot_longer(all_of(marker_levels), names_to = "Marker", values_to = "value") %>%
  filter(value == 1) %>%
  mutate(Marker = factor(Marker, levels = marker_levels)) %>%
  pivot_wider(names_from = Marker, values_from = value, values_fill = 0)

# ---------- Plot ----------
plots <- AMRgen::amr_upset(
  binary_matrix = df_upset,
  min_set_size = 1,
  
  # (3) order columns by assay value
  order = "value",
  
  # (2) add set-size bar + printed counts
  plot_set_size = TRUE,
  print_set_size = TRUE,
  
  # (4) remove the middle stacked category barplot
  plot_category = FALSE,
  
  # keep or disable category counts (irrelevant if plot_category=FALSE)
  print_category_counts = FALSE,
  
  assay = "disk",
  
  # small fix: argument name in AMRgen is boxplot_col (not boxplot_colour)
  boxplot_col = "grey"
)

p_caz <- plots$plot +
  ggtitle("Resistance mechanism combinations (CAZ focus)") +
  labs(caption = "EUCAST: I counted I with S. Determinants: CTX-M (15/14/156/3) by name; SHV/TEM by mutation only (SHV-12 excluded); AmpC plasmid; OmpK35/OmpK36 alterations.")

print(p_caz)

ggsave(
  filename = "CAZ_amr_upset_with_AmpC_Porin.pdf",
  plot = p_caz,
  width = 11, height = 8, units = "in"
)


names(plots)
# show all data-frame-ish objects inside plots
lapply(plots, function(x) if (is.data.frame(x)) head(x, 3) else NULL)

plots$summary %>%
  dplyr::filter(marker_count == 1,
                marker_list %in% c("CTX_M_15", "CTX_M_14")) %>%
  dplyr::select(marker_list, n, median, q25, q75)

plots$summary

write.xlsx(plots$summary, "upset_summary_CAZ.xlsx")


#----- new AMRgen plot CTX
disk_ctx <- suppressWarnings(as.numeric(dat0$cefotaxime_5))  # CTX disk diffusion

df_markers <- tibble(
  id    = seq_len(nrow(dat0)),
  pheno = ctx,   # <-- use CTX phenotype here (or pheno if already CTX)
  disk  = disk_ctx,
  R     = dat0$ctx_r,
  NWT   = 0L,
  
  # ESBL determinants
  CTX_M_15  = as.integer(has_CTXM_15),
  CTX_M_14  = as.integer(has_CTXM_14),
  CTX_M_156 = as.integer(has_CTXM_156),
  CTX_M_3   = as.integer(has_CTXM_3),
  SHV_mut   = as.integer(SHV_mut_effective),
  TEM_mut   = as.integer(dat0$tem_is_esbl),
  
  # AmpC & Porins
  AmpC_plasmid = as.integer(AmpC_plasmid),
  OmpK36_mut   = as.integer(Porin_OmpK36),
  OmpK35_mut   = as.integer(Porin_OmpK35)
)

marker_levels <- c(
  "CTX_M_15","CTX_M_14","CTX_M_156","CTX_M_3",
  "SHV_mut","TEM_mut",
  "AmpC_plasmid","OmpK36_mut","OmpK35_mut"
)

df_upset %>%
  unite("combo", all_of(marker_levels), sep="_", remove=FALSE) %>%
  count(combo, sort=TRUE) %>%
  filter(grepl("_1_", combo))   # combos with at least one "1" somewhere

df_upset <- df_markers %>%
  select(id, pheno, disk, R, NWT, all_of(marker_levels)) %>%
  pivot_longer(all_of(marker_levels), names_to = "Marker", values_to = "value") %>%
  filter(value == 1) %>%
  mutate(Marker = factor(Marker, levels = marker_levels)) %>%
  pivot_wider(names_from = Marker, values_from = value, values_fill = 0)

plots <- AMRgen::amr_upset(
  binary_matrix = df_upset,
  min_set_size = 1,
  order = "value",
  plot_set_size = TRUE,
  print_set_size = TRUE,
  plot_category = FALSE,
  print_category_counts = FALSE,
  assay = "disk",
  boxplot_col = "grey"
)

p_ctx <- plots$plot +
  ggtitle("Resistance mechanism combinations (CTX focus)") +
  labs(caption = "EUCAST: I counted I with S. Determinants: CTX-M (15/14/156/3) by name; SHV/TEM by mutation only (SHV-12 excluded); AmpC plasmid; OmpK35/OmpK36 alterations.")

print(p_ctx)

ggsave(
  filename = "CTX_amr_upset_with_AmpC_Porin.pdf",
  plot = p_ctx,
  width = 11, height = 8, units = "in"
)


names(plots)
# show all data-frame-ish objects inside plots
lapply(plots, function(x) if (is.data.frame(x)) head(x, 3) else NULL)

plots$summary %>%
  dplyr::filter(marker_count == 1,
                marker_list %in% c("CTX_M_15", "CTX_M_14")) %>%
  dplyr::select(marker_list, n, median, q25, q75)

plots$summary

write.xlsx(plots$summary, "upset_summary_CTX.xlsx")

### PPV
# CTX
geno_ctx <- df_markers %>%
  select(id, all_of(marker_levels)) %>%
  pivot_longer(
    cols = all_of(marker_levels),
    names_to = "marker",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  transmute(
    id = id,
    marker = marker,
    drug_class = "Cephalosporins",
    drug_agent = "Cefotaxime"
  )

pheno_ctx <- df_markers %>%
  transmute(
    id = id,
    drug_agent = "Cefotaxime",
    pheno = case_when(
      pheno == "R" ~ "R",
      pheno %in% c("S", "I") ~ "S",
      TRUE ~ NA_character_
    ),
    disk = disk,
    mic = NA_real_
  ) %>%
  filter(!is.na(pheno))

solo_ctx <- AMRgen::solo_ppv_analysis(
  geno_table = geno_ctx,
  pheno_table = pheno_ctx,
  antibiotic = "Cefotaxime",
  drug_class_list = c("Cephalosporins"),
  geno_sample_col = "id",
  pheno_sample_col = "id",
  sir_col = "pheno",
  marker_col = "marker",
  keep_assay_values = FALSE,
  min = 1
)

names(solo_ctx)
solo_ctx$solo_stats
print(solo_ctx$combined_plot)

write.xlsx(solo_ctx$solo_stats, "CTX_solo_ppv_stats.xlsx")

ggsave(
  filename = "CTX_solo_PPV_AMRgen.pdf",
  plot = solo_ctx$combined_plot,
  width = 10, height = 7, units = "in"
)


### -- convergent characteristics--
# Minimal deps
library(dplyr)
library(stringr)
library(tidyr)
library(scales)

# Helper: pretty n/N (%)
nN_fmt <- function(n, N) sprintf("%d/%d (%.1f%%)", n, N, 100 * ifelse(N > 0, n/N, NA_real_))

# Helper: ICU vs Non-ICU from hospital_ward free text
derive_ward_group <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    str_detect(tolower(x), "\\bicu\\b|intensive|critical|\\bccu\\b|\\bhdu\\b|\\bnicu\\b|\\bpicu\\b") ~ "ICU",
    TRUE ~ "Non-ICU"
  )
}

# --- Canonicalize fields in BOTH frames (keep names exactly as you have) ------
df_conv_flag <- df_conv_flag %>%
  mutate(
    sl            = as.character(sl),
    hospital      = as.character(site),
    specimen_type = as.character(specimen_type),
    patient_type  = as.character(patient_type),
    year          = suppressWarnings(as.integer(year)),
    ward_group    = derive_ward_group(as.character(hospital_ward))
  )

conv_only <- conv_only %>%
  mutate(
    sl            = as.character(sl),
    hospital      = as.character(site),
    specimen_type = as.character(specimen_type),
    patient_type  = as.character(patient_type),
    year          = suppressWarnings(as.integer(year)),
    ward_group    = derive_ward_group(as.character(hospital_ward))
  )

# --- 1) n/N (%) by SL ---------------------------------------------------------
iuc_by_sl <- df_conv_flag %>%
  count(sl, name = "N") %>%
  left_join(conv_only %>% count(sl, name = "n"), by = "sl") %>%
  mutate(n = replace_na(n, 0L),
         pct = ifelse(N > 0, n/N, NA_real_),
         `n/N (%)` = nN_fmt(n, N)) %>%
  arrange(desc(pct), desc(n))
print(iuc_by_sl)

# --- 2) n/N (%) by hospital ---------------------------------------------------
iuc_by_hospital <- df_conv_flag %>%
  count(hospital, name = "N") %>%
  left_join(conv_only %>% count(hospital, name = "n"), by = "hospital") %>%
  mutate(n = replace_na(n, 0L),
         pct = ifelse(N > 0, n/N, NA_real_),
         `n/N (%)` = nN_fmt(n, N)) %>%
  arrange(desc(pct), desc(n))
print(iuc_by_hospital)

# --- 3) n/N (%) by specimen_type ---------------------------------------------
iuc_by_specimen <- df_conv_flag %>%
  count(specimen_type, name = "N") %>%
  left_join(conv_only %>% count(specimen_type, name = "n"), by = "specimen_type") %>%
  mutate(n = replace_na(n, 0L),
         pct = ifelse(N > 0, n/N, NA_real_),
         `n/N (%)` = nN_fmt(n, N)) %>%
  arrange(desc(pct), desc(n))
print(iuc_by_specimen)

# --- 4) ICU vs Non-ICU: X% vs Y% + Fisher’s exact test -----------------------
icu_tab <- df_conv_flag %>%
  filter(ward_group %in% c("ICU","Non-ICU")) %>%
  count(ward_group, name = "N_total") %>%
  left_join(conv_only %>% filter(ward_group %in% c("ICU","Non-ICU")) %>%
              count(ward_group, name = "n_conv"),
            by = "ward_group") %>%
  mutate(n_conv = replace_na(n_conv, 0L),
         n_nonconv = N_total - n_conv)

xt_icu <- as.matrix(icu_tab %>%
                      arrange(match(ward_group, c("ICU","Non-ICU"))) %>%
                      select(n_nonconv, n_conv))
rownames(xt_icu) <- icu_tab$ward_group

ft_icu <- fisher.test(xt_icu)

icu_summary <- icu_tab %>%
  mutate(pct = ifelse(N_total > 0, n_conv / N_total, NA_real_),
         `n/N (%)` = nN_fmt(n_conv, N_total)) %>%
  select(ward_group, n_conv, N_total, `n/N (%)`, pct)
print(icu_summary)
cat(sprintf("Fisher (ICU vs Non-ICU): OR=%.3f, 95%% CI [%.3f, %.3f], p=%.4f\n",
            unname(ft_icu$estimate),
            unname(ft_icu$conf.int[1]),
            unname(ft_icu$conf.int[2]),
            ft_icu$p.value))

# --- 5) HAI vs CAI (patient_type): X% vs Y% + Fisher’s exact test ------------
pt_tab <- df_conv_flag %>%
  filter(!is.na(patient_type)) %>%
  count(patient_type, name = "N_total") %>%
  left_join(conv_only %>%
              filter(!is.na(patient_type)) %>%
              count(patient_type, name = "n_conv"),
            by = "patient_type") %>%
  mutate(n_conv = replace_na(n_conv, 0L),
         n_nonconv = N_total - n_conv)

# keep a clean 2x2 if you have HAI/CAI; otherwise Fisher will still run
pt_tab2 <- if (all(c("HAI","CAI") %in% pt_tab$patient_type)) {
  pt_tab %>% filter(patient_type %in% c("HAI","CAI")) %>%
    arrange(match(patient_type, c("HAI","CAI")))
} else pt_tab

xt_pt <- as.matrix(pt_tab2 %>% select(n_nonconv, n_conv))
rownames(xt_pt) <- pt_tab2$patient_type

ft_pt <- fisher.test(xt_pt)

pt_summary <- pt_tab %>%
  mutate(pct = ifelse(N_total > 0, n_conv / N_total, NA_real_),
         `n/N (%)` = nN_fmt(n_conv, N_total)) %>%
  select(patient_type, n_conv, N_total, `n/N (%)`, pct)
print(pt_summary)
cat(sprintf("Fisher (patient_type): OR=%.3f, 95%% CI [%.3f, %.3f], p=%.4f\n",
            unname(ft_pt$estimate),
            unname(ft_pt$conf.int[1]),
            unname(ft_pt$conf.int[2]),
            ft_pt$p.value))

# --- 6) iuc prevalence by year (n/N %) ---------------------------------------
# 6) iuc prevalence by year + Cochran–Armitage test
iuc_by_year <- df_conv_flag %>%
  filter(!is.na(year)) %>%
  count(year, name = "N_total") %>%
  left_join(conv_only %>% filter(!is.na(year)) %>% count(year, name = "n_conv"),
            by = "year") %>%
  mutate(n_conv = replace_na(n_conv, 0L),
         pct = ifelse(N_total > 0, n_conv / N_total, NA_real_),
         `n/N (%)` = nN_fmt(n_conv, N_total)) %>%
  arrange(year)

# Cochran–Armitage test for trend
xt_year <- as.matrix(iuc_by_year %>% select(n_conv, N_total))
trend_p <- NA_real_
if (nrow(iuc_by_year) >= 3) {
  trend_test <- DescTools::CochranArmitageTest(
    matrix(c(iuc_by_year$n_conv,
             iuc_by_year$N_total - iuc_by_year$n_conv), ncol = 2))
  trend_p <- trend_test$p.value
}
iuc_by_year$Trend_p <- ifelse(!is.na(trend_p),
                              sprintf("Cochran–Armitage p=%.4f", trend_p),
                              NA_character_)

# --- Save everything ---
write_xlsx(
  list(
    "iuc_by_SL"         = iuc_by_sl,
    "iuc_by_hospital"   = iuc_by_hosp,
    "iuc_by_specimen"   = iuc_by_specimen,
    "ICU_vs_NonICU"     = icu_summary,
    "PatientType_HAI_CAI" = pt_summary,
    "iuc_by_year"       = iuc_by_year
  ),
  path = "iuc_convergence_summary.xlsx"
)

cat("\nAll summaries written to 'iuc_convergence_summary.xlsx'\n")


### ============================================================
### Mechanism classification (CTX-M + SHV mutations + TEM status)
### AmpC is ALIGNED to your "AmpC phenotype check (code that works)"
###   - Token-safe DHA/CMY/LAP matching
###   - Option to use bla_acquired only (matches manual counts, e.g. DHA=14)
###   - Exports ampc_dha/cmy/lap for audit
### Carb fix:
###   - Includes OXA-204 in carbapenemase regex
### Outputs categories exactly:
###   CTX-M-only;
###   CTX-M + SHV mutation ± TEM mutation;
###   porin mutation only;
###   AmpC beta-lactamase;
###   carbapenemase + porin mutation;
###   carbapenemase + porin mutation + AmpC;
###   carbapenemase only;
###   none.
### ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
  library(janitor)
  library(writexl)
  library(stringi)
})

# =========================
# USER SWITCH (IMPORTANT)
# =========================
# TRUE  -> AmpC calls (esp DHA) match what you manually counted in bla_acquired (e.g. DHA=14)
# FALSE -> AmpC calls from all r4 sources (can yield DHA=15 if another source contains DHA)
use_bla_acquired_only_for_ampc <- TRUE


# ---------------- helpers ----------------
core_id <- function(x){
  s <- as.character(x)
  s <- gsub("^PID-\\d+-", "", s, ignore.case = TRUE)
  s <- gsub("_S\\d+.*$", "", s, ignore.case = TRUE)
  trimws(s)
}

coalesce_cols_chr <- function(df, candidates){
  cols <- intersect(candidates, names(df))
  if (!length(cols)) return(rep(NA_character_, nrow(df)))
  out <- as.character(df[[cols[1]]])
  if (length(cols) > 1) {
    for (nm in cols[-1]) out <- dplyr::coalesce(out, as.character(df[[nm]]))
  }
  out
}

# robust "flag" coercion for 1/0, TRUE/FALSE, yes/no, etc.
as_flag <- function(x){
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  z <- tolower(trimws(as.character(x)))
  z[is.na(z)] <- ""
  z %in% c("1","true","t","yes","y","pos","positive","esbl")
}

# Unicode-robust normaliser for gene strings
norm_ctx_text <- function(v){
  x <- as.character(v); x[is.na(x)] <- ""
  x <- stringi::stri_trans_nfkc(x)
  x <- stringi::stri_replace_all_regex(x, "[\\p{Cf}\\p{Cc}]", "")
  x <- stringi::stri_replace_all_regex(x, "\\p{Z}+", " ")
  x <- stringi::stri_replace_all_regex(x, "[\\p{Pd}]", "-")
  x <- chartr("\u039C\u00B5\u0421\u0422\u0425\u041C", "MMCTXM", x)
  x <- stringi::stri_replace_all_regex(x, "[*^?]", "")
  toupper(x)
}

# Extract clean CTX-M tokens (e.g. "CTX-M-103;CTX-M-114"), empty "" if none
extract_ctxm_tokens <- function(v){
  x <- norm_ctx_text(v)
  hits <- stringi::stri_extract_all_regex(x, "(?iu)CTX[-_ ]?M[-_ ]?\\d+")
  vapply(hits, function(h){
    h <- unique(h[!is.na(h)])
    if (length(h)) paste(h, collapse = ";") else ""
  }, character(1))
}

# CTX-M presence detector (NA-safe & Unicode-robust)
has_ctxm <- function(v){
  x <- as.character(v); x[is.na(x)] <- ""
  x <- stringi::stri_trans_nfkc(x)
  x <- stringi::stri_replace_all_regex(x, "[\\p{Cf}\\p{Cc}]", "")
  x <- stringi::stri_replace_all_regex(x, "\\p{Z}+", " ")
  x <- stringi::stri_replace_all_regex(x, "[\\p{Pd}]", "-")
  x <- chartr("\u039C\u00B5\u0421\u0422\u0425\u041C", "MMCTXM", x)
  x <- stringi::stri_replace_all_regex(x, "[*^?]", "")
  stringi::stri_detect_regex(x, "(?iu)CTX\\s*[-_]?\\s*M")
}

# Paste + normalise multiple columns into one text field (for regex detection)
paste_cols_norm <- function(d, cols){
  cols <- intersect(cols, names(d))
  if (!length(cols)) return(rep(NA_character_, nrow(d)))
  m <- as.data.frame(lapply(d[cols], norm_ctx_text))
  txt <- do.call(paste, c(m, list(sep=";")))
  txt[txt %in% c("", "NA")] <- NA
  txt
}

# collapse free-text mutation columns into one uppercase string (audit only)
collapse_mut_cols <- function(d, cols){
  cols <- intersect(cols, names(d))
  if (!length(cols)) return(rep(NA_character_, nrow(d)))
  mat <- as.data.frame(lapply(d[cols], as.character))
  out <- apply(mat, 1, function(row){
    row <- row[!is.na(row) & nzchar(row)]
    if (!length(row)) return(NA_character_)
    out <- paste(row, collapse = ";")
    out <- stringi::stri_trans_nfkc(out)
    str_squish(toupper(out))
  })
  out
}

# tolerant matchers for SHV mutation tokens (e.g. "E240K" or "240K")
norm_mut_text <- function(x) toupper(str_squish(ifelse(is.na(x), "", x)))
mut_to <- function(txt, pos, aa){
  z <- norm_mut_text(txt)
  pat <- paste0("\\b", pos, "\\s*[-_:/]?\\s*", aa, "\\b|\\b[A-Z]{1,3}\\s*", pos, "\\s*[-_:/]?\\s*", aa, "\\b")
  str_detect(z, regex(pat))
}
mut_any <- function(txt, pos){
  z <- norm_mut_text(txt)
  pat <- paste0("\\b", pos, "\\s*[A-Z]\\b|\\b[A-Z]{1,3}\\s*", pos, "\\s*[A-Z]\\b")
  str_detect(z, regex(pat))
}

# ---- token-safe gene hit (used for AmpC alignment) ----
hit <- function(txt, pattern){
  !is.na(txt) & str_detect(txt, regex(pattern, ignore_case = TRUE))
}
tok <- function(gene){
  # matches DHA, DHA-1, blaDHA-1, BLADHA-1, etc (prevents substring junk)
  paste0("(^|[^A-Z0-9])(?:BLA)?", gene, "(?:\\s*[-_ ]\\s*\\d+)?($|[^A-Z0-9])")
}


# ---------------- load / clean ----------------
infile <- "/Users/danait/Dropbox/WGS_analysis/WGS_data_analysed/dat_withrescued_22.xlsx"
dat_withrescued <- read_xlsx(infile, guess_max = 200000) %>% clean_names()

# ---------------- unified ID + dedup BEFORE rules ----------------
dat0 <- dat_withrescued %>%
  mutate(
    id_unified      = coalesce_cols_chr(., c("id_export","id","genome_id","sample","strain","isolate","isolate_id")),
    id_core_unified = core_id(id_unified)
  )

has_rescue_col <- "rescue_all_genes_genes" %in% names(dat0)

dat1 <- dat0 %>%
  arrange(
    desc(if (has_rescue_col) !is.na(.data[["rescue_all_genes_genes"]]) else FALSE),
    id_core_unified, id_unified
  ) %>%
  distinct(id_core_unified, .keep_all = TRUE)

message("N isolates (deduped) = ", nrow(dat1))

# ---------------- rule sources ----------------
r1 <- c("bla_esbl_acquired","truncated_resistance_hits","spurious_resistance_hits",
        "rescue_esbl_list","rescue_all_genes_genes")   # ESBL = CTX-M evidence sources

r2 <- c("bla_carb_acquired","rescue_carb_list","rescue_carb_family",
        "rescue_klebsiella_pneumo_complex_amr_bla_carb_acquired")

r3 <- c("omp_mutations","rescue_klebsiella_pneumo_complex_amr_omp_mutations")

r4 <- c("bla_acquired","klebsiella_pneumo_complex__amr__bla_acquired","beta_lactam_genes_final",
        "rescue_beta_lactam_genes","rescue_all_genes_genes",
        "rescue_klebsiella_pneumo_complex_amr_bla_acquired")

# ---------------- CTX-M detection (robust + audit tokens) ----------------
dat1 <- dat1 %>%
  mutate(
    ctxm_trunc_tok_raw = if ("truncated_resistance_hits" %in% names(.)) extract_ctxm_tokens(.data[["truncated_resistance_hits"]]) else "",
    ctxm_spur_tok_raw  = if ("spurious_resistance_hits"  %in% names(.)) extract_ctxm_tokens(.data[["spurious_resistance_hits"]])  else "",
    
    hit_bla_esbl_acquired         = if ("bla_esbl_acquired" %in% names(.))           has_ctxm(.data[["bla_esbl_acquired"]])           else FALSE,
    hit_truncated_resistance_hits = if ("truncated_resistance_hits" %in% names(.))  has_ctxm(.data[["truncated_resistance_hits"]])  else FALSE,
    hit_spurious_resistance_hits  = if ("spurious_resistance_hits" %in% names(.))   has_ctxm(.data[["spurious_resistance_hits"]])   else FALSE,
    hit_rescue_esbl_list          = if ("rescue_esbl_list" %in% names(.))            has_ctxm(.data[["rescue_esbl_list"]])            else FALSE,
    hit_rescue_all_genes_genes    = if ("rescue_all_genes_genes" %in% names(.))      has_ctxm(.data[["rescue_all_genes_genes"]])      else FALSE,
    
    ctxm_truncated_tokens = ifelse(ctxm_trunc_tok_raw == "" & hit_truncated_resistance_hits, "CTX-M", ctxm_trunc_tok_raw),
    ctxm_spurious_tokens  = ifelse(ctxm_spur_tok_raw  == "" & hit_spurious_resistance_hits,  "CTX-M", ctxm_spur_tok_raw)
  )

esbl_hits_mat <- cbind(
  dat1$hit_bla_esbl_acquired,
  dat1$hit_truncated_resistance_hits,
  dat1$hit_spurious_resistance_hits,
  dat1$hit_rescue_esbl_list,
  dat1$hit_rescue_all_genes_genes
)
dat1$esbl_ctxm_any <- rowSums(as.matrix(esbl_hits_mat), na.rm = TRUE) > 0

# ---------------- SHV mutations (ESBL by mutation) + TEM status (from tem_is_esbl) ----------------
shv_mut_cols <- c("shv_mutations",
                  "rescue_klebsiella_pneumo_complex_amr_shv_mutations",
                  "shv_mutations.x","shv_mutations.y")

tem_mut_cols <- c("tem_mutations",
                  "rescue_klebsiella_pneumo_complex_amr_tem_mutations",
                  "tem_mutations.x","tem_mutations.y")

dat1 <- dat1 %>%
  mutate(
    shv_mut_text = collapse_mut_cols(., shv_mut_cols),
    tem_mut_text = collapse_mut_cols(., tem_mut_cols)
  )

# SHV ESBL by mutation: 238X OR 179X OR 169X OR 148X OR (240K AND 35Q)
shv_238X <- mut_any(dat1$shv_mut_text, 238)
shv_179X <- mut_any(dat1$shv_mut_text, 179)
shv_169X <- mut_any(dat1$shv_mut_text, 169)
shv_148X <- mut_any(dat1$shv_mut_text, 148)
shv_240K <- mut_to(dat1$shv_mut_text, 240, "K")
shv_35Q  <- mut_to(dat1$shv_mut_text,  35, "Q")

dat1$esbl_shv_mut_any <- (!is.na(dat1$shv_mut_text)) &
  (shv_238X | shv_179X | shv_169X | shv_148X | (shv_240K & shv_35Q))

# TEM ESBL status: use tem_is_esbl directly (no recomputation)
dat1$esbl_tem_any <- if ("tem_is_esbl" %in% names(dat1)) as_flag(dat1$tem_is_esbl) else FALSE

# ---------------- other rules (carbapenemase / porin / AmpC) ----------------
txt_r2 <- paste_cols_norm(dat1, r2)
txt_r3 <- paste_cols_norm(dat1, r3)
txt_r4 <- paste_cols_norm(dat1, r4)

# Carbapenemase: include OXA-204 (FIX)
pat_carb  <- regex("(NDM|KPC|VIM|IMP|OXA[-_ ]?(48|181|232|204|244))", ignore_case = TRUE)

pat_porin <- regex("OMPK?3[56]|OMPK?36GD|OMPK?35[-_ ]?\\d+%|OMPK?36[-_ ]?\\d+%", ignore_case = TRUE)

carb_any  <- !is.na(txt_r2) & str_detect(txt_r2, pat_carb)
porin_txt <- !is.na(txt_r3) & str_detect(txt_r3, pat_porin)

# ---------------- AmpC (ALIGNED TO YOUR WORKING AmpC CODE) ----------------
# choose basis
txt_ampc_basis <- if (use_bla_acquired_only_for_ampc && "bla_acquired" %in% names(dat1)) {
  message("AmpC detection basis: bla_acquired only")
  norm_ctx_text(dat1$bla_acquired)
} else {
  message("AmpC detection basis: ALL r4 sources (txt_r4)")
  txt_r4
}

ampc_dha <- hit(txt_ampc_basis, tok("DHA"))
ampc_cmy <- hit(txt_ampc_basis, tok("CMY"))
ampc_lap <- hit(txt_ampc_basis, tok("LAP"))

ampc_any <- ampc_dha | ampc_cmy | ampc_lap

# ---------------- porin “flag columns” if present (often 1/TRUE) ----------------
omp_flag_cols <- intersect(
  c("rescue_omp_k35_trunc","rescue_omp_k35_fs",
    "rescue_omp_k36_trunc","rescue_omp_k36_fs",
    "rescue_omp_k36_l3_g_dins"),
  names(dat1)
)

porin_flags <- if (length(omp_flag_cols)) {
  rowSums(as.data.frame(lapply(dat1[omp_flag_cols], function(z){
    vz <- tolower(as.character(z))
    (vz %in% c("1","true","t","yes","y")) | suppressWarnings(as.numeric(vz) == 1)
  })), na.rm = TRUE) > 0
} else FALSE

porin_any <- porin_txt | porin_flags


# ---------------- final mechanism categories (exact requested labels) ----------------
mech_levels <- c(
  "CTX-M-only",
  "CTX-M + SHV mutation ± TEM mutation",
  "porin mutation only",
  "AmpC beta-lactamase",
  "carbapenemase + porin mutation",
  "carbapenemase + porin mutation + AmpC",
  "carbapenemase only",
  "none"
)

dat_rules <- dat1 %>%
  mutate(
    # audit-friendly CTX-M tokens
    ctxm_from_trunc_or_spur = na_if(paste0(
      ifelse(ctxm_truncated_tokens != "", ctxm_truncated_tokens, ""),
      ifelse(ctxm_truncated_tokens != "" & ctxm_spurious_tokens != "", ";", ""),
      ifelse(ctxm_spurious_tokens  != "", ctxm_spurious_tokens, "")
    ), ""),
    
    carb_any  = carb_any,
    porin_any = porin_any,
    
    # AmpC flags (now aligned)
    ampc_any  = ampc_any,
    ampc_dha  = ampc_dha,
    ampc_cmy  = ampc_cmy,
    ampc_lap  = ampc_lap,
    
    # ESBL flags
    esbl_ctxm_any     = esbl_ctxm_any,
    esbl_shv_mut_any  = esbl_shv_mut_any,
    esbl_tem_any      = esbl_tem_any,
    
    # mutually exclusive final grouping (priority order matters)
    mechanism_group = case_when(
      carb_any & porin_any & ampc_any ~ "carbapenemase + porin mutation + AmpC",
      carb_any & porin_any            ~ "carbapenemase + porin mutation",
      carb_any                        ~ "carbapenemase only",
      
      esbl_ctxm_any & esbl_shv_mut_any ~ "CTX-M + SHV mutation ± TEM mutation",
      esbl_ctxm_any                    ~ "CTX-M-only",
      
      # porin-only must truly be porin without other classes
      porin_any & !(carb_any | ampc_any | esbl_ctxm_any | esbl_shv_mut_any | esbl_tem_any) ~ "porin mutation only",
      
      ampc_any ~ "AmpC beta-lactamase",
      
      TRUE ~ "none"
    ),
    mechanism_group = factor(mechanism_group, levels = mech_levels),
    
    # optional audit string (helps debugging)
    mech_rules = {
      labs <- cbind(esbl_ctxm_any, esbl_shv_mut_any, esbl_tem_any, carb_any, porin_any, ampc_any)
      apply(labs, 1, function(v){
        nm <- c("ESBL_CTXM","ESBL_SHVmut","ESBL_TEMflag","CARB","PORIN","AMPC")[which(v)]
        if (length(nm) == 0) "none" else paste(nm, collapse = " + ")
      })
    }
  )

# ---- ensure map_tarss is present (FOUND / NOT FOUND) if available ----
tarss_col <- grep("tarss", names(dat_rules), value = TRUE)[1]
dat_rules <- dat_rules %>%
  mutate(
    map_tarss = if (!is.na(tarss_col)) toupper(trimws(as.character(.data[[tarss_col]]))) else NA_character_
  )

# ---------------- export (Excel) ----------------
mechanism_flags <- dat_rules %>%
  mutate(id = dplyr::coalesce(id_export, id, genome_id)) %>%
  select(
    id_core_unified, id,
    map_tarss,
    site = any_of("site"),
    year = any_of("year"),
    cefotaxime_5,
    ceftazidime_10,
    pheno_3gc,
    mechanism_group,
    
    esbl_ctxm_any, esbl_shv_mut_any, esbl_tem_any,
    carb_any, porin_any,
    
    # AmpC audit columns (NEW)
    ampc_any, ampc_dha, ampc_cmy, ampc_lap,
    
    mech_rules,
    ctxm_from_trunc_or_spur,
    shv_mut_text, tem_mut_text,
    tem_is_esbl = any_of("tem_is_esbl"),
    starts_with("hit_")
  )

not_matching_full <- dat_rules %>% filter(mechanism_group == "none")

outfile <- "Mechanism_rules_CTXM_SHVmut_TEMflag_dedup.xlsx"
write_xlsx(
  list(
    Mechanism_flags       = mechanism_flags,
    Not_matching_FULLDATA = not_matching_full
  ),
  outfile
)

# ---------------- QC messages ----------------
message("Wrote: ", outfile)
message("N isolates (deduped) = ", nrow(dat_rules))

message("\nMechanism_group counts:")
print(table(dat_rules$mechanism_group, useNA = "ifany"))

message("\nKey flags:")
message("CTX-M any = ", sum(dat_rules$esbl_ctxm_any, na.rm = TRUE))
message("SHV ESBL-mutation any = ", sum(dat_rules$esbl_shv_mut_any, na.rm = TRUE))
message("TEM ESBL flag any = ", sum(dat_rules$esbl_tem_any, na.rm = TRUE))

message("\nAmpC QC (aligned):")
message("Any AmpC = ", sum(dat_rules$ampc_any, na.rm = TRUE))
message("DHA      = ", sum(dat_rules$ampc_dha, na.rm = TRUE))
message("CMY      = ", sum(dat_rules$ampc_cmy, na.rm = TRUE))
message("LAP      = ", sum(dat_rules$ampc_lap, na.rm = TRUE))


#---- plot
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(janitor)
  library(scales)
  library(readr)
})

# ----- input -----
mech_df <- readxl::read_excel(
  "/Users/danait/Mechanism_rules_CTXM_SHVmut_TEMflag_dedup.xlsx",
  sheet = "Mechanism_flags"
) %>%
  janitor::clean_names()

# ---- phenotype levels shown in the figure ----
lvl_pheno <- c("CAZ-only", "CTX+CAZ", "CTX-only")

# ---- mechanism levels (MATCH your new mechanism_group labels exactly) ----
lvl_mech <- c(
  "CTX-M-only",
  "CTX-M + SHV mutation ± TEM mutation",
  "porin mutation only",
  "AmpC beta-lactamase",
  "carbapenemase + porin mutation",
  "carbapenemase + porin mutation + AmpC",
  "carbapenemase only",
  "none"
)

# ---- palette (named to match lvl_mech) ----
# (Colour-blind friendly, stable; adjust if you have house style.)
pal_mech <- c(
  "CTX-M-only"                              = "#0072B2", # blue
  "CTX-M + SHV mutation ± TEM mutation"     = "#6A3D9A", # violet
  "porin mutation only"                     = "#CC79A7", # purple
  "AmpC beta-lactamase"                     = "#009E73", # green
  "carbapenemase + porin mutation"          = "#D55E00", # vermillion
  "carbapenemase + porin mutation + AmpC"   = "#E69F00", # orange
  "carbapenemase only"                      = "#56B4E9", # sky blue
  "none"                                    = "#999999"  # grey
)

# ---- basic column checks (helpful reviewer-proofing) ----
stopifnot("pheno_3gc" %in% names(mech_df))
stopifnot("mechanism_group"   %in% names(mech_df))

# Warn if there are unexpected mechanism labels (prevents silent NAs)
unknown_mech <- setdiff(unique(na.omit(as.character(mech_df$mechanism_group))), lvl_mech)
if (length(unknown_mech) > 0) {
  message("WARNING: Found unexpected mechanism_group labels:\n  - ",
          paste(unknown_mech, collapse = "\n  - "))
}

# ---- clean factors for plotting ----
mech_df_clean <- mech_df %>%
  mutate(
    pheno_3gc = factor(pheno_3gc, levels = c(lvl_pheno, "Neither")),
    
    # Ensure blanks/NA become "none" (optional but pragmatic)
    mechanism_group = ifelse(is.na(mechanism_group) | mechanism_group == "", "none", mechanism_group),
    mechanism_group = factor(mechanism_group, levels = lvl_mech)
  )

# ---- count first, then complete to include missing combos with n = 0 ----
mech_by_pheno <- mech_df_clean %>%
  filter(pheno_3gc %in% lvl_pheno) %>%  # drop "Neither"
  count(pheno_3gc, mechanism_group, name = "n") %>%
  complete(
    pheno_3gc = factor(lvl_pheno, levels = lvl_pheno),
    mechanism_group = factor(lvl_mech, levels = lvl_mech),
    fill = list(n = 0)
  ) %>%
  group_by(pheno_3gc) %>%
  mutate(
    denom    = sum(n),
    pct_frac = ifelse(denom > 0, n / denom, 0)
  ) %>%
  ungroup()

# Optional: export a supplement-ready table
readr::write_csv(mech_by_pheno, "mechanisms_by_phenotype.csv")

# ---- plot ----
p <- ggplot(mech_by_pheno,
            aes(x = pheno_3gc, y = pct_frac, fill = mechanism_group)) +
  geom_col(color = "white", linewidth = 0.25) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = pal_mech, drop = FALSE, name = "Mechanism") +
  labs(
    x = NULL,
    y = "Percent within phenotype",
    title = "Mechanism combinations by 3GCR phenotype"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold")
  )

# Add stacked counts (only for non-zero slices)
p1 <- p +
  geom_text(
    data = mech_by_pheno %>% filter(n > 0),
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3,
    colour = "white"
  )

print(p1)

ggsave("Fig_mechanism_by_phenotype.pdf", p1,
       width = 210, height = 90, units = "mm")


# ===========================
# Transmission analysis toolkit for K. pneumoniae
# ===========================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(igraph)
  library(ggraph)
})

# ---------------------------
# 0) User inputs (EDIT THESE)
# ---------------------------

# Path to SNP difference matrix (Excel; first column = isolate IDs; square matrix thereafter)
SNP_XLSX  <- "/Users/danait/Downloads/kpn_difference_matrix_matched_to_metadata.csv"

# Metadata source:
# - If you already have `dat` in your environment, set USE_DAT <- TRUE
# - Or point to a CSV/Excel with at least: isolate_id, hospital, hospital_ward, year/month/day (or collection_date), and ESBL flag.
meta_raw <- read_xlsx("/Users/danait/Dropbox/WGS_analysis/WGS_data_analysed/combined_metadata_qc.xlsx")  # e.g., "/path/to/meta.csv"

# Which columns in your metadata are the true IDs and fields
ID_COL_CANDIDATES <- c("id_export","id","genome_id","id_join","id")  # one of these must exist
HOSP_COL          <- "site...150"
WARD_COL          <- "hospital_ward"
ESBL_COL          <- "Any_ESBL"        # logical/0/1
DATE_COL          <- NULL              # if NULL, will build from year/month/day

YEAR_COL          <- "year"
MONTH_COL         <- "month"
DAY_COL           <- "day"             # if 'day' missing, will default to 15th of month

# Harmonize hospital names (optional)

# Analysis grid
SNP_GRID <- c(10, 12, 15, 18, 20, 25, 30, 35, 40)
WIN_GRID <- c(28, 56, 84)   # 4, 8, 12 weeks

# Primary cut-off for networks
PRIMARY_SNP <- 15
PRIMARY_WIN <- 56

Secondary_SNP <-20

# Output directory
OUTDIR <- "wgs_transmission_outputs"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR)

# ------------------------------------------------
# 1) Helpers
# ------------------------------------------------

pick_first_existing <- function(df, candidates) {
  ok <- candidates[candidates %in% names(df)]
  if (length(ok) == 0) stop("None of the ID candidate columns exist in metadata: ", paste(candidates, collapse=", "))
  ok[1]
}

to_logical01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  x <- as.character(x)
  case_when(
    x %in% c("1","TRUE","True","true","yes","Yes","Y") ~ 1L,
    x %in% c("0","FALSE","False","false","no","No","N") ~ 0L,
    TRUE ~ NA_integer_
  )
}

build_date <- function(df) {
  if (!is.null(DATE_COL) && DATE_COL %in% names(df)) {
    d <- suppressWarnings(as.Date(df[[DATE_COL]]))
  } else {
    y <- df[[YEAR_COL]]
    m <- if (MONTH_COL %in% names(df)) df[[MONTH_COL]] else 6
    dday <- if (DAY_COL %in% names(df)) df[[DAY_COL]] else 15
    d <- suppressWarnings(make_date(year = as.integer(y), month = as.integer(m), day = as.integer(dday)))
  }
  d
}

upper_triangle_pairs <- function(M) {
  M <- as.matrix(M)
  stopifnot(nrow(M) == ncol(M))
  ids <- rownames(M)
  if (is.null(ids) || anyDuplicated(ids)) stop("SNP matrix must have unique rownames as isolate IDs.")
  M[lower.tri(M, diag = TRUE)] <- NA
  coords <- which(!is.na(M), arr.ind = TRUE)
  tibble(
    from = ids[coords[, "row"]],
    to   = ids[coords[, "col"]],
    snp  = M[coords]
  )
}

components_table <- function(g) {
  mem <- components(g)$membership
  tibble(id = names(mem), cl = as.integer(mem)) %>% count(cl, name = "size")
}

score_threshold <- function(pair_df, meta, snp_thr, win_days) {
  e <- pair_df %>% 
    filter(snp <= snp_thr) %>%
    left_join(meta, by = c("from"="isolate_id")) %>% 
    rename(hosp1 = hospital, ward1 = ward, date1 = date) %>%
    left_join(meta, by = c("to"="isolate_id")) %>% 
    rename(hosp2 = hospital, ward2 = ward, date2 = date) %>%
    mutate(dt = abs(as.numeric(date1 - date2))) %>%
    filter(!is.na(dt), dt <= win_days)
  
  if (nrow(e) == 0) {
    return(tibble(
      snp = snp_thr, win = win_days, 
      n_pairs = 0, n_clusters = 0, n_iso_clustered = 0, med_size = NA_real_,
      pairs_same_ward = NA_real_, pairs_same_hosp = NA_real_,
      pct_clustered = 0, attributable = 0
    ))
  }
  
  n_pairs <- nrow(e)
  g <- graph_from_data_frame(e %>% select(from, to), directed = FALSE)
  cl_sizes <- components_table(g)
  n_iso <- sum(cl_sizes$size)
  med_size <- median(cl_sizes$size)
  
  # Concordance
  pairs_same_ward <- mean(e$hosp1 == e$hosp2 & e$ward1 == e$ward2, na.rm = TRUE)
  pairs_same_hosp <- mean(e$hosp1 == e$hosp2, na.rm = TRUE)
  
  # % clustered (of all isolates in metadata)
  pct_clustered <- n_iso / nrow(meta)
  
  # approximate "attributable to transmission": sum(size-1)/N
  attributable <- sum(pmax(cl_sizes$size - 1, 0)) / nrow(meta)
  
  tibble(
    snp = snp_thr, win = win_days, 
    n_pairs = n_pairs,
    n_clusters = nrow(cl_sizes),
    n_iso_clustered = n_iso,
    med_size = med_size,
    pairs_same_ward = pairs_same_ward,
    pairs_same_hosp = pairs_same_hosp,
    pct_clustered = pct_clustered,
    attributable = attributable
  )
}

# --- per-(snp,win) cluster stats incl. duration & mono-ward ---
stats_one <- function(pair_df, meta, snp_thr, win_days){
  e <- pair_df %>%
    filter(snp <= snp_thr) %>%
    left_join(meta, by = c("from"="isolate_id")) %>%
    rename(hosp1=hospital, ward1=ward, date1=date) %>%
    left_join(meta, by = c("to"  ="isolate_id")) %>%
    rename(hosp2=hospital, ward2=ward, date2=date) %>%
    mutate(dt = abs(as.numeric(date1 - date2))) %>%
    filter(!is.na(dt), dt <= win_days)
  
  if (nrow(e) == 0) {
    return(tibble(
      snp=snp_thr, win=win_days, n_pairs=0, n_clusters=0, n_iso_clustered=0,
      med_size=NA_real_, pct_clustered=0, attributable=0,
      med_span=NA_real_, iqr_span=NA_real_, mono_ward_pct=NA_real_
    ))
  }
  
  # graph → cluster membership
  g  <- graph_from_data_frame(e %>% select(from, to), directed=FALSE)
  mem <- components(g)$membership
  cl_sizes <- tibble(id = names(mem), cl = as.integer(mem)) %>% count(cl, name="size")
  
  n_iso <- sum(cl_sizes$size)
  med_size <- median(cl_sizes$size)
  pct_clustered <- n_iso / nrow(meta)
  attributable  <- sum(pmax(cl_sizes$size - 1, 0)) / nrow(meta)
  
  # per-cluster span & mono-ward
  per_cl <- tibble(id = names(mem), cl = as.integer(mem)) %>%
    left_join(meta %>% select(isolate_id, date, ward), by = c("id" = "isolate_id")) %>%
    group_by(cl) %>%
    summarise(span = as.numeric(max(date, na.rm=TRUE) - min(date, na.rm=TRUE)),
              n_wards = n_distinct(ward), .groups="drop")
  
  tibble(
    snp = snp_thr, win = win_days,
    n_pairs = nrow(e),
    n_clusters = nrow(cl_sizes),
    n_iso_clustered = n_iso,
    med_size = med_size,
    pct_clustered = pct_clustered,
    attributable = attributable,
    med_span = median(per_cl$span, na.rm=TRUE),
    iqr_span = IQR(per_cl$span, na.rm=TRUE),
    mono_ward_pct = mean(per_cl$n_wards == 1, na.rm=TRUE)
  )
}

permute_ward_time_once <- function(meta) {
  # shuffle ward within hospital × year × month strata (preserves margins)
  meta %>%
    mutate(yy = year(date), mm = month(date)) %>%
    group_by(hospital, yy, mm) %>%
    mutate(ward = sample(ward)) %>%
    ungroup() %>%
    select(-yy, -mm)
}

score_concordance_perm <- function(pair_df, meta, snp_thr, win_days, n_perm = 500, seed = 7) {
  set.seed(seed)
  # Observed
  obs <- score_threshold(pair_df, meta, snp_thr, win_days)
  obs_c <- obs$pairs_same_ward
  
  # Permutations
  perm_vals <- replicate(n_perm, {
    m2 <- permute_ward_time_once(meta)
    score_threshold(pair_df, m2, snp_thr, win_days)$pairs_same_ward
  })
  tibble(
    snp = snp_thr, win = win_days,
    obs_pairs_same_ward = obs_c,
    null_mean = mean(perm_vals, na.rm=TRUE),
    null_q95  = quantile(perm_vals, 0.95, na.rm=TRUE),
    null_q99  = quantile(perm_vals, 0.99, na.rm=TRUE)
  )
}

clusters_at <- function(pair_df, meta, snp_thr, win_days) {
  e <- pair_df %>% 
    filter(snp <= snp_thr) %>%
    left_join(meta, by = c("from"="isolate_id")) %>% 
    rename(date1 = date) %>%
    left_join(meta, by = c("to"="isolate_id")) %>% 
    rename(date2 = date) %>%
    mutate(dt = abs(as.numeric(date1 - date2))) %>%
    filter(!is.na(dt), dt <= win_days)
  if (nrow(e) == 0) return(tibble(id = character(), cl = integer()))
  g <- graph_from_data_frame(e %>% select(from, to), directed = FALSE)
  mem <- components(g)$membership
  tibble(id = names(mem), cl = as.integer(mem))
}

membership_agreement <- function(A, B) {
  # crude per-isolate agreement (same cluster label after aligning overlapping isolates)
  ids <- intersect(A$id, B$id)
  if (length(ids) == 0) return(NA_real_)
  A2 <- A %>% filter(id %in% ids) %>% arrange(id)
  B2 <- B %>% filter(id %in% ids) %>% arrange(id)
  mean(A2$cl == B2$cl)
}

# ------------------------------------------------
# 2) Load SNP matrix
# ------------------------------------------------

snp_xl <- read_csv(SNP_XLSX)
colnames(snp_xl)[1] <- "isolate_id"

# Keep unique rownames, coerce numeric
snp_mat <- snp_xl %>%
  filter(!duplicated(isolate_id)) %>%
  column_to_rownames("isolate_id") %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  as.matrix()

# Build pair list (upper triangle)
pair_df <- upper_triangle_pairs(snp_mat)

# ------------------------------------------------
# 3) Load / prepare metadata
# ------------------------------------------------

ID_COL <- pick_first_existing(meta_raw, ID_COL_CANDIDATES)

meta <- meta_raw %>%
  transmute(
    isolate_id = .data[[ID_COL]],
    hospital   = .data[[HOSP_COL]],
    ward       = .data[[WARD_COL]]
  ) %>%
  mutate(
    hospital = recode(hospital, !!!HOSP_RECODE)
  )

# Add date
meta$date <- build_date(meta_raw)
stopifnot(!is.null(meta$date))

# Clean and keep only isolates present in the SNP matrix
meta <- meta %>%
  filter(!is.na(isolate_id), isolate_id %in% rownames(snp_mat)) %>%
  mutate(
    isolate_id = as.character(isolate_id),
    ward = str_squish(as.character(ward)),
    hospital = str_squish(as.character(hospital))
  )

# ------------------------------------------------
# 4) Concordance grid evaluation
# ------------------------------------------------

grid <- expand_grid(snp = SNP_GRID, win = WIN_GRID)
scores <- pmap_dfr(list(grid$snp, grid$win), ~ score_threshold(pair_df, meta, ..1, ..2))

# Save
write_csv(scores, file.path(OUTDIR, "concordance_grid_scores.csv"))

# Quick “elbow” view (optional)
p_sameward <- ggplot(scores, aes(snp, pairs_same_ward, color = factor(win))) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = SNP_GRID) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(y = "% Same-ward", color = "Window (days)")

ggsave("same_ward_pairs.png")

# ---- grid of supervisor-requested metrics ----
grid <- expand_grid(snp = SNP_GRID, win = WIN_GRID)
grid_stats <- grid %>%
  pmap(~ stats_one(pair_df, meta, ..1, ..2)) %>%  # returns a list of tibbles
  list_rbind()                                    # purrr::list_rbind


write_csv(grid_stats, file.path(OUTDIR, "grid_cluster_stats.csv"))
library(scales)

# A1) % isolates in clusters
p_frac <- ggplot(grid_stats, aes(snp, pct_clustered, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "% isolates in clusters", x = "SNP threshold", color = "Window (days)")

p_frac

ggsave(file.path(OUTDIR, "plot_pct_clustered.png"), p_frac, width=7, height=5, dpi=300)

# A2) % attributable to transmission
p_attr <- ggplot(grid_stats, aes(snp, attributable, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "% attributable to transmission", x = "SNP threshold", color = "Window (days)")

p_attr

ggsave(file.path(OUTDIR, "plot_pct_attributable.png"), p_attr, width=7, height=5, dpi=300)

# B) Number of clusters & median cluster size
p_nclu <- ggplot(grid_stats, aes(snp, n_clusters, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "Number of clusters", x = "SNP threshold", color = "Window (days)")

p_nclu
ggsave(file.path(OUTDIR, "plot_n_clusters.png"), p_nclu, width=7, height=5, dpi=300)

p_medsize <- ggplot(grid_stats, aes(snp, med_size, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "Median cluster size", x = "SNP threshold", color = "Window (days)")

p_medsize

ggsave(file.path(OUTDIR, "plot_med_cluster_size.png"), p_medsize, width=7, height=5, dpi=300)

# C) Cluster duration (median ± IQR ribbon)
p_span <- ggplot(grid_stats, aes(snp, med_span, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  geom_ribbon(aes(ymin = med_span - iqr_span/2, ymax = med_span + iqr_span/2, fill = factor(win)),
              alpha = 0.15, color = NA) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "Median cluster span (days) ± IQR", x = "SNP threshold",
       color = "Window (days)", fill = "Window (days)")

p_span
ggsave(file.path(OUTDIR, "plot_cluster_span.png"), p_span, width=7, height=5, dpi=300)

# D) Mono-ward share (cluster-level specificity)
p_mono <- ggplot(grid_stats, aes(snp, mono_ward_pct, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "% clusters single-ward", x = "SNP threshold", color = "Window (days)")

p_mono

ggsave(file.path(OUTDIR, "plot_mono_ward_pct.png"), p_mono, width=7, height=5, dpi=300)


# Collect legends across all plots and set a single title/subtitle
# --- Compute the grid of stats (if you haven't already) ---
library(purrr); library(tidyr); library(scales); library(patchwork)

grid_stats <- expand_grid(snp = SNP_GRID, win = WIN_GRID) %>%
  mutate(dat = map2(snp, win, ~ stats_one(pair_df, meta, .x, .y))) %>%
  unnest(dat)

# --- 1) Cluster fraction (% isolates in clusters) ---
p_frac <- ggplot(grid_stats, aes(snp, pct_clustered, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "% isolates in clusters", x = "SNP threshold", color = "Window (days)")

# --- 2) Number of clusters ---
p_nclu <- ggplot(grid_stats, aes(snp, n_clusters, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "Number of clusters", x = "SNP threshold", color = "Window (days)")

# --- 3) Cluster duration (median span ± IQR) ---
p_span <- ggplot(grid_stats, aes(snp, med_span, color = factor(win), group = factor(win))) +
  geom_line(size=1) + geom_point(size=2) +
  geom_ribbon(aes(ymin = med_span - iqr_span/2, ymax = med_span + iqr_span/2, fill = factor(win)),
              alpha = 0.15, color = NA) +
  scale_x_continuous(breaks = SNP_GRID) +
  labs(y = "Median cluster span (days) ± IQR", x = "SNP threshold",
       color = "Window (days)", fill = "Window (days)")

# --- Put them together & save ---
combo <- (p_frac / p_nclu | p_span/ p_sameward) +
  plot_layout(guides = "collect", heights = c(1,1,1.2)) &
  theme(legend.position = "bottom")

combo
ggsave(file.path(OUTDIR, "threshold_window_summary.png"),
       combo, width = 10, height = 12, dpi = 300)

ggsave(file.path(OUTDIR, "threshold_window_summary.pdf"),
       combo, width = 10, height = 12, device = cairo_pdf)


# ------------------------------------------------
# 5) Permutation null at the primary threshold (ward/time concordance)
# ------------------------------------------------

perm_summary <- score_concordance_perm(pair_df, meta, PRIMARY_SNP, PRIMARY_WIN, n_perm = 500, seed = 7)
write_csv(perm_summary, file.path(OUTDIR, "perm_concordance_primary.csv"))

# ------------------------------------------------
# 6) Cluster stability between adjacent thresholds (at 8 weeks / 56 days)
# ------------------------------------------------

stab_df <- tibble()
thr_seq <- sort(SNP_GRID)
for (i in seq_len(length(thr_seq) - 1)) {
  A <- clusters_at(pair_df, meta, thr_seq[i], PRIMARY_WIN)
  B <- clusters_at(pair_df, meta, thr_seq[i+1], PRIMARY_WIN)
  stab_df <- bind_rows(stab_df, tibble(
    snp_a = thr_seq[i],
    snp_b = thr_seq[i+1],
    agreement = membership_agreement(A, B)
  ))
}
write_csv(stab_df, file.path(OUTDIR, "cluster_stability_pairs.csv"))

# ------------------------------------------------
# 7) Ward-level transmission networks per hospital (primary threshold)
# ------------------------------------------------

# Build edges satisfying PRIMARY_SNP/PRIMARY_WIN
edges_primary <- pair_df %>%
  filter(snp <= PRIMARY_SNP) %>%
  left_join(meta, by = c("from"="isolate_id")) %>% 
  rename(hosp1 = hospital, ward1 = ward, date1 = date, esbl1 = esbl) %>%
  left_join(meta, by = c("to"="isolate_id")) %>% 
  rename(hosp2 = hospital, ward2 = ward, date2 = date, esbl2 = esbl) %>%
  mutate(dt = abs(as.numeric(date1 - date2))) %>%
  filter(!is.na(dt), dt <= PRIMARY_WIN) %>%
  drop_na(hosp1, ward1, hosp2, ward2)

if (nrow(edges_primary) > 0) {
  # collapse to within-hospital, ward-to-ward edges (undirected)
  edges_by_hosp <- edges_primary %>%
    filter(hosp1 == hosp2) %>%
    transmute(hospital = hosp1,
              w1 = pmap_chr(list(ward1, ward2), ~ sort(c(..1, ..2))[1]),
              w2 = pmap_chr(list(ward1, ward2), ~ sort(c(..1, ..2))[2])) %>%
    group_by(hospital, w1, w2) %>%
    summarise(n_links = n(), .groups = "drop")
  
  # node attributes per ward
  clustered_ids <- unique(c(edges_primary$from, edges_primary$to))
  nodes_by_hosp <- meta %>%
    filter(isolate_id %in% clustered_ids) %>%
    group_by(hospital, ward) %>%
    summarise(
      n_clustered = n(),
      esbl_pct = mean(esbl, na.rm = TRUE),
      .groups = "drop"
    )
  
  out_net_dir <- file.path(OUTDIR, "networks_by_ward")
  if (!dir.exists(out_net_dir)) dir.create(out_net_dir)
  
  hosp_list <- sort(unique(edges_by_hosp$hospital))
  for (h in hosp_list) {
    e_h <- edges_by_hosp %>% filter(hospital == h) %>%
      mutate(w1 = str_squish(w1), w2 = str_squish(w2)) %>%
      filter(w1 != "", w2 != "") %>%
      distinct(w1, w2, .keep_all = TRUE)
    
    n_h <- nodes_by_hosp %>% filter(hospital == h) %>%
      transmute(ward = str_squish(ward), n_clustered, esbl_pct) %>%
      filter(ward != "") %>%
      distinct(ward, .keep_all = TRUE)
    
    wards_all <- sort(unique(c(e_h$w1, e_h$w2, n_h$ward)))
    n_h <- tibble(ward = wards_all) %>%
      left_join(n_h, by = "ward") %>%
      mutate(n_clustered = coalesce(n_clustered, 0L),
             esbl_pct = coalesce(esbl_pct, 0))
    
    stopifnot(!anyDuplicated(n_h$ward), !any(is.na(n_h$ward)))
    
    g <- graph_from_data_frame(
      d = e_h %>% select(from = w1, to = w2, n_links),
      vertices = n_h %>% rename(name = ward),
      directed = FALSE
    )
    
    p <- ggraph(g, layout = "fr") +
      geom_edge_link(aes(width = n_links), alpha = 0.35, color = "grey60") +
      geom_node_point(aes(size = n_clustered, fill = esbl_pct), shape = 21, color = "grey25", stroke = 0.4) +
      geom_node_text(aes(label = name), size = 3, repel = TRUE) +
      scale_size_continuous(name = "# clustered isolates", range = c(3, 14)) +
      scale_fill_gradientn(name = "ESBL+ proportion", colours = c("#F2F2F2", "#B22222"), limits = c(0,1), labels = scales::percent) +
      scale_edge_width_continuous(name = "# ward links", range = c(0.3, 2.5)) +
      labs(
        title = paste0("Ward-level transmission network: ", h, " (≤", PRIMARY_SNP, " SNP; ≤", PRIMARY_WIN, " days)"),
        subtitle = "Node size = # clustered isolates; fill = ESBL+ proportion; edge width = number of ward-to-ward links"
      ) +
      theme_void() +
      theme(legend.position = "right",
            plot.title = element_text(face="bold"))
    
    ggsave(file.path(out_net_dir, paste0("network_", h, "_SNP", PRIMARY_SNP, "_WIN", PRIMARY_WIN, ".png")),
           p, width = 10, height = 8, dpi = 300)
  }
}

# ---------------------------
# 8) Final console summaries
# ---------------------------

message("\nSaved:")
message(" - Concordance grid: ", file.path(OUTDIR, "concordance_grid_scores.csv"))
message(" - Permutation summary (primary): ", file.path(OUTDIR, "perm_concordance_primary.csv"))
message(" - Cluster stability pairs: ", file.path(OUTDIR, "cluster_stability_pairs.csv"))
message(" - Ward networks (PNGs): ", file.path(OUTDIR, "networks_by_ward"))

# ---------------------------
# 9) How to read outputs (quick tips)
# ---------------------------
# * 'concordance_grid_scores.csv' → pick the elbow: sharp rise 10→15, smaller gains after.
# * 'perm_concordance_primary.csv' → show observed ward-concordance >> null_q95 → your cut-off is not random proximity.
# * 'cluster_stability_pairs.csv' → high agreement 15↔20 means stable clusters around your chosen threshold.


########### The final code ################# this code doesnt include SHV or TEM 
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(readxl); library(janitor); library/writexl; library(stringi)
})

# ---------------- helpers ----------------
core_id <- function(x){
  s <- as.character(x)
  s <- gsub("^PID-\\d+-", "", s, ignore.case = TRUE)
  s <- gsub("_S\\d+.*$", "", s, ignore.case = TRUE)
  trimws(s)
}
coalesce_cols_chr <- function(df, candidates){
  cols <- intersect(candidates, names(df))
  if (!length(cols)) return(rep(NA_character_, nrow(df)))
  out <- as.character(df[[cols[1]]])
  if (length(cols) > 1) for (nm in cols[-1]) out <- dplyr::coalesce(out, as.character(df[[nm]]))
  out
}

# Unicode-robust normaliser for gene strings (keep only one definition)
norm_ctx_text <- function(v){
  x <- as.character(v); x[is.na(x)] <- ""
  x <- stringi::stri_trans_nfkc(x)                                # canonicalise
  x <- stringi::stri_replace_all_regex(x, "[\\p{Cf}\\p{Cc}]", "") # strip format/control (ZWSP, SHY…)
  x <- stringi::stri_replace_all_regex(x, "\\p{Z}+", " ")         # collapse spaces
  x <- stringi::stri_replace_all_regex(x, "[\\p{Pd}]", "-")       # any dash -> '-'
  x <- chartr("\u039C\u00B5\u0421\u0422\u0425\u041C", "MMCTXM", x) # Greek μ/M, Cyrillic C/T/X/M
  x <- stringi::stri_replace_all_regex(x, "[*^?]", "")            # strip ornaments
  toupper(x)
}

# Extract clean CTX-M tokens (e.g. "CTX-M-103;CTX-M-114"), empty "" if none
extract_ctxm_tokens <- function(v){
  x <- norm_ctx_text(v)
  hits <- stringi::stri_extract_all_regex(x, "(?iu)CTX[-_ ]?M[-_ ]?\\d+")
  vapply(hits, function(h){
    h <- unique(h[!is.na(h)])
    if (length(h)) paste(h, collapse = ";") else ""
  }, character(1))
}

# NA-safe, Unicode/dash-robust CTX-M presence detector
has_ctxm <- function(v){
  x <- as.character(v); x[is.na(x)] <- ""
  x <- stringi::stri_trans_nfkc(x)
  x <- stringi::stri_replace_all_regex(x, "[\\p{Cf}\\p{Cc}]", "")
  x <- stringi::stri_replace_all_regex(x, "\\p{Z}+", " ")
  x <- stringi::stri_replace_all_regex(x, "[\\p{Pd}]", "-")
  x <- chartr("\u039C\u00B5\u0421\u0422\u0425\u041C", "MMCTXM", x)
  x <- stringi::stri_replace_all_regex(x, "[*^?]", "")
  stringi::stri_detect_regex(x, "(?iu)CTX\\s*[-_]?\\s*M")
}

sir_norm <- function(x){
  y <- toupper(trimws(as.character(x)))
  y <- gsub("[\\*\\^\\?\\+]+$", "", y)
  y[!y %in% c("S","I","R")] <- NA_character_
  y
}

# Paste & normalise multiple columns
paste_cols_norm <- function(d, cols){
  cols <- intersect(cols, names(d))
  if (!length(cols)) return(rep(NA_character_, nrow(d)))
  m <- as.data.frame(lapply(d[cols], norm_ctx_text))
  txt <- do.call(paste, c(m, list(sep=";")))
  txt[txt %in% c("", "NA")] <- NA
  txt
}

# ---------------- unified ID + dedup BEFORE rules ----------------
dat0 <- dat_withrescued %>%
  mutate(
    id_unified      = coalesce_cols_chr(., c("id_export","id","genome_id","sample","strain","isolate","isolate_id")),
    id_core_unified = core_id(id_unified)
  )

has_rescue_col <- "rescue_all_genes_genes" %in% names(dat0)
dat1 <- dat0 %>%
  arrange(desc(if (has_rescue_col) !is.na(.data[["rescue_all_genes_genes"]]) else FALSE),
          id_core_unified, id_unified) %>%
  distinct(id_core_unified, .keep_all = TRUE)

# ---------------- rule sources ----------------
r1 <- c("bla_esbl_acquired","truncated_resistance_hits","spurious_resistance_hits",
        "rescue_esbl_list","rescue_all_genes_genes")   # ESBL = CTX-M only
r2 <- c("bla_carb_acquired","rescue_carb_list","rescue_carb_family",
        "rescue_klebsiella_pneumo_complex_amr_bla_carb_acquired")
r3 <- c("omp_mutations","rescue_klebsiella_pneumo_complex_amr_omp_mutations")
r4 <- c("bla_acquired","klebsiella_pneumo_complex__amr__bla_acquired","beta_lactam_genes_final",
        "rescue_beta_lactam_genes","rescue_all_genes_genes",
        "rescue_klebsiella_pneumo_complex_amr_bla_acquired")

# ---------------- CTX-M detection (NA-safe & unified) ----------------
dat1 <- dat1 %>%
  mutate(
    # tokens (for reporting)
    ctxm_trunc_tok_raw = extract_ctxm_tokens(.data[["truncated_resistance_hits"]]),
    ctxm_spur_tok_raw  = extract_ctxm_tokens(.data[["spurious_resistance_hits"]]),
    
    # robust boolean hits per column
    hit_bla_esbl_acquired      = if ("bla_esbl_acquired" %in% names(.))      has_ctxm(.data[["bla_esbl_acquired"]])           else FALSE,
    hit_truncated_resistance_hits = if ("truncated_resistance_hits" %in% names(.)) has_ctxm(.data[["truncated_resistance_hits"]]) else FALSE,
    hit_spurious_resistance_hits  = if ("spurious_resistance_hits" %in% names(.))  has_ctxm(.data[["spurious_resistance_hits"]])  else FALSE,
    hit_rescue_esbl_list       = if ("rescue_esbl_list" %in% names(.))       has_ctxm(.data[["rescue_esbl_list"]])            else FALSE,
    hit_rescue_all_genes_genes = if ("rescue_all_genes_genes" %in% names(.)) has_ctxm(.data[["rescue_all_genes_genes"]])      else FALSE,
    
    # final tokens (use generic 'CTX-M' if digits absent but presence detected)
    ctxm_truncated_tokens = ifelse(ctxm_trunc_tok_raw == "" & hit_truncated_resistance_hits, "CTX-M", ctxm_trunc_tok_raw),
    ctxm_spurious_tokens  = ifelse(ctxm_spur_tok_raw  == "" & hit_spurious_resistance_hits,  "CTX-M", ctxm_spur_tok_raw)
  )

# NA-safe rowwise OR across all CTX-M sources
esbl_hits_mat <- cbind(
  dat1$hit_bla_esbl_acquired,
  dat1$hit_truncated_resistance_hits,
  dat1$hit_spurious_resistance_hits,
  dat1$hit_rescue_esbl_list,
  dat1$hit_rescue_all_genes_genes
)
dat1$esbl_ctxm_any <- rowSums(as.matrix(esbl_hits_mat), na.rm = TRUE) > 0

# ---------------- other rules ----------------
txt_r2 <- paste_cols_norm(dat1, r2)
txt_r3 <- paste_cols_norm(dat1, r3)
txt_r4 <- paste_cols_norm(dat1, r4)

pat_carb  <- regex("(NDM|KPC|VIM|IMP|OXA[-_ ]?(48|181|232|244))", ignore_case = TRUE)
pat_porin <- regex("OMPK?3[56]|OMPK?36GD|OMPK?35[-_ ]?\\d+%|OMPK?36[-_ ]?\\d+%", ignore_case = TRUE)
pat_ampc  <- regex("\\b(DHA|CMY|LAP)\\b", ignore_case = TRUE)

carb_any  <- !is.na(txt_r2) & str_detect(txt_r2, pat_carb)
porin_txt <- !is.na(txt_r3) & str_detect(txt_r3, pat_porin)

omp_flag_cols <- intersect(c("rescue_omp_k35_trunc","rescue_omp_k35_fs",
                             "rescue_omp_k36_trunc","rescue_omp_k36_fs",
                             "rescue_omp_k36_l3_g_dins"), names(dat1))

# robust coercion for porin flag columns that may be 1/TRUE/"TRUE"/"1"/yes
porin_flags <- if (length(omp_flag_cols)) {
  rowSums(as.data.frame(lapply(dat1[omp_flag_cols], function(z){
    vz <- tolower(as.character(z))
    (vz %in% c("1","true","t","yes","y")) | suppressWarnings(as.numeric(vz) == 1)
  })), na.rm = TRUE) > 0
} else FALSE

porin_any <- porin_txt | porin_flags
ampc_any  <- !is.na(txt_r4) & str_detect(txt_r4, pat_ampc)

# ---------------- assemble ----------------
dat_rules <- dat1 %>%
  mutate(
    # report normalised CTX-M tokens from TRUNC/SPUR
    ctxm_from_trunc_or_spur = na_if(paste0(
      ifelse(ctxm_truncated_tokens != "", ctxm_truncated_tokens, ""),
      ifelse(ctxm_truncated_tokens != "" & ctxm_spurious_tokens != "", ";", ""),
      ifelse(ctxm_spurious_tokens  != "", ctxm_spurious_tokens, "")
    ), ""),
    
    carb_any      = carb_any,
    porin_any     = porin_any,
    ampc_any      = ampc_any,
    
    mech_rules = apply(cbind(esbl_ctxm_any, carb_any, porin_any, ampc_any), 1, function(v){
      labs <- c("ESBL_CTXM","CARB","PORIN","AMPC")[which(v)]
      if (length(labs) == 0) "none" else paste(labs, collapse = " + ")
    })
  )


# ---------------- export ----------------
mechanism_flags <- dat_rules %>%
  mutate(id = dplyr::coalesce(id_export, id, genome_id)) %>%
  select(id_core_unified, id, site, year = dplyr::any_of("year"),
         pheno_3gc,
         esbl_ctxm_any, carb_any, porin_any, ampc_any, mech_rules,
         starts_with("hit_"), ctxm_from_trunc_or_spur)

not_matching_full <- dat_rules %>% filter(mech_rules == "none")

writexl::write_xlsx(
  list(
    Mechanism_flags       = mechanism_flags,
    Not_matching_FULLDATA = not_matching_full
  ),
  "Mechanism_rules_CTXM_NO_THRESHOLD_dedup.xlsx"
)

message("Mechanism CTX-M total = ", sum(dat_rules$esbl_ctxm_any, na.rm=TRUE), " of ", nrow(dat_rules))
message("CTX-M hits — bla_esbl=", sum(dat_rules$hit_bla_esbl_acquired,      na.rm=TRUE),
        ", truncated=",            sum(dat_rules$hit_truncated_resistance_hits, na.rm=TRUE),
        ", spurious=",             sum(dat_rules$hit_spurious_resistance_hits,  na.rm=TRUE),
        ", rescue_esbl=",          sum(dat_rules$hit_rescue_esbl_list,          na.rm=TRUE),
        ", rescue_all=",           sum(dat_rules$hit_rescue_all_genes_genes,    na.rm=TRUE))




