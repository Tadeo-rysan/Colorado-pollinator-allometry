# ============================================================================
# COLORADO 2025: KOMPLETNI SPOJENY PIPELINE ANALYZ H1 & H2 (BASE R + lme4)
# ============================================================================
# Tento jeden skript propojuje vsechny dilci casti do jednoho uceleneho behu:
#   CAST 1: Nacteni dat H1/H2 a prirazeni celedi a barvy (Colorado_celedi1.xlsx)
#   CAST 2: Globalni modely H1 (alometrie b vs. izometrie b = 1)
#   CAST 3: Samostatne modely pro modelove druhy (n_obs >= 12, n_pos >= 7, visits >= 15)
#   CAST 4: Modely po celedich (n_species >= 3, n_obs >= 30, + Apiaceae) + krivky y(x) a p = y/x
#   CAST 5: Hypoteza H2 - Prekryv spekter opylovacu (sim_order, sim_fine)
#   CAST 6: Hypoteza H2 - Vliv stejne vs. jine celedi (interakce log_other * prop_same_fam)
#   CAST 7: Hypoteza H2 - Vliv stejne vs. jine barvy kvetu (color human)
#   CAST 8: Hypoteza H2 - Kombinace 4 skupin (stejna/jina celed x stejna/jina barva)
#
# Vse psano v cistem Base R + readxl + lme4 (bez tidyverse / glmmTMB).
# ============================================================================

options(stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------
# 0. NASTAVENI CEST, SLOZEK A BALICKU
# ----------------------------------------------------------------------------
BASE <- "C:/Users/Tadea/Documents/colorado"
if (!dir.exists(BASE)) {
  message("Vyber libovolny soubor ze slozky projektu Colorado:")
  BASE <- dirname(file.choose())
}
OUT <- file.path(BASE, "VYSLEDKY_LMER_BASE_R")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("readxl", "lme4")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(readxl)
library(lme4)

# Pomocne funkce pro testovani sklonu a export koeficientu
test_slope <- function(model, term, b0 = 1) {
  cf <- summary(model)$coefficients
  est <- cf[term, "Estimate"]
  se  <- cf[term, "Std. Error"]
  t_val <- (est - b0) / se
  df_approx <- nrow(model@frame) - nrow(cf)
  p_val <- 2 * pt(abs(t_val), df = df_approx, lower.tail = FALSE)
  data.frame(
    est = round(est, 4),
    SE  = round(se, 4),
    t   = round(t_val, 3),
    df  = df_approx,
    p   = signif(p_val, 4)
  )
}

extract_coefs <- function(model, model_name) {
  cf <- as.data.frame(summary(model)$coefficients)
  cf$term <- rownames(cf)
  cf$model <- model_name
  df_approx <- nrow(model@frame) - nrow(cf)
  cf$p_value <- signif(2 * pt(abs(cf[["t value"]]), df = df_approx, lower.tail = FALSE), 4)
  rownames(cf) <- NULL
  cf[, c("model", "term", "Estimate", "Std. Error", "t value", "p_value")]
}

find_file <- function(filename, prompt_label) {
  p1 <- file.path(OUT, filename)
  if (file.exists(p1)) return(p1)
  p2 <- file.path(BASE, filename)
  if (file.exists(p2)) return(p2)
  message("Soubor '", filename, "' nenalezen automaticky. ", prompt_label)
  file.choose()
}

# ============================================================================
# CAST 1: NACTENI DAT A VYPOCET OKOLNICH PROMENNYCH (CELED x BARVA)
# ============================================================================
cat("\n[1/8] Nacitam data H1, H2 a tabulku celedi + barev...\n")

h1_path <- find_file("01_H1_data_pro_lmer.csv", "Vyber 01_H1_data_pro_lmer.csv:")
H1_dat <- read.csv(h1_path)
H1_dat$Plant  <- trimws(as.character(H1_dat$Plant))
H1_dat$Site   <- factor(H1_dat$Site)
H1_dat$week_f <- factor(H1_dat$week)

if (!("visits_per_run" %in% names(H1_dat))) {
  H1_dat$visits_per_run <- H1_dat$visits / H1_dat$n_runs
}
H1_dat$log_rate  <- log(H1_dat$visits_per_run + 1)
H1_dat$log_own   <- log(H1_dat$flowers)
H1_dat$log_other <- log1p(H1_dat$other_flowers)

celedi_path <- find_file("Colorado_celedi1.xlsx", "Vyber Colorado_celedi1.xlsx:")
celedi <- as.data.frame(read_excel(celedi_path, sheet = "Plant_Family"))
names(celedi) <- trimws(names(celedi))

col_idx <- grep("^color", names(celedi), ignore.case = TRUE)
if (length(col_idx) == 0) {
  stop("V tabulce celedi nebyl nalezen sloupec s barvou ('color human').")
}
celedi$Plant       <- trimws(as.character(celedi$Plant))
celedi$Family      <- trimws(as.character(celedi$Family))
celedi$color_human <- trimws(as.character(celedi[[col_idx[1]]]))

H1_dat <- merge(H1_dat, celedi[, c("Plant", "Family", "color_human")], by = "Plant", all.x = TRUE)
if (any(is.na(H1_dat$Family)) || any(is.na(H1_dat$color_human))) {
  warning("Nektere druhy nemaji prirazenou celed nebo barvu!")
}

# Vypocet okolnich kvetu podle celedi, barvy a jejich 4 vzajemne se vylucujicich kombinaci
tot_fam     <- ave(H1_dat$flowers, paste(H1_dat$block, H1_dat$Family, sep = "__"), FUN = sum)
tot_col     <- ave(H1_dat$flowers, paste(H1_dat$block, H1_dat$color_human, sep = "__"), FUN = sum)
tot_fam_col <- ave(H1_dat$flowers, paste(H1_dat$block, H1_dat$Family, H1_dat$color_human, sep = "__"), FUN = sum)

# 1) Podil stejne vs. jine celedi
H1_dat$same_fam_other_flowers <- pmax(0, tot_fam - H1_dat$flowers)
H1_dat$diff_fam_flowers       <- pmax(0, H1_dat$other_flowers - H1_dat$same_fam_other_flowers)
H1_dat$prop_same_fam          <- ifelse(H1_dat$other_flowers > 0,
                                        H1_dat$same_fam_other_flowers / H1_dat$other_flowers, 0)
H1_dat$log_same_fam           <- log1p(H1_dat$same_fam_other_flowers)
H1_dat$log_diff_fam           <- log1p(H1_dat$diff_fam_flowers)

# 2) Podil stejne vs. jine barvy
H1_dat$same_col_other_flowers <- pmax(0, tot_col - H1_dat$flowers)
H1_dat$diff_col_flowers       <- pmax(0, H1_dat$other_flowers - H1_dat$same_col_other_flowers)
H1_dat$prop_same_col          <- ifelse(H1_dat$other_flowers > 0,
                                        H1_dat$same_col_other_flowers / H1_dat$other_flowers, 0)
H1_dat$log_same_col           <- log1p(H1_dat$same_col_other_flowers)
H1_dat$log_diff_col           <- log1p(H1_dat$diff_col_flowers)

# 3) 4 kombinace (stejna/jina celed x stejna/jina barva)
H1_dat$same_fam_same_col <- pmax(0, tot_fam_col - H1_dat$flowers)
H1_dat$same_fam_diff_col <- pmax(0, tot_fam - tot_fam_col)
H1_dat$diff_fam_same_col <- pmax(0, tot_col - tot_fam_col)
H1_dat$diff_fam_diff_col <- pmax(0, H1_dat$other_flowers - H1_dat$same_fam_same_col -
                                    H1_dat$same_fam_diff_col - H1_dat$diff_fam_same_col)

for (v in c("same_fam_same_col", "same_fam_diff_col", "diff_fam_same_col", "diff_fam_diff_col")) {
  H1_dat[[paste0("prop_", v)]] <- ifelse(H1_dat$other_flowers > 0, H1_dat[[v]] / H1_dat$other_flowers, 0)
  H1_dat[[paste0("log_", v)]]  <- log1p(H1_dat[[v]])
}

H1_dat$Plant <- factor(H1_dat$Plant)
write.csv(H1_dat, file.path(OUT, "01b_kompletni_matice_H1_celed_barva.csv"), row.names = FALSE)

# ============================================================================
# CAST 2: GLOBALNI MODELY H1 (VSECHNY DRUHY DOHROMADY)
# ============================================================================
cat("[2/8] Pocitam globalni modely H1...\n")

m1_own <- lmer(log_rate ~ log_own + (1 | week_f) + (1 | Site) + (1 | Plant),
               data = H1_dat, REML = TRUE)
m2_own_other <- lmer(log_rate ~ log_own + log_other + (1 | week_f) + (1 | Site) + (1 | Plant),
                     data = H1_dat, REML = TRUE)

t_g1 <- test_slope(m1_own,       "log_own",   b0 = 1)
t_g2 <- test_slope(m2_own_other, "log_own",   b0 = 1)
t_go <- test_slope(m2_own_other, "log_other", b0 = 0)

tab_H1_global <- rbind(
  cbind(model = "m1_own",       term = "log_own",   H0 = 1, t_g1),
  cbind(model = "m2_own_other", term = "log_own",   H0 = 1, t_g2),
  cbind(model = "m2_own_other", term = "log_other", H0 = 0, t_go)
)
write.csv(tab_H1_global, file.path(OUT, "02_H1_testy_alometrie_b.csv"), row.names = FALSE)

# ============================================================================
# CAST 3: SAMOSTATNE MODELY PRO MODELOVE DRUHY
# ============================================================================
cat("[3/8] Pocitam modely pro jednotlive modelove druhy...\n")

sp_tab <- aggregate(visits ~ Plant + Family + color_human, data = H1_dat, FUN = function(v) c(
  n_obs  = length(v),
  n_pos  = sum(v > 0),
  visits = sum(v)
))
sp_summary <- data.frame(
  Plant       = as.character(sp_tab$Plant),
  Family      = sp_tab$Family,
  color_human = sp_tab$color_human,
  n_obs       = sp_tab$visits[, "n_obs"],
  n_pos       = sp_tab$visits[, "n_pos"],
  visits      = sp_tab$visits[, "visits"]
)
sp_summary <- sp_summary[order(-sp_summary$n_pos, -sp_summary$n_obs), ]
write.csv(sp_summary, file.path(OUT, "06a_prehled_druhu_pocty_dat.csv"), row.names = FALSE)

vybrane_druhy <- subset(sp_summary, n_obs >= 12 & n_pos >= 7 & visits >= 15)$Plant

sp_models_list <- list()
for (sp in vybrane_druhy) {
  sub <- subset(H1_dat, as.character(Plant) == sp)
  re_parts <- c()
  if (length(unique(sub$week_f)) > 1) re_parts <- c(re_parts, "(1 | week_f)")
  if (length(unique(sub$Site))   > 1) re_parts <- c(re_parts, "(1 | Site)")
  re_str <- paste(re_parts, collapse = " + ")

  m1 <- lmer(as.formula(paste("log_rate ~ log_own +", re_str)), data = sub, REML = TRUE)
  m2 <- lmer(as.formula(paste("log_rate ~ log_own + log_other +", re_str)), data = sub, REML = TRUE)

  t1_own <- test_slope(m1, "log_own",   b0 = 1)
  t2_own <- test_slope(m2, "log_own",   b0 = 1)
  t2_oth <- test_slope(m2, "log_other", b0 = 0)

  sp_models_list[[sp]] <- data.frame(
    Plant       = sp,
    Family      = sub$Family[1],
    color_human = sub$color_human[1],
    n_obs       = nrow(sub),
    n_pos       = sum(sub$visits > 0),
    visits      = sum(sub$visits),
    m1_int      = round(fixef(m1)["(Intercept)"], 4),
    m1_b_own    = t1_own$est,
    m1_SE       = t1_own$SE,
    m1_p_b1     = t1_own$p,
    m2_b_own    = t2_own$est,
    m2_SE_own   = t2_own$SE,
    m2_p_b1     = t2_own$p,
    m2_b_oth    = t2_oth$est,
    m2_SE_oth   = t2_oth$SE,
    m2_p_oth0   = t2_oth$p
  )
}
sp_results <- do.call(rbind, sp_models_list)
write.csv(sp_results, file.path(OUT, "06_modely_jednotlive_druhy.csv"), row.names = FALSE)

# ============================================================================
# CAST 4: MODELY PO CELEDICH + PDF KRIVEK (y(x) A p = y/x)
# ============================================================================
cat("[4/8] Pocitam modely po celedich a vykresluji alometricke krivky...\n")

fam_tab <- aggregate(visits ~ Family, data = H1_dat, FUN = function(v) c(
  n_obs  = length(v),
  n_pos  = sum(v > 0),
  visits = sum(v)
))
fam_sp <- aggregate(as.character(Plant) ~ Family, data = H1_dat, FUN = function(p) length(unique(p)))
names(fam_sp)[2] <- "n_species"

fam_summary <- merge(fam_sp, data.frame(
  Family = fam_tab$Family,
  n_obs  = fam_tab$visits[, "n_obs"],
  n_pos  = fam_tab$visits[, "n_pos"],
  visits = fam_tab$visits[, "visits"]
), by = "Family")
fam_summary <- fam_summary[order(-fam_summary$n_obs), ]
write.csv(fam_summary, file.path(OUT, "07a_prehled_celedi_pocty_dat.csv"), row.names = FALSE)

vybrane_celedi <- subset(fam_summary, (n_species >= 3 & n_obs >= 30 & n_pos >= 10) | Family == "Apiaceae")$Family

fam_models_list <- list()
for (fam in vybrane_celedi) {
  sub <- subset(H1_dat, Family == fam)
  re_parts <- c()
  if (length(unique(sub$week_f)) > 1) re_parts <- c(re_parts, "(1 | week_f)")
  if (length(unique(sub$Site))   > 1) re_parts <- c(re_parts, "(1 | Site)")
  if (length(unique(sub$Plant))  > 1) re_parts <- c(re_parts, "(1 | Plant)")
  re_str <- paste(re_parts, collapse = " + ")

  m1 <- lmer(as.formula(paste("log_rate ~ log_own +", re_str)), data = sub, REML = TRUE)
  m2 <- lmer(as.formula(paste("log_rate ~ log_own + log_other +", re_str)), data = sub, REML = TRUE)

  t1_own <- test_slope(m1, "log_own",   b0 = 1)
  t2_own <- test_slope(m2, "log_own",   b0 = 1)
  t2_oth <- test_slope(m2, "log_other", b0 = 0)

  fam_models_list[[fam]] <- data.frame(
    Family    = fam,
    n_sp      = length(unique(sub$Plant)),
    n_obs     = nrow(sub),
    n_pos     = sum(sub$visits > 0),
    visits    = sum(sub$visits),
    m1_int    = round(fixef(m1)["(Intercept)"], 4),
    m1_b_own  = t1_own$est,
    m1_SE     = t1_own$SE,
    m1_p_b1   = t1_own$p,
    m2_b_own  = t2_own$est,
    m2_SE_own = t2_own$SE,
    m2_p_b1   = t2_own$p,
    m2_b_oth  = t2_oth$est,
    m2_SE_oth = t2_oth$SE,
    m2_p_oth0 = t2_oth$p
  )
}
fam_results <- do.call(rbind, fam_models_list)
write.csv(fam_results, file.path(OUT, "07_modely_celedi.csv"), row.names = FALSE)

# Vykresleni PDF s krivkami y(x) a p = y/x podle skript-JH.R
pdf(file.path(OUT, "08_krivky_druhy_a_celedi.pdf"), width = 15, height = 7)
par(mfrow = c(2, 5), mar = c(4, 4, 3.2, 1))
for (i in seq_len(min(10, nrow(sp_results)))) {
  row <- sp_results[i, ]
  sub <- subset(H1_dat, as.character(Plant) == row$Plant)
  x <- seq(1, max(20, max(sub$flowers)), length.out = 200)
  y <- pmax(0, (exp(row$m1_int) * ((x + 1)^row$m1_b_own)) - 1)
  plot(sub$flowers, sub$visits_per_run, pch = 16, col = rgb(0.2, 0.4, 0.7, 0.45),
       xlab = "pocet kvetu (x)", ylab = "navstevy / pruchod (y)",
       main = paste0(row$Plant, "\n(", row$Family, ", b = ", row$m1_b_own, ")"))
  lines(x, y, col = "firebrick", lwd = 2)
}

par(mfrow = c(2, 5), mar = c(4, 4, 3.2, 1))
for (i in seq_len(min(10, nrow(sp_results)))) {
  row <- sp_results[i, ]
  sub <- subset(H1_dat, as.character(Plant) == row$Plant)
  x <- seq(1, max(20, max(sub$flowers)), length.out = 200)
  y <- pmax(0, (exp(row$m1_int) * ((x + 1)^row$m1_b_own)) - 1)
  p <- y / x
  plot(x, p, type = "l", col = "darkgreen", lwd = 2,
       xlab = "pocet kvetu (x)", ylab = "navstevy na 1 kvet (y/x)",
       main = paste0(row$Plant, "\n(b = ", row$m1_b_own, ", p_b1 = ", row$m1_p_b1, ")"))
}

par(mfrow = c(2, 4), mar = c(4, 4, 3.2, 1))
for (i in seq_len(nrow(fam_results))) {
  row <- fam_results[i, ]
  sub <- subset(H1_dat, Family == row$Family)
  x <- seq(1, max(20, max(sub$flowers)), length.out = 200)
  y <- pmax(0, (exp(row$m1_int) * ((x + 1)^row$m1_b_own)) - 1)
  plot(sub$flowers, sub$visits_per_run, pch = 16, col = rgb(0.2, 0.5, 0.3, 0.4),
       xlab = "pocet kvetu (x)", ylab = "navstevy / pruchod (y)",
       main = paste0(row$Family, " (", row$n_sp, " dr.)\n(b = ", row$m1_b_own, ")"))
  lines(x, y, col = "firebrick", lwd = 2)
}
dev.off()

# ============================================================================
# CAST 5: HYPOTEZA H2 - PREKRYV SPEKTER OPYLOVACU (pokud existuje 03_H2_data)
# ============================================================================
h2_candidate1 <- file.path(OUT, "03_H2_data_pro_lmer.csv")
h2_candidate2 <- file.path(BASE, "03_H2_data_pro_lmer.csv")
h2_path <- if (file.exists(h2_candidate1)) h2_candidate1 else if (file.exists(h2_candidate2)) h2_candidate2 else ""

if (nzchar(h2_path)) {
  cat("[5/8] Pocitam modely H2 s prekryvem spekter opylovacu (sim_order, sim_fine)...\n")
  H2_dat <- read.csv(h2_path)
  H2_dat$Plant  <- factor(trimws(as.character(H2_dat$Plant)))
  H2_dat$Site   <- factor(H2_dat$Site)
  H2_dat$week_f <- factor(H2_dat$week)

  m_H2_order <- lmer(log_rate ~ log_own + log_other + sim_order + cov_tax +
                                (1 | week_f) + (1 | Site) + (1 | Plant),
                     data = H2_dat, REML = TRUE)
  m_H2_fine  <- lmer(log_rate ~ log_own + log_other + sim_fine + cov_tax +
                                (1 | week_f) + (1 | Site) + (1 | Plant),
                     data = H2_dat, REML = TRUE)

  tab_H2_poll <- rbind(
    extract_coefs(m_H2_order, "H2_prekryv_order"),
    extract_coefs(m_H2_fine,  "H2_prekryv_fine")
  )
  write.csv(tab_H2_poll, file.path(OUT, "04_H2_koeficienty_lmer.csv"), row.names = FALSE)
} else {
  cat("[5/8] Soubor 03_H2_data_pro_lmer.csv nenalezen, preskakuji cast 5.\n")
}

# ============================================================================
# CAST 6: HYPOTEZA H2 - PODIL STEJNE vs. JINE CELEDI + GRAF
# ============================================================================
cat("[6/8] Pocitam modely H2 pro podil stejne vs. jine celedi...\n")

m_fam_add   <- lmer(log_rate ~ log_own + log_other + prop_same_fam +
                               (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)
m_fam_int   <- lmer(log_rate ~ log_own + log_other * prop_same_fam +
                               (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)
m_fam_split <- lmer(log_rate ~ log_own + log_same_fam + log_diff_fam +
                               (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)

tab_fam_H2 <- rbind(
  extract_coefs(m_fam_add,   "M1_aditivni_prop_same_fam"),
  extract_coefs(m_fam_int,   "M2_interakce_log_other_x_prop_same_fam"),
  extract_coefs(m_fam_split, "M3_rozdeleni_same_vs_diff_fam")
)
write.csv(tab_fam_H2, file.path(OUT, "09_H2_podil_stejne_vs_jine_celedi.csv"), row.names = FALSE)

pdf(file.path(OUT, "10_graf_H2_podil_stejne_celedi.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
cf2 <- fixef(m_fam_int)
x_oth <- seq(min(H1_dat$log_other), max(H1_dat$log_other), length.out = 100)
m_own <- mean(H1_dat$log_own)
cols  <- c("steelblue", "forestgreen", "darkorange", "firebrick")
props <- c(0, 0.25, 0.50, 0.75)

plot(NULL, xlim = range(x_oth), ylim = c(0, 1.0),
     xlab = "log(pocet kvetu ostatnich druhu + 1)",
     ylab = "Predikovana log-navstevnost",
     main = "A) Interakce abundance a podilu stejne celedi")
for (k in seq_along(props)) {
  p_val <- props[k]
  y_hat <- cf2["(Intercept)"] + cf2["log_own"] * m_own + cf2["log_other"] * x_oth +
           cf2["prop_same_fam"] * p_val + cf2["log_other:prop_same_fam"] * x_oth * p_val
  lines(x_oth, y_hat, col = cols[k], lwd = 2)
}
legend("topright", legend = paste0("Stejna celed: ", props * 100, " %"),
       col = cols, lwd = 2, bty = "n")

p_grid <- seq(0, 1, length.out = 100)
b_marg <- cf2["log_other"] + cf2["log_other:prop_same_fam"] * p_grid
plot(p_grid * 100, b_marg, type = "l", lwd = 2, col = "navy",
     xlab = "Podil kvetu ze stejne celedi v okoli (%)",
     ylab = "Sklon vlivu okolnich kvetu (b_other)",
     main = "B) Prechod ke kompetici pri vyssim podilu stejne celedi")
abline(h = 0, lty = 2)
dev.off()

# ============================================================================
# CAST 7: HYPOTEZA H2 - PODIL STEJNE vs. JINE BARVY KVETU
# ============================================================================
cat("[7/8] Pocitam modely H2 pro barvu kvetu (color human)...\n")

m_C1_add   <- lmer(log_rate ~ log_own + log_other + prop_same_col +
                              (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)
m_C2_int   <- lmer(log_rate ~ log_own + log_other * prop_same_col +
                              (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)
m_C3_split <- lmer(log_rate ~ log_own + log_same_col + log_diff_col +
                              (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)
m_C4_both  <- lmer(log_rate ~ log_own + log_other * prop_same_col + log_other * prop_same_fam +
                              (1 | week_f) + (1 | Site) + (1 | Plant), data = H1_dat, REML = TRUE)

tab_col <- rbind(
  extract_coefs(m_C1_add,   "C1_aditivni_prop_same_col"),
  extract_coefs(m_C2_int,   "C2_interakce_log_other_x_prop_same_col"),
  extract_coefs(m_C3_split, "C3_rozdeleni_same_vs_diff_col"),
  extract_coefs(m_C4_both,  "C4_spolecny_model_barva_vs_celed")
)
write.csv(tab_col, file.path(OUT, "11_H2_podil_stejne_vs_jine_barvy.csv"), row.names = FALSE)

# ============================================================================
# CAST 8: HYPOTEZA H2 - 4 KOMBINACE (STEJNA/JINA CELED x STEJNA/JINA BARVA)
# ============================================================================
cat("[8/8] Pocitam modely H2 pro 4 kombinace celed x barva...\n")

m_F1 <- lmer(log_rate ~ log_own + log_other +
                        prop_same_fam_same_col + prop_same_fam_diff_col + prop_diff_fam_same_col +
                        (1 | week_f) + (1 | Site) + (1 | Plant),
             data = H1_dat, REML = TRUE)

m_F2 <- lmer(log_rate ~ log_own +
                        log_other * (prop_same_fam_same_col + prop_same_fam_diff_col + prop_diff_fam_same_col) +
                        (1 | week_f) + (1 | Site) + (1 | Plant),
             data = H1_dat, REML = TRUE)

tab_4 <- rbind(
  extract_coefs(m_F1, "F1_aditivni_4_skupiny"),
  extract_coefs(m_F2, "F2_interakce_4_skupiny")
)
write.csv(tab_4, file.path(OUT, "13_H2_kombinace_celed_x_barva.csv"), row.names = FALSE)

cat("\n====================================================================\n")
cat("HOTOVO! Vsech 8 casti probehlo uspesne.\n")
cat("Vsechny tabulky a grafy najdes ve slozce:\n  ", OUT, "\n")
cat("====================================================================\n")