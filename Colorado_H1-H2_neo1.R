# ============================================================================
# COLORADO: ANALYZA H1 A H2 V CISTEM BASE R + lme4::lmer (LOG-LOG SKALA)
# ============================================================================
# Vstupy:
#   1) Colorado_pollinators.xlsx     (list "List 1")
#   2) Colorado_phenology_check.xlsx (list "Daily_phenology")
# ============================================================================

options(stringsAsFactors = FALSE)

# 1. NASTAVENI CEST A BALICKU (pouze readxl a lme4)
# ----------------------------------------------------------------------------
BASE <- "C:/Users/Tadea/Documents/colorado"
if (!dir.exists(BASE)) {
  message("Vyber slozku se vstupnimi Excel soubory.")
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

POLL_FILE <- file.path(BASE, "Colorado_pollinators.xlsx")
PHEN_FILE <- file.path(BASE, "Colorado_phenology_check.xlsx")

# Pomocne funkce v cistem base R
clean_str <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  return(x)
}

is_valid_name <- function(x) {
  bad <- c("", "na", "n/a", "nic", "none", "unknown", "nezname")
  return(!(tolower(trimws(x)) %in% bad))
}

parse_excel_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  num <- suppressWarnings(as.numeric(as.character(x)))
  d <- as.Date(num, origin = "1899-12-30")
  miss <- is.na(d)
  if (any(miss)) d[miss] <- suppressWarnings(as.Date(as.character(x[miss])))
  return(d)
}

parse_time_min <- function(x) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    val <- x[i]
    if (is.na(val)) next
    if (inherits(val, "POSIXt")) {
      lt <- as.POSIXlt(val)
      out[i] <- lt$hour * 60 + lt$min + lt$sec / 60
      next
    }
    s <- trimws(as.character(val))
    num <- suppressWarnings(as.numeric(gsub(",", ".", s, fixed = TRUE)))
    if (is.finite(num)) {
      out[i] <- if (num < 1) num * 1440 else if (num >= 24) (num - floor(num)) * 1440 else num * 60
      next
    }
    m <- regmatches(s, regexpr("[0-2]?[0-9]:[0-5][0-9](:[0-5][0-9])?", s))
    if (length(m) == 1 && nzchar(m)) {
      parts <- as.numeric(strsplit(m, ":", fixed = TRUE)[[1]])
      out[i] <- parts[1] * 60 + parts[2] + ifelse(length(parts) > 2, parts[3] / 60, 0)
    }
  }
  return(out)
}

get_iso_week <- function(d) as.integer(format(d, "%V"))
get_iso_year <- function(d) as.integer(format(d, "%G"))

# Test alometrickeho exponentu b proti libovolne nulove hypoteze (napr. b0 = 1)
test_slope <- function(model, term, b0 = 1) {
  cf <- summary(model)$coefficients
  if (!(term %in% rownames(cf))) stop("Promenna nenalezena v modelu: ", term)
  est <- cf[term, "Estimate"]
  se  <- cf[term, "Std. Error"]
  t_val <- (est - b0) / se
  df_approx <- nrow(model@frame) - nrow(cf)
  p_val <- 2 * pt(abs(t_val), df = df_approx, lower.tail = FALSE)
  ci_low  <- est - qt(0.975, df = df_approx) * se
  ci_high <- est + qt(0.975, df = df_approx) * se
  data.frame(
    term     = term,
    H0_slope = b0,
    estimate = round(est, 4),
    SE       = round(se, 4),
    CI_95_lo = round(ci_low, 4),
    CI_95_hi = round(ci_high, 4),
    t_value  = round(t_val, 3),
    df       = df_approx,
    p_value  = signif(p_val, 4)
  )
}

# ============================================================================
# 2. NACTENI NAVSTEV OPYLOVACU A VYPOCET TERENNIHO USILI (POCET PRUCHODU)
# ============================================================================
hdr <- names(read_excel(POLL_FILE, sheet = "List 1", n_max = 0))
col_t <- rep("guess", length(hdr))
col_t[hdr %in% c("od", "do")] <- "text"
raw_poll <- as.data.frame(read_excel(POLL_FILE, sheet = "List 1", col_types = col_t))

P <- data.frame(
  Date   = parse_excel_date(raw_poll$Datum),
  Site   = clean_str(raw_poll$Site),
  run    = clean_str(raw_poll$run),
  side   = toupper(clean_str(raw_poll[["R/L"]])),
  Plant  = clean_str(raw_poll$Plant),
  Order  = clean_str(raw_poll$order),
  Group  = clean_str(raw_poll$Group1),
  visits = suppressWarnings(as.numeric(raw_poll$quantity)),
  start  = parse_time_min(raw_poll$od),
  end    = parse_time_min(raw_poll$do)
)

# Ponechame jen radky s platnym datem, lokalitou, pruchodem a stranou R/L
P <- subset(P, !is.na(Date) & Site != "" & run != "" & side %in% c("R", "L"))
P$end <- ifelse(is.finite(P$end) & is.finite(P$start) & P$end < P$start, P$end + 1440, P$end)
P$year  <- get_iso_year(P$Date)
P$week  <- get_iso_week(P$Date)
P$block <- paste(P$Site, P$year, P$week, sep = "__")

# Kontrola kompletnosti jednotlivych pruchodu (musi byt pritomny obe strany R i L)
P$run_id <- paste(P$Site, P$Date, P$run, sep = "__")
side_counts <- tapply(P$side, P$run_id, function(s) length(unique(s)))

min_start <- tapply(P$start, P$run_id, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
max_end   <- tapply(P$end,   P$run_id, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
run_dur   <- max_end - min_start

valid_runs <- names(side_counts)[
  side_counts == 2 & is.finite(run_dur) & run_dur > 0 & run_dur <= 480
]

# Pocet kompletnich pruchodu v kazdem bloku (Site x rok x tyden)
runs_df <- unique(P[P$run_id %in% valid_runs, c("Site", "year", "week", "block", "run_id")])
E_week <- aggregate(run_id ~ Site + year + week + block, data = runs_df, FUN = length)
names(E_week)[names(E_week) == "run_id"] <- "n_runs"

# Navstevy pouze v kompletnich pruchodech a na platnych rostlinach
P_valid <- subset(P, run_id %in% valid_runs & is_valid_name(Plant) &
                     is.finite(visits) & visits > 0)

V_block <- aggregate(visits ~ Site + year + week + block + Plant, data = P_valid, FUN = sum)

# ============================================================================
# 3. NACTENI FENOLOGIE (ABUNDANCE KVETU) A SPOJENI S NAVSTEVAMI
# ============================================================================
raw_phen <- as.data.frame(read_excel(PHEN_FILE, sheet = "Daily_phenology"))
F_dat <- data.frame(
  Site    = clean_str(raw_phen$Site),
  Plant   = clean_str(raw_phen$Plant),
  Date    = parse_excel_date(raw_phen$Datum),
  flowers = suppressWarnings(as.numeric(raw_phen$abundance_interpolated))
)
F_dat <- subset(F_dat, !is.na(Date) & Site != "" & is_valid_name(Plant) &
                       is.finite(flowers) & flowers >= 0)
F_dat$year  <- get_iso_year(F_dat$Date)
F_dat$week  <- get_iso_week(F_dat$Date)
F_dat$block <- paste(F_dat$Site, F_dat$year, F_dat$week, sep = "__")

# Sezonni maximum kvetu pro kazdy druh na dane lokalite
peaks <- aggregate(flowers ~ Site + Plant + year, data = F_dat, FUN = max)
names(peaks)[names(peaks) == "flowers"] <- "peak_flowers"

# Tydenni prumer kvetu v tydnech, kdy probihalo sbirani opylovacu
F_week <- aggregate(flowers ~ Site + Plant + year + week + block, data = F_dat, FUN = mean)
F_week <- merge(F_week, E_week, by = c("Site", "year", "week", "block"))
F_week <- merge(F_week, peaks,  by = c("Site", "Plant", "year"))
F_week <- subset(F_week, flowers > 0 & n_runs > 0)

F_week$relative_peak <- F_week$flowers / F_week$peak_flowers

# Celkovy pocet kvetu vsech OSTATNICH druhu na lokalite v danem tydnu
block_totals <- ave(F_week$flowers, F_week$block, FUN = sum)
F_week$other_flowers <- block_totals - F_week$flowers

# Spojeni s poctem navstev (chybejici navstevy = 0)
H1_dat <- merge(F_week, V_block[, c("block", "Plant", "visits")],
                by = c("block", "Plant"), all.x = TRUE)
H1_dat$visits[is.na(H1_dat$visits)] <- 0

# ============================================================================
# 4. LOG-LOG TRANSFORMACE A FITOVANI MODELU PRES lme4::lmer
# ============================================================================
# OPRAVA 1: log((visits / n_runs) + 1) zajisti, ze nulove navstevy maji vzdy
#           hodnotu 0 bez ohledu na pocet pruchodu n_runs.
# OPRAVA 2: Pouzivame fixni efekty Site + week_f (misto block), aby nebyl
#           mezidruhovy prediktor log_other matematicky svazany s konstantnim
#           souctem kvetu v bloku.
H1_dat$visits_per_run <- H1_dat$visits / H1_dat$n_runs
H1_dat$log_rate       <- log(H1_dat$visits_per_run + 1)
H1_dat$log_own        <- log(H1_dat$flowers)
H1_dat$log_other      <- log1p(H1_dat$other_flowers)
H1_dat$Plant          <- factor(H1_dat$Plant)
H1_dat$Site           <- factor(H1_dat$Site)
H1_dat$week_f         <- factor(H1_dat$week)
H1_dat$block          <- factor(H1_dat$block)

write.csv(H1_dat, file.path(OUT, "01_H1_data_pro_lmer.csv"), row.names = FALSE)

# Model 1: Vnitrodruhovy efekt (vlastni kvety) + kontrola fenologie, lokality a tydne
m1_own <- lmer(log_rate ~ log_own + relative_peak + Site + week_f + (1 | Plant),
               data = H1_dat, REML = TRUE)

# Model 2: Vnitrodruhovy + mezidruhovy efekt (ostatni kvety po odfiltrovani vlastnich)
m2_own_other <- lmer(log_rate ~ log_own + log_other + relative_peak + Site + week_f + (1 | Plant),
                     data = H1_dat, REML = TRUE)

cat("\n====================================================================\n")
cat("SHRNUTI MODELU H1 (lmer v log-log skale)\n")
cat("====================================================================\n")
print(summary(m2_own_other)$coefficients[c("(Intercept)", "log_own", "log_other", "relative_peak"), ])

# ============================================================================
# 5. TESTOVANI ALOMETRIE vs. IZOMETRIE (b = 1)
# ============================================================================
# Izometrie: b = 1 (navstevnost na 1 kvet nezavisi na celkovem poctu kvetu)
# Alometrie: b != 1 (b < 1 kompetice/saturace, b > 1 facilitace)
allometry_tests <- rbind(
  cbind(model = "m1_own (Site + week)",       test_slope(m1_own,       "log_own",   b0 = 1)),
  cbind(model = "m2_own_other (Site + week)", test_slope(m2_own_other, "log_own",   b0 = 1)),
  cbind(model = "m2_own_other (Site + week)", test_slope(m2_own_other, "log_other", b0 = 0)),
  cbind(model = "m2_own_other (Site + week)", test_slope(m2_own_other, "log_other", b0 = 1))
)

cat("\n====================================================================\n")
cat("TESTY ALOMETRICKÉHO EXPONENTU b (H0: b = 1 pro vlastni, b = 0/1 pro okolni)\n")
cat("====================================================================\n")
print(allometry_tests, row.names = FALSE)
write.csv(allometry_tests, file.path(OUT, "02_H1_testy_alometrie_b.csv"), row.names = FALSE)

# ============================================================================
# 6. H2: TAXONOMICKA PODOBNOST SOUSEDU (ORDER vs. FINE vs. SIM_GAP) V lmer
# ============================================================================
P_tax <- subset(P_valid, is_valid_name(Order) & is_valid_name(Group))
P_tax$fine_group <- paste(P_tax$Order, P_tax$Group, sep = "__")

calc_similarity_vector <- function(df_all, df_tax, level = c("Order", "fine_group"), min_visits = 5) {
  level <- match.arg(level)
  sim_out <- rep(NA_real_, nrow(df_all))
  cov_out <- rep(NA_real_, nrow(df_all))

  for (wk in unique(df_all$week)) {
    tr <- df_tax[df_tax$week != wk, ]
    if (nrow(tr) == 0) next
    tab <- tapply(tr$visits, list(tr$Plant, tr[[level]]), sum, default = 0)
    tot <- rowSums(tab)
    prop_mat <- sweep(tab, 1, pmax(tot, 1), "/")

    idx_wk <- which(df_all$week == wk)
    for (blk in unique(df_all$block[idx_wk])) {
      rows <- idx_wk[df_all$block[idx_wk] == blk]
      if (length(rows) < 2) next
      for (r in rows) {
        p_foc <- as.character(df_all$Plant[r])
        if (!(p_foc %in% rownames(prop_mat)) || tot[p_foc] < min_visits) next

        others <- rows[rows != r & df_all$flowers[rows] > 0]
        if (length(others) == 0) next
        w <- df_all$flowers[others]
        p_oth <- as.character(df_all$Plant[others])

        known <- (p_oth %in% rownames(prop_mat)) & (tot[p_oth] >= min_visits)
        cov_out[r] <- sum(w[known]) / sum(w)
        if (!any(known)) next

        overlaps <- sapply(p_oth[known], function(po) {
          sum(pmin(prop_mat[p_foc, ], prop_mat[po, ]))
        })
        sim_out[r] <- weighted.mean(overlaps, w[known])
      }
    }
  }
  return(list(sim = sim_out, cov = cov_out))
}

res_order <- calc_similarity_vector(H1_dat, P_tax, level = "Order",      min_visits = 5)
res_fine  <- calc_similarity_vector(H1_dat, P_tax, level = "fine_group", min_visits = 5)

H1_dat$sim_order <- res_order$sim
H1_dat$sim_fine  <- res_fine$sim
H1_dat$sim_gap   <- H1_dat$sim_order - H1_dat$sim_fine
H1_dat$cov_tax   <- res_fine$cov

H2_dat <- subset(H1_dat, is.finite(sim_order) & is.finite(sim_fine) & cov_tax >= 0.50)
write.csv(H2_dat, file.path(OUT, "03_H2_data_pro_lmer.csv"), row.names = FALSE)

m_H2_order <- lmer(log_rate ~ sim_order + log_own + log_other + relative_peak + cov_tax +
                              Site + week_f + (1 | Plant),
                   data = H2_dat, REML = FALSE)

m_H2_fine  <- lmer(log_rate ~ sim_fine  + log_own + log_other + relative_peak + cov_tax +
                              Site + week_f + (1 | Plant),
                   data = H2_dat, REML = FALSE)

m_H2_gap   <- lmer(log_rate ~ sim_fine + sim_gap + log_own + log_other + relative_peak + cov_tax +
                              Site + week_f + (1 | Plant),
                   data = H2_dat, REML = FALSE)

cat("\n====================================================================\n")
cat("SROVNANI MODELU H2 (ANOVA / LRT V lmer)\n")
cat("====================================================================\n")
print(anova(m_H2_fine, m_H2_gap))

h2_coefs <- rbind(
  cbind(model = "m_H2_order", as.data.frame(summary(m_H2_order)$coefficients)),
  cbind(model = "m_H2_fine",  as.data.frame(summary(m_H2_fine)$coefficients)),
  cbind(model = "m_H2_gap",   as.data.frame(summary(m_H2_gap)$coefficients))
)
h2_coefs$term <- rownames(h2_coefs)
write.csv(h2_coefs, file.path(OUT, "04_H2_koeficienty_lmer.csv"), row.names = FALSE)

# ============================================================================
# 7. PREHLEDNY GRAF ALOMETRICKÉHO VZTAHU V LOG-LOG SKALE (BASE R GRAPHICS)
# ============================================================================
pdf(file.path(OUT, "05_graf_alometrie_log_log.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Panel A: Vnitrodruhovy vztah (log_own vs log_rate)
plot(H1_dat$log_own, H1_dat$log_rate,
     pch = 16, col = rgb(0.2, 0.4, 0.7, 0.4),
     xlab = "log(pocet vlastnich kvetu)",
     ylab = "log((navstevy / pocet pruchodu) + 1)",
     main = "A) Vnitrodruhovy vztah (y = e^a * x^b)")
fit_a <- lm(log_rate ~ log_own, data = H1_dat)
abline(fit_a, col = "firebrick", lwd = 2)
abline(a = mean(H1_dat$log_rate) - 1 * mean(H1_dat$log_own), b = 1,
       col = "black", lty = 2, lwd = 1.5)
legend("topleft", legend = c(sprintf("Odhad sklonu b = %.2f", coef(fit_a)[2]), "Izometrie (b = 1)"),
       col = c("firebrick", "black"), lty = c(1, 2), lwd = c(2, 1.5), bty = "n")

# Panel B: Mezidruhovy vztah (log_other vs rezidua po odfiltrovani vlastnich kvetu)
res_own <- resid(lm(log_rate ~ log_own + relative_peak, data = H1_dat))
plot(H1_dat$log_other, res_own,
     pch = 16, col = rgb(0.2, 0.6, 0.4, 0.4),
     xlab = "log(pocet kvetu ostatních druhu + 1)",
     ylab = "Rezidua log-navstevnosti (po odfiltrovani vlastnich kvetu)",
     main = "B) Mezidruhovy vztah (okolni kvety)")
abline(lm(res_own ~ H1_dat$log_other), col = "darkgreen", lwd = 2)
abline(h = 0, col = "black", lty = 2)

dev.off()
message("Hotovo! Vsechny vysledky ulozeny do: ", OUT)
