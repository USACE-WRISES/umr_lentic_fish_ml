# ── Paths ──────────────────────────────────────────────────────────────────
DATA_PATH   <- "data/cleaned/fish_all_cv.RData"
OUTPUT_PATH <- "output/"

# ── Species ────────────────────────────────────────────────────────────────
# Populated automatically from data in 01_data_prep.R
# but can be subset here for testing
SP_SUBSET <- NULL  # NULL = run all species

# ── Variable sets ──────────────────────────────────────────────────────────
VARS_REST <- c("occ", "cluster", "period",
               "temp", "depth", "current", "do",
               "vegd", "substrt", "snag")

VARS_BLGL <- c("occ", "cluster",
               "temp", "depth", "current", "do")

CONTINUOUS_VARS  <- c("temp", "depth", "current", "do")
CATEGORICAL_VARS <- c("period", "vegd", "substrt", "snag")

# ── CV settings ────────────────────────────────────────────────────────────
N_FOLDS         <- 6       # pool-level spatial CV
N_CLUSTERS      <- 21      # patch clusters for within-pool CV
N_REPEATS       <- 30      # repeats for patch-cluster CV

# ── Model settings ─────────────────────────────────────────────────────────
GRID_SIZE  <- 25
N_TREES    <- 1000
IMPORTANCE <- "permutation"

# ── Threshold settings ─────────────────────────────────────────────────────
THRESH_METHOD <- "dist_to_01"   # "dist_to_01", "tss", or "sens_spec_diff"
THRESH_SEQ    <- seq(0, 1, by = 0.01)

# ── Model type labels ──────────────────────────────────────────────────────
# Used to name output files and tag results
MODEL_TYPES <- list(
  spatial_rest   = list(cv = "spatial",  vars = "rest"),
  spatial_blgl   = list(cv = "spatial",  vars = "blgl"),
  cluster_rest   = list(cv = "cluster",  vars = "rest"),
  cluster_blgl   = list(cv = "cluster",  vars = "blgl"),
  random_rest    = list(cv = "random",   vars = "rest"),
  random_blgl    = list(cv = "random",   vars = "blgl")
)