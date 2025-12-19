# ============================================
# Lineament extraction run template (PowerShell)
# Matches: lineaments_from_strikes_3.py + hp_bp_coincidence.py
# ============================================

$ErrorActionPreference = "Stop"

# -----------------------------
# User settings (EDIT THESE)
# -----------------------------
$PYTHON_BIN = "python"

# Script paths (edit to full paths if needed)
$LINEAMENT_SCRIPT  = "lineaments_from_strikes_3.py"
$COINCIDENCE_SCRIPT = "hp_bp_coincidence.py"

# Input files (run twice: HP and BP strike picks)
$HP_INPUT = "Test Data/hp_input_points.csv"
$BP_INPUT = "Test Data/bp_input_points.csv"

# Output directory (timestamped)
$RUN_TAG = Get-Date -Format "yyyyMMdd_HHmmss"
$OUT_DIR = "outputs\$RUN_TAG"
New-Item -ItemType Directory -Force -Path $OUT_DIR | Out-Null

# -----------------------------
# Input format options
# -----------------------------
# "csv" or "xyz"
$FORMAT = "csv"

# CSV column names (used when FORMAT="csv")
$COL_LINE   = "LineNumber"
$COL_E      = "Easting"
$COL_N      = "Northing"
$COL_STRIKE = "Strike"

# XYZ indices (used when FORMAT="xyz")
$IDX_LINE   = 0
$IDX_E      = 1
$IDX_N      = 2
$IDX_STRIKE = 3
$XYZ_HAS_HEADER = $false  # add --xyz_has_header if true

# -----------------------------
# Lineament algorithm parameters
# -----------------------------
$CONE_HALF_ANGLE_DEG   = 10
$STRIKE_MATCH_TOL_DEG  = 10
$MAX_LINK_DISTANCE     = 300   # <-- YOU SET THIS (map units, typically metres)
$CELL_SIZE             = 300   # typically set = MAX_LINK_DISTANCE
$MIN_POINTS            = 5

# Controls
$REQUIRE_RECIPROCAL = $false   # adds --require_reciprocal if true
$ALLOW_BACKWARD     = $false   # adds --allow_backward if true (treat strike as bidirectional axis-line)

# Skip-one-line bridging
$SKIP_LINES           = 1      # 0=off; 1 recommended
$SKIP_DISTANCE_FACTOR = 1.6
$SKIP_ANGLE_BONUS_DEG = 2.0

# Debug candidate output (set empty string to disable)
$DEBUG_OUT_HP = "$OUT_DIR\hp_debug_candidates.csv"
$DEBUG_OUT_BP = "$OUT_DIR\bp_debug_candidates.csv"

# -----------------------------
# HP/BP coincidence parameters
# -----------------------------
$MATCH_DIST       = 1500
$MATCH_STRIKE_TOL = 15
$MATCH_GRID_CELL  = 2000

# Output names for coincidence tool
$OUT_MATCHES       = "$OUT_DIR\hp_bp_matches.csv"
$OUT_CLASSIFIED    = "$OUT_DIR\hp_bp_classified.csv"
$OUT_INTERSECTIONS = "$OUT_DIR\hp_bp_intersections.csv"

# -----------------------------
# Build argument arrays
# -----------------------------
$COMMON_ARGS = @(
  "--format", $FORMAT,

  "--cone_half_angle_deg", $CONE_HALF_ANGLE_DEG,
  "--strike_match_tol_deg", $STRIKE_MATCH_TOL_DEG,
  "--max_link_distance", $MAX_LINK_DISTANCE,
  "--cell_size", $CELL_SIZE,

  "--min_points", $MIN_POINTS,

  "--skip_lines", $SKIP_LINES,
  "--skip_distance_factor", $SKIP_DISTANCE_FACTOR,
  "--skip_angle_bonus_deg", $SKIP_ANGLE_BONUS_DEG
)

# Input format specifics
$FORMAT_ARGS = @()
if ($FORMAT -eq "csv") {
  $FORMAT_ARGS = @(
    "--col_line", $COL_LINE,
    "--col_easting", $COL_E,
    "--col_northing", $COL_N,
    "--col_strike", $COL_STRIKE
  )
}
elseif ($FORMAT -eq "xyz") {
  $FORMAT_ARGS = @(
    "--idx_line", $IDX_LINE,
    "--idx_easting", $IDX_E,
    "--idx_northing", $IDX_N,
    "--idx_strike", $IDX_STRIKE
  )
  if ($XYZ_HAS_HEADER) {
    $FORMAT_ARGS += "--xyz_has_header"
  }
}
else {
  throw "FORMAT must be 'csv' or 'xyz'"
}

# Optional boolean flags
$BOOL_ARGS = @()
if ($REQUIRE_RECIPROCAL) { $BOOL_ARGS += "--require_reciprocal" }
if ($ALLOW_BACKWARD)     { $BOOL_ARGS += "--allow_backward" }

# Debug args (optional)
$DEBUG_ARGS_HP = @()
if ($DEBUG_OUT_HP -ne "") { $DEBUG_ARGS_HP = @("--out_debug_candidates", $DEBUG_OUT_HP) }

$DEBUG_ARGS_BP = @()
if ($DEBUG_OUT_BP -ne "") { $DEBUG_ARGS_BP = @("--out_debug_candidates", $DEBUG_OUT_BP) }

# -----------------------------
# Run HP
# -----------------------------
Write-Host "=== Running HP lineaments ==="

& $PYTHON_BIN $LINEAMENT_SCRIPT $HP_INPUT `
  $COMMON_ARGS `
  $FORMAT_ARGS `
  $BOOL_ARGS `
  --out_points   "$OUT_DIR\hp_lineaments_points.csv" `
  --out_segments "$OUT_DIR\hp_lineaments_segments.csv" `
  --out_summary  "$OUT_DIR\hp_lineaments_summary.csv" `
  $DEBUG_ARGS_HP

# -----------------------------
# Run BP
# -----------------------------
Write-Host "=== Running BP lineaments ==="

& $PYTHON_BIN $LINEAMENT_SCRIPT $BP_INPUT `
  $COMMON_ARGS `
  $FORMAT_ARGS `
  $BOOL_ARGS `
  --out_points   "$OUT_DIR\bp_lineaments_points.csv" `
  --out_segments "$OUT_DIR\bp_lineaments_segments.csv" `
  --out_summary  "$OUT_DIR\bp_lineaments_summary.csv" `
  $DEBUG_ARGS_BP

# -----------------------------
# Run coincidence + intersections
# -----------------------------
Write-Host "=== Running HP/BP coincidence + intersections ==="

& $PYTHON_BIN $COINCIDENCE_SCRIPT `
  --hp_summary   "$OUT_DIR\hp_lineaments_summary.csv" `
  --hp_segments  "$OUT_DIR\hp_lineaments_segments.csv" `
  --bp_summary   "$OUT_DIR\bp_lineaments_summary.csv" `
  --bp_segments  "$OUT_DIR\bp_lineaments_segments.csv" `
  --match_dist   $MATCH_DIST `
  --strike_tol   $MATCH_STRIKE_TOL `
  --cell         $MATCH_GRID_CELL `
  --out_matches       $OUT_MATCHES `
  --out_classified    $OUT_CLASSIFIED `
  --out_intersections $OUT_INTERSECTIONS

Write-Host "=== DONE ==="
Write-Host "Outputs: $OUT_DIR"
