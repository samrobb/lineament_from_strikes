#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User settings (edit these)
# -----------------------------

# Python executable (set to the one you use with QGIS / your environment)
PYTHON_BIN="python"

# Paths to scripts
LINEAMENT_SCRIPT="lineaments_from_strikes_3.py"
COINCIDENCE_SCRIPT="hp_bp_coincidence.py"

# Input files (run twice: HP and BP strike picks)
HP_INPUT="Test Data/hp_input_points.csv"
BP_INPUT="Test Data/bp_input_points.csv"

# Output folder (auto timestamped)
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="outputs/${RUN_TAG}"
mkdir -p "${OUT_DIR}"

# -----------------------------
# Input format options
# -----------------------------
# Choose FORMAT: "csv" or "xyz"
FORMAT="csv"

# If CSV: specify column names
COL_LINE="LineNumber"
COL_E="Easting"
COL_N="Northing"
COL_STRIKE="Strike"

# If XYZ: specify column indices (0-based) and whether it has a header
XYZ_HAS_HEADER="false"   # true/false
IDX_LINE=0
IDX_E=1
IDX_N=2
IDX_STRIKE=3

# -----------------------------
# Lineament algorithm parameters
# -----------------------------
CONE_HALF_ANGLE_DEG=15
STRIKE_MATCH_TOL_DEG=15

# Link/search geometry (tune these for your survey)
LINE_SPACING_M=200              # for your survey: 100-200 typical
LINK_FACTOR=2.0                 # max_link_distance = line_spacing * link_factor
MAX_LINK_DISTANCE=""            # optional override (leave "" to auto-derive)
CELL_SIZE=""                    # optional override (leave "" to auto-derive)
MIN_POINTS=4

# Controls
REQUIRE_RECIPROCAL="false"      # true/false (stricter)
ALLOW_BACKWARD="false"          # true/false (strike axis: often leave false; set true if needed)

# Skip-one-line bridging
SKIP_LINES=1                    # 0 to disable, 1 recommended
SKIP_DISTANCE_FACTOR=1.6
SKIP_ANGLE_BONUS_DEG=2.0

# Debug candidates output (empty to disable)
DEBUG_OUT_HP="${OUT_DIR}/hp_debug_candidates.csv"
DEBUG_OUT_BP="${OUT_DIR}/bp_debug_candidates.csv"

# -----------------------------
# HP/BP coincidence parameters
# -----------------------------
MATCH_DIST=1500                 # meters (center-to-center match threshold)
MATCH_STRIKE_TOL=15             # degrees (axis strike diff tolerance)
MATCH_GRID_CELL=2000            # for speed (meters)

# -----------------------------
# Helper: build common args
# -----------------------------
COMMON_ARGS=(
  "--cone_half_angle_deg" "${CONE_HALF_ANGLE_DEG}"
  "--strike_match_tol_deg" "${STRIKE_MATCH_TOL_DEG}"
  "--line_spacing" "${LINE_SPACING_M}"
  "--link_factor" "${LINK_FACTOR}"
  "--min_points" "${MIN_POINTS}"
  "--skip_lines" "${SKIP_LINES}"
  "--skip_distance_factor" "${SKIP_DISTANCE_FACTOR}"
  "--skip_angle_bonus_deg" "${SKIP_ANGLE_BONUS_DEG}"
)

# Optional boolean flags
BOOL_ARGS=()
if [[ "${REQUIRE_RECIPROCAL}" == "true" ]]; then
  BOOL_ARGS+=("--require_reciprocal")
fi
if [[ "${ALLOW_BACKWARD}" == "true" ]]; then
  BOOL_ARGS+=("--allow_backward")
fi

# Optional overrides
OVERRIDE_ARGS=()
if [[ -n "${MAX_LINK_DISTANCE}" ]]; then
  OVERRIDE_ARGS+=("--max_link_distance" "${MAX_LINK_DISTANCE}")
fi
if [[ -n "${CELL_SIZE}" ]]; then
  OVERRIDE_ARGS+=("--cell_size" "${CELL_SIZE}")
fi

# Input format args
FORMAT_ARGS=("--format" "${FORMAT}")
if [[ "${FORMAT}" == "csv" ]]; then
  FORMAT_ARGS+=(
    "--col_line" "${COL_LINE}"
    "--col_easting" "${COL_E}"
    "--col_northing" "${COL_N}"
    "--col_strike" "${COL_STRIKE}"
  )
elif [[ "${FORMAT}" == "xyz" ]]; then
  FORMAT_ARGS+=(
    "--idx_line" "${IDX_LINE}"
    "--idx_easting" "${IDX_E}"
    "--idx_northing" "${IDX_N}"
    "--idx_strike" "${IDX_STRIKE}"
  )
  if [[ "${XYZ_HAS_HEADER}" == "true" ]]; then
    FORMAT_ARGS+=("--xyz_has_header")
  fi
else
  echo "ERROR: FORMAT must be 'csv' or 'xyz'"
  exit 1
fi

# -----------------------------
# Run HP extraction
# -----------------------------
echo "=== Running HP lineaments ==="
"${PYTHON_BIN}" "${LINEAMENT_SCRIPT}" "${HP_INPUT}" \
  "${FORMAT_ARGS[@]}" \
  "${COMMON_ARGS[@]}" \
  "${OVERRIDE_ARGS[@]}" \
  "${BOOL_ARGS[@]}" \
  --out_points "${OUT_DIR}/hp_lineaments_points.csv" \
  --out_segments "${OUT_DIR}/hp_lineaments_segments.csv" \
  --out_summary "${OUT_DIR}/hp_lineaments_summary.csv" \
  --out_debug_candidates "${DEBUG_OUT_HP}"

# -----------------------------
# Run BP extraction
# -----------------------------
echo "=== Running BP lineaments ==="
"${PYTHON_BIN}" "${LINEAMENT_SCRIPT}" "${BP_INPUT}" \
  "${FORMAT_ARGS[@]}" \
  "${COMMON_ARGS[@]}" \
  "${OVERRIDE_ARGS[@]}" \
  "${BOOL_ARGS[@]}" \
  --out_points "${OUT_DIR}/bp_lineaments_points.csv" \
  --out_segments "${OUT_DIR}/bp_lineaments_segments.csv" \
  --out_summary "${OUT_DIR}/bp_lineaments_summary.csv" \
  --out_debug_candidates "${DEBUG_OUT_BP}"

# -----------------------------
# Run coincidence + intersections
# -----------------------------
echo "=== Running HP/BP coincidence + intersections ==="
"${PYTHON_BIN}" "${COINCIDENCE_SCRIPT}" \
  --hp_summary "${OUT_DIR}/hp_lineaments_summary.csv" \
  --hp_segments "${OUT_DIR}/hp_lineaments_segments.csv" \
  --bp_summary "${OUT_DIR}/bp_lineaments_summary.csv" \
  --bp_segments "${OUT_DIR}/bp_lineaments_segments.csv" \
  --match_dist "${MATCH_DIST}" \
  --strike_tol "${MATCH_STRIKE_TOL}" \
  --cell "${MATCH_GRID_CELL}" \
  --out_prefix "${OUT_DIR}/hp_bp"

echo "=== Done ==="
echo "Outputs in: ${OUT_DIR}"
echo "HP summary: ${OUT_DIR}/hp_lineaments_summary.csv"
echo "BP summary: ${OUT_DIR}/bp_lineaments_summary.csv"
echo "HP/BP classified: ${OUT_DIR}/hp_bp_classified.csv"
echo "HP/BP intersections: ${OUT_DIR}/hp_bp_intersections.csv"
