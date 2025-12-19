#!/usr/bin/env python3
"""
Lineament extraction from strike picks on neighbouring survey lines.

Core idea:
- Each point (E,N,Strike,Line) can link to at most one point on the next/prev line.
- Linking is allowed if the candidate lies inside a strike-directed cone (±cone_half_angle),
  within max_link_distance, and has similar strike (±strike_match_tol).
- Chains are grown forward and backward; only chains with >= min_points are kept.
- Optional skip-one-line bridging can connect across a single missing/poor line.
- Optional debug output can dump the candidate evaluation for tuning.

Outputs:
- lineaments_points.csv
- lineaments_segments.csv
- lineaments_summary.csv (center point + mean strike + length + quality metrics + confidence score)
"""

import math
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Tuple, Set, Optional, Iterable, Any


# ----------------------------
# Data model
# ----------------------------

@dataclass(frozen=True)
class Point:
    idx: int
    line: int
    e: float
    n: float
    strike: float  # degrees, axis in [0,180)


# ----------------------------
# Angle & geometry helpers
# ----------------------------

def normalize_axis_deg(a: float) -> float:
    """Normalize an axis direction to [0,180). Strike is an axis: 0 and 180 are the same."""
    a = a % 180.0
    if a < 0:
        a += 180.0
    return a


def axis_angle_diff_deg(a: float, b: float) -> float:
    """Smallest difference between two axis directions (both in [0,180))."""
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


def vector_axis_azimuth_deg(e1: float, n1: float, e2: float, n2: float) -> float:
    """
    Azimuth of the vector from p1->p2 as an axis direction [0,180).
    (0 ~ North, 90 ~ East)
    """
    de = e2 - e1
    dn = n2 - n1
    az = math.degrees(math.atan2(de, dn))  # atan2(x=east, y=north)
    az = az % 360.0
    if az >= 180.0:
        az -= 180.0
    return az


def dist(e1: float, n1: float, e2: float, n2: float) -> float:
    return math.hypot(e2 - e1, n2 - n1)


def ray_perp_distance(p: Point, strike_axis_deg: float, q: Point) -> Tuple[float, float]:
    """
    Distance metrics relative to a ray from p in strike direction (axis):
      - along: signed distance along the ray direction (must be >0 to be 'forward')
      - perp: perpendicular distance to the ray line (always positive)

    Uses a unit direction vector for strike_axis_deg interpreted as azimuth-from-North.
    """
    a = normalize_axis_deg(strike_axis_deg)
    rad = math.radians(a)
    ux = math.sin(rad)  # East component
    uy = math.cos(rad)  # North component

    dx = q.e - p.e
    dy = q.n - p.n

    along = dx * ux + dy * uy
    perp = abs(dx * (-uy) + dy * (ux))
    return along, perp


# ----------------------------
# Circular statistics on AXIS data (0-180)
# ----------------------------

def mean_axis_deg(values_deg: List[float]) -> float:
    """
    Mean of axis values in [0,180), using the 2-theta method.
    """
    if not values_deg:
        return float("nan")
    s2 = 0.0
    c2 = 0.0
    for a in values_deg:
        a = normalize_axis_deg(a)
        r = math.radians(2.0 * a)
        s2 += math.sin(r)
        c2 += math.cos(r)
    if abs(s2) < 1e-12 and abs(c2) < 1e-12:
        return float("nan")
    mean2 = math.degrees(math.atan2(s2, c2))
    return normalize_axis_deg(0.5 * mean2)


def axis_circular_std_deg(values_deg: List[float]) -> float:
    """
    Circular standard deviation for axis data (0-180), using 2-theta resultant length:
      R = |sum(exp(i*2θ))| / n
      sigma_2 = sqrt(-2 ln R)  (radians)   [for circular]
      sigma_axis = 0.5 * sigma_2 (converted back to degrees)
    """
    n = len(values_deg)
    if n < 2:
        return 0.0
    s2 = 0.0
    c2 = 0.0
    for a in values_deg:
        r = math.radians(2.0 * normalize_axis_deg(a))
        s2 += math.sin(r)
        c2 += math.cos(r)
    R = math.hypot(s2, c2) / max(1, n)
    # Guard against log(0)
    R = max(1e-12, min(1.0, R))
    sigma2 = math.sqrt(max(0.0, -2.0 * math.log(R)))  # radians
    return 0.5 * math.degrees(sigma2)


# ----------------------------
# Input parsing
# ----------------------------

def parse_line_number(s: str) -> int:
    """
    Extract integer from line labels like 'L12', 'L0012', 'l12', 'Line L12'.
    Falls back to plain int if the whole string is numeric.
    """
    s = str(s).strip()
    m = re.search(r'(\d+)', s)
    if m:
        return int(m.group(1))
    raise ValueError(f"Could not parse LineNumber from: {s!r}")


def read_points_csv(csv_path: str,
                    col_line: str = "LineNumber",
                    col_e: str = "Easting",
                    col_n: str = "Northing",
                    col_strike: str = "Strike") -> List[Point]:
    points: List[Point] = []
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {col_line, col_e, col_n, col_strike}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

        for i, row in enumerate(reader):
            line = parse_line_number(row[col_line])
            e = float(row[col_e])
            n = float(row[col_n])
            strike = normalize_axis_deg(float(row[col_strike]))
            points.append(Point(idx=i, line=line, e=e, n=n, strike=strike))
    return points


def read_points_xyz(path: str,
                    idx_line: int,
                    idx_e: int,
                    idx_n: int,
                    idx_strike: int,
                    has_header: bool = False) -> List[Point]:
    points: List[Point] = []
    with open(path, "r", encoding="utf-8-sig") as f:
        if has_header:
            next(f, None)
        for i, line in enumerate(f):
            if not line.strip():
                continue
            parts = line.split()
            line_num = parse_line_number(parts[idx_line])
            e = float(parts[idx_e])
            n = float(parts[idx_n])
            strike = normalize_axis_deg(float(parts[idx_strike]))
            points.append(Point(idx=len(points), line=line_num, e=e, n=n, strike=strike))
    return points


# ----------------------------
# Spatial indexing
# ----------------------------

def nearby_candidates(grid: Dict[Tuple[int, int], List[int]],
                      points: List[Point],
                      e: float, n: float,
                      radius: float,
                      cell_size: float) -> List[int]:
    cx = int(math.floor(e / cell_size))
    cy = int(math.floor(n / cell_size))
    r_cells = int(math.ceil(radius / cell_size))
    out: List[int] = []
    for dx in range(-r_cells, r_cells + 1):
        for dy in range(-r_cells, r_cells + 1):
            out.extend(grid.get((cx + dx, cy + dy), []))
    return out


# ----------------------------
# Candidate evaluation (debuggable)
# ----------------------------

def iter_candidates_on_line(points: List[Point],
                            grid: Dict[Tuple[int, int], List[int]],
                            p_idx: int,
                            target_line: int,
                            max_link_distance: float,
                            cell_size: float,
                            cone_half_angle_deg: float,
                            strike_match_tol_deg: float,
                            *,
                            enforce_forward_ray: bool = True) -> Iterable[Dict[str, Any]]:
    """
    Yield candidate evaluation dicts for every point on target_line in the nearby grid.

    Each yielded dict contains:
      - q_idx, q_e, q_n, q_strike
      - d, az_axis, ang_misfit, strike_misfit, along, perp
      - passed (bool) and fail_reason (str or '')
      - score tuple (used for ranking) when passed
    """
    p = points[p_idx]
    cand_idxs = nearby_candidates(grid, points, p.e, p.n, max_link_distance, cell_size)

    for j in cand_idxs:
        q = points[j]
        if q.line != target_line:
            continue

        d = dist(p.e, p.n, q.e, q.n)
        if d == 0:
            yield {"p_idx": p_idx, "q_idx": j, "passed": False, "fail_reason": "zero_distance"}
            continue
        if d > max_link_distance:
            yield {"p_idx": p_idx, "q_idx": j, "passed": False, "fail_reason": "too_far", "d": d}
            continue

        az_axis = vector_axis_azimuth_deg(p.e, p.n, q.e, q.n)
        ang_misfit = axis_angle_diff_deg(az_axis, p.strike)
        if ang_misfit > cone_half_angle_deg:
            yield {"p_idx": p_idx, "q_idx": j, "passed": False, "fail_reason": "outside_cone",
                   "d": d, "az_axis": az_axis, "ang_misfit": ang_misfit}
            continue

        strike_misfit = axis_angle_diff_deg(q.strike, p.strike)
        if strike_misfit > strike_match_tol_deg:
            yield {"p_idx": p_idx, "q_idx": j, "passed": False, "fail_reason": "strike_mismatch",
                   "d": d, "az_axis": az_axis, "ang_misfit": ang_misfit, "strike_misfit": strike_misfit}
            continue

        along, perp = ray_perp_distance(p, p.strike, q)
        if enforce_forward_ray and along <= 0:
            yield {"p_idx": p_idx, "q_idx": j, "passed": False, "fail_reason": "behind_ray",
                   "d": d, "az_axis": az_axis, "ang_misfit": ang_misfit, "strike_misfit": strike_misfit,
                   "along": along, "perp": perp}
            continue

        score = (ang_misfit, perp, d, strike_misfit)
        yield {"p_idx": p_idx, "q_idx": j, "q_e": q.e, "q_n": q.n, "q_strike": q.strike,
               "d": d, "az_axis": az_axis, "ang_misfit": ang_misfit, "strike_misfit": strike_misfit,
               "along": along, "perp": perp,
               "passed": True, "fail_reason": "", "score": score}


def best_match_on_line(points: List[Point],
                       grid: Dict[Tuple[int, int], List[int]],
                       p_idx: int,
                       target_line: int,
                       max_link_distance: float,
                       cell_size: float,
                       cone_half_angle_deg: float,
                       strike_match_tol_deg: float,
                       *,
                       debug_rows: Optional[List[Dict[str, Any]]] = None,
                       debug_tag: str = "",
                       enforce_forward_ray: bool = True) -> Optional[int]:
    """
    Return the single best q on target_line, with optional debug row capture.
    Tie-break: smallest ang_misfit, then smallest perp, then smallest distance, then strike misfit.
    """
    best: Optional[Tuple[Tuple[float, float, float, float], int]] = None
    for rec in iter_candidates_on_line(
        points, grid, p_idx, target_line,
        max_link_distance, cell_size,
        cone_half_angle_deg, strike_match_tol_deg,
        enforce_forward_ray=enforce_forward_ray
    ):
        if debug_rows is not None:
            # decorate with context
            p = points[p_idx]
            rec2 = dict(rec)
            rec2.update({
                "tag": debug_tag,
                "p_line": p.line, "p_e": p.e, "p_n": p.n, "p_strike": p.strike,
                "target_line": target_line,
                "max_link_distance": max_link_distance,
                "cone_half_angle_deg": cone_half_angle_deg,
                "strike_match_tol_deg": strike_match_tol_deg,
            })
            debug_rows.append(rec2)

        if not rec.get("passed", False):
            continue
        score = rec["score"]
        j = rec["q_idx"]
        if best is None or score < best[0]:
            best = (score, j)

    return None if best is None else best[1]


# ----------------------------
# Lineament building (with optional skip-one-line bridging)
# ----------------------------

def build_lineaments(points: List[Point],
                     cone_half_angle_deg: float = 10.0,
                     strike_match_tol_deg: float = 10.0,
                     max_link_distance: float = 200.0,
                     min_points: int = 5,
                     cell_size: float = 200.0,
                     require_reciprocal: bool = True,
                     skip_lines: int = 0,
                     skip_distance_factor: float = 1.6,
                     skip_angle_bonus_deg: float = 2.0,
                     debug_rows: Optional[List[Dict[str, Any]]] = None,
                     enforce_forward_ray: bool = True
                     ) -> Tuple[List[List[int]], List[Tuple[int, int, int]]]:
    """
    Returns:
      - lineaments: list of chains (ordered point indices)
      - segments: list of (u,v,skipped) where skipped=0 normal, 1 bridged a missing line
    """
    by_line: Dict[int, List[int]] = defaultdict(list)
    for p in points:
        by_line[p.line].append(p.idx)

    sorted_lines = sorted(by_line.keys())
    if len(sorted_lines) < 2:
        return [], []

    next_line: Dict[int, Optional[int]] = {}
    prev_line: Dict[int, Optional[int]] = {}
    for i, ln in enumerate(sorted_lines):
        prev_line[ln] = sorted_lines[i - 1] if i > 0 else None
        next_line[ln] = sorted_lines[i + 1] if i < len(sorted_lines) - 1 else None

    # per-line grids
    grid_by_line: Dict[int, Dict[Tuple[int, int], List[int]]] = {}
    for ln, idxs in by_line.items():
        grid: Dict[Tuple[int, int], List[int]] = defaultdict(list)
        for idx in idxs:
            pp = points[idx]
            cx = int(math.floor(pp.e / cell_size))
            cy = int(math.floor(pp.n / cell_size))
            grid[(cx, cy)].append(pp.idx)
        grid_by_line[ln] = grid

    used: Set[int] = set()
    lineaments: List[List[int]] = []
    kept_segments: List[Tuple[int, int, int]] = []

    def try_link(p_idx: int,
                 target_ln: int,
                 *,
                 tag: str,
                 max_d: float,
                 cone: float,
                 strike_tol: float) -> Optional[int]:
        if target_ln not in grid_by_line:
            return None
        return best_match_on_line(
            points, grid_by_line[target_ln], p_idx, target_ln,
            max_d, cell_size, cone, strike_tol,
            debug_rows=debug_rows, debug_tag=tag,
            enforce_forward_ray=enforce_forward_ray
        )

    def link_forward(p_idx: int) -> Tuple[Optional[int], int]:
        """
        Returns (q_idx, skipped) where skipped=0 normal next line, skipped=1 bridged.
        """
        p = points[p_idx]
        nl = next_line[p.line]
        if nl is None:
            return None, 0

        q_idx = try_link(
            p_idx, nl, tag="fwd",
            max_d=max_link_distance,
            cone=cone_half_angle_deg,
            strike_tol=strike_match_tol_deg
        )
        if q_idx is None or q_idx in used:
            # try skip-one-line (or more, but default is 1)
            for _ in range(skip_lines):
                nnl = next_line.get(nl) if nl is not None else None
                if nnl is None:
                    break
                q2 = try_link(
                    p_idx, nnl, tag="skip_fwd",
                    max_d=max_link_distance * skip_distance_factor,
                    cone=cone_half_angle_deg + skip_angle_bonus_deg,
                    strike_tol=strike_match_tol_deg + skip_angle_bonus_deg
                )
                if q2 is not None and q2 not in used:
                    q_idx = q2
                    nl = nnl
                    return q_idx, 1
                nl = nnl
            return None, 0

        if require_reciprocal:
            pl = prev_line[points[q_idx].line]
            if pl is None:
                return None, 0
            back = try_link(
                q_idx, pl, tag="recip_back",
                max_d=max_link_distance,
                cone=cone_half_angle_deg,
                strike_tol=strike_match_tol_deg
            )
            if back != p_idx:
                return None, 0

        return q_idx, 0

    def link_backward(p_idx: int) -> Tuple[Optional[int], int]:
        p = points[p_idx]
        pl = prev_line[p.line]
        if pl is None:
            return None, 0

        q_idx = try_link(
            p_idx, pl, tag="bwd",
            max_d=max_link_distance,
            cone=cone_half_angle_deg,
            strike_tol=strike_match_tol_deg
        )
        if q_idx is None or q_idx in used:
            for _ in range(skip_lines):
                ppl = prev_line.get(pl) if pl is not None else None
                if ppl is None:
                    break
                q2 = try_link(
                    p_idx, ppl, tag="skip_bwd",
                    max_d=max_link_distance * skip_distance_factor,
                    cone=cone_half_angle_deg + skip_angle_bonus_deg,
                    strike_tol=strike_match_tol_deg + skip_angle_bonus_deg
                )
                if q2 is not None and q2 not in used:
                    q_idx = q2
                    pl = ppl
                    return q_idx, 1
                pl = ppl
            return None, 0

        if require_reciprocal:
            nl = next_line[points[q_idx].line]
            if nl is None:
                return None, 0
            fwd = try_link(
                q_idx, nl, tag="recip_fwd",
                max_d=max_link_distance,
                cone=cone_half_angle_deg,
                strike_tol=strike_match_tol_deg
            )
            if fwd != p_idx:
                return None, 0

        return q_idx, 0

    # seed order: by line then position
    all_idxs = sorted(range(len(points)), key=lambda i: (points[i].line, points[i].e, points[i].n))

    for seed in all_idxs:
        if seed in used:
            continue

        # backward chain
        bchain = [seed]
        bseg_skip: List[int] = []  # parallel to segments between bchain points (after reverse handled below)
        cur = seed
        while True:
            prev_idx, skipped = link_backward(cur)
            if prev_idx is None or prev_idx in bchain:
                break
            bchain.append(prev_idx)
            bseg_skip.append(skipped)
            cur = prev_idx
        # reverse bchain and skips to match order
        bchain = list(reversed(bchain))
        bseg_skip = list(reversed(bseg_skip))

        # forward chain
        fchain: List[int] = []
        fseg_skip: List[int] = []
        cur = seed
        while True:
            nxt, skipped = link_forward(cur)
            if nxt is None or nxt in fchain:
                break
            fchain.append(nxt)
            fseg_skip.append(skipped)
            cur = nxt

        # stitch
        chain = bchain + fchain
        seg_skips = bseg_skip + fseg_skip  # length = len(chain)-1

        # enforce one point per line (belt + suspenders)
        seen_lines: Set[int] = set()
        dedup_chain: List[int] = []
        dedup_skips: List[int] = []
        for i, idx in enumerate(chain):
            ln = points[idx].line
            if ln in seen_lines:
                continue
            seen_lines.add(ln)
            if dedup_chain:
                # keep the skip flag from the segment that led to this point
                # find corresponding segment index in original chain
                # (approx: use previous kept index position)
                # simplest: recompute later; here we just drop duplicates which are rare.
                pass
            dedup_chain.append(idx)

        # recompute skip flags for dedup_chain by checking line gaps
        dedup_skips = []
        for a, b in zip(dedup_chain[:-1], dedup_chain[1:]):
            la = points[a].line
            lb = points[b].line
            # if not adjacent in acquisition list, it was a skip bridge
            # (difference in *index* of line in sorted_lines)
            ia = sorted_lines.index(la)
            ib = sorted_lines.index(lb)
            dedup_skips.append(1 if abs(ib - ia) > 1 else 0)

        if len(dedup_chain) >= min_points:
            for idx in dedup_chain:
                used.add(idx)
            lineaments.append(dedup_chain)

            for (a, b), sk in zip(zip(dedup_chain[:-1], dedup_chain[1:]), dedup_skips):
                u, v = (a, b) if a < b else (b, a)
                kept_segments.append((u, v, sk))

    return lineaments, kept_segments


# ----------------------------
# Ordering + metrics + outputs
# ----------------------------

def order_by_line(points: List[Point], idxs: List[int]) -> List[int]:
    return sorted(idxs, key=lambda i: (points[i].line, points[i].e, points[i].n))


def polyline_length(points: List[Point], ordered_idxs: List[int]) -> float:
    if len(ordered_idxs) < 2:
        return 0.0
    total = 0.0
    for a, b in zip(ordered_idxs[:-1], ordered_idxs[1:]):
        pa, pb = points[a], points[b]
        total += dist(pa.e, pa.n, pb.e, pb.n)
    return total


def segment_axis_azimuths(points: List[Point], ordered_idxs: List[int]) -> List[float]:
    azs: List[float] = []
    for a, b in zip(ordered_idxs[:-1], ordered_idxs[1:]):
        pa, pb = points[a], points[b]
        azs.append(vector_axis_azimuth_deg(pa.e, pa.n, pb.e, pb.n))
    return azs


def segment_lengths(points: List[Point], ordered_idxs: List[int]) -> List[float]:
    ds: List[float] = []
    for a, b in zip(ordered_idxs[:-1], ordered_idxs[1:]):
        pa, pb = points[a], points[b]
        ds.append(dist(pa.e, pa.n, pb.e, pb.n))
    return ds


def mean_abs_turn_deg(azs: List[float]) -> float:
    if len(azs) < 2:
        return 0.0
    diffs = [axis_angle_diff_deg(a2, a1) for a1, a2 in zip(azs[:-1], azs[1:])]
    return sum(diffs) / len(diffs)


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def confidence_score(length: float,
                     npts: int,
                     strike_sigma: float,
                     dir_sigma: float,
                     turn_mean: float,
                     step_cv: float,
                     *,
                     max_link_distance: float,
                     min_points: int,
                     strike_match_tol_deg: float,
                     cone_half_angle_deg: float) -> float:
    """
    Conservative 0–100 score for cartography/QA and quick filtering (not a geologic truth metric).
    """
    if npts <= 1:
        return 0.0

    # Length support: how close average step is to the allowed max (short hop chains can be noisy)
    avg_step = length / max(1, (npts - 1))
    wL = clamp(avg_step / max(1e-9, max_link_distance), 0.0, 1.0)

    # Strike coherence: penalize scatter relative to tolerance
    wS = clamp(1.0 - (strike_sigma / max(1e-9, strike_match_tol_deg)), 0.0, 1.0)

    # Direction coherence (segment azimuth scatter): relative to cone
    wD = clamp(1.0 - (dir_sigma / max(1e-9, cone_half_angle_deg)), 0.0, 1.0)

    # Smoothness: mean absolute turn should be small
    wT = clamp(1.0 - (turn_mean / max(1e-9, 2.0 * cone_half_angle_deg)), 0.0, 1.0)

    # Continuity: step length variability
    wC = clamp(1.0 - (step_cv / 0.75), 0.0, 1.0)  # 0.75 is a practical cutoff

    # Support (points): ramp up after min_points
    wN = clamp((npts - min_points + 1) / max(1.0, float(min_points)), 0.0, 1.0)

    score = 100.0 * (wL ** 0.75) * (wS ** 1.0) * (wD ** 0.75) * (wT ** 0.75) * (wC ** 0.5) * (wN ** 0.5)
    return clamp(score, 0.0, 100.0)


def write_outputs(points: List[Point],
                  lineaments: List[List[int]],
                  segments: List[Tuple[int, int, int]],
                  out_points_csv: str = "lineaments_points.csv",
                  out_segments_csv: str = "lineaments_segments.csv",
                  out_summary_csv: str = "lineaments_summary.csv",
                  *,
                  max_link_distance: float,
                  min_points: int,
                  strike_match_tol_deg: float,
                  cone_half_angle_deg: float) -> None:
    """
    Writes:
      - Points CSV (ordered points per lineament)
      - Segments CSV (one row per segment, with skip flag)
      - Summary CSV with confidence metrics for plotting/QA
    """
    # stable ordering (by line)
    ordered_components: List[List[int]] = [order_by_line(points, comp) for comp in lineaments]

    # map point idx -> lineament id
    idx_to_lid: Dict[int, int] = {}
    for lid, comp in enumerate(ordered_components, start=1):
        for idx in comp:
            idx_to_lid[idx] = lid

    # Points export
    with open(out_points_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["LineamentID", "Seq", "PointIndex", "LineNumber", "Easting", "Northing", "Strike"])
        for lid, comp in enumerate(ordered_components, start=1):
            for seq, idx in enumerate(comp, start=1):
                p = points[idx]
                w.writerow([lid, seq, p.idx, p.line, p.e, p.n, p.strike])

    # Segments export
    with open(out_segments_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "LineamentID", "U_PointIndex", "V_PointIndex",
            "U_Easting", "U_Northing", "V_Easting", "V_Northing",
            "VectorAxisAzimuth", "SegmentLength", "SkippedLine"
        ])
        for u, v, sk in segments:
            lu = idx_to_lid.get(u)
            lv = idx_to_lid.get(v)
            if lu is None or lv is None or lu != lv:
                continue
            pu, pv = points[u], points[v]
            az = vector_axis_azimuth_deg(pu.e, pu.n, pv.e, pv.n)
            d = dist(pu.e, pu.n, pv.e, pv.n)
            w.writerow([lu, u, v, pu.e, pu.n, pv.e, pv.n, az, d, sk])

    # Summary export (center point + mean strike + length + metrics)
    with open(out_summary_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "LineamentID", "NumPoints",
            "CenterEasting", "CenterNorthing",
            "MeanStrike", "StrikeSigma",
            "Length",
            "MeanSegAzimuth", "SegAzSigma",
            "MeanAbsTurn",
            "StepMean", "StepStd", "StepCV",
            "Confidence"
        ])

        for lid, comp in enumerate(ordered_components, start=1):
            es = [points[i].e for i in comp]
            ns = [points[i].n for i in comp]
            ss = [points[i].strike for i in comp]

            center_e = sum(es) / len(es)
            center_n = sum(ns) / len(ns)
            mean_strike = mean_axis_deg(ss)
            strike_sigma = axis_circular_std_deg(ss)

            length = polyline_length(points, comp)
            seg_az = segment_axis_azimuths(points, comp)
            seg_sigma = axis_circular_std_deg(seg_az) if seg_az else 0.0
            mean_seg_az = mean_axis_deg(seg_az) if seg_az else float("nan")
            turn_mean = mean_abs_turn_deg(seg_az)

            steps = segment_lengths(points, comp)
            if steps:
                step_mean = sum(steps) / len(steps)
                step_var = sum((d - step_mean) ** 2 for d in steps) / max(1, len(steps))
                step_std = math.sqrt(step_var)
                step_cv = (step_std / step_mean) if step_mean > 0 else 0.0
            else:
                step_mean = step_std = step_cv = 0.0

            conf = confidence_score(
                length=length,
                npts=len(comp),
                strike_sigma=strike_sigma,
                dir_sigma=seg_sigma,
                turn_mean=turn_mean,
                step_cv=step_cv,
                max_link_distance=max_link_distance,
                min_points=min_points,
                strike_match_tol_deg=strike_match_tol_deg,
                cone_half_angle_deg=cone_half_angle_deg
            )

            w.writerow([
                lid, len(comp),
                center_e, center_n,
                mean_strike, strike_sigma,
                length,
                mean_seg_az, seg_sigma,
                turn_mean,
                step_mean, step_std, step_cv,
                conf
            ])


def write_debug_candidates(debug_rows: List[Dict[str, Any]], out_csv: str) -> None:
    if not debug_rows:
        return
    # unify keys
    keys: List[str] = sorted({k for r in debug_rows for k in r.keys()})
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in debug_rows:
            w.writerow(r)


# ----------------------------
# CLI
# ----------------------------

def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Build lineaments by linking strike-consistent points on neighbouring lines.")

    # Input format
    ap.add_argument("input_path", help="Input CSV or XYZ file")
    ap.add_argument("--format", choices=["csv", "xyz"], default="csv", help="Input file format")

    # CSV columns
    ap.add_argument("--col_line", default="LineNumber")
    ap.add_argument("--col_easting", default="Easting")
    ap.add_argument("--col_northing", default="Northing")
    ap.add_argument("--col_strike", default="Strike")

    # XYZ indices (0-based)
    ap.add_argument("--idx_line", type=int, default=0)
    ap.add_argument("--idx_easting", type=int, default=1)
    ap.add_argument("--idx_northing", type=int, default=2)
    ap.add_argument("--idx_strike", type=int, default=3)
    ap.add_argument("--xyz_has_header", action="store_true")

    # Tolerances
    ap.add_argument("--cone_half_angle_deg", type=float, default=10.0, help="Angular cone half-angle around strike (deg)")
    ap.add_argument("--strike_match_tol_deg", type=float, default=10.0, help="Strike agreement tolerance between points (deg)")
    ap.add_argument("--max_link_distance", type=float, default=200.0, help="Max distance to link points (map units)")
    ap.add_argument("--cell_size", type=float, default=200.0, help="Spatial hash cell size (map units)")

    ap.add_argument("--min_points", type=int, default=5, help="Minimum points required to save a lineament (>= this)")
    ap.add_argument("--require_reciprocal", action="store_true", help="Require mutual best-match on adjacent lines")

    # Skip-one-line bridging
    ap.add_argument("--skip_lines", type=int, default=0, help="Allow bridging across this many missing lines (0=off; 1 recommended)")
    ap.add_argument("--skip_distance_factor", type=float, default=1.6, help="Distance multiplier when bridging (default 1.6)")
    ap.add_argument("--skip_angle_bonus_deg", type=float, default=2.0, help="Add this many degrees to cone and strike tol when bridging")

    # Ray vs axis-line behaviour
    ap.add_argument("--allow_backward", action="store_true",
                    help="If set, do NOT enforce along>0 (treat strike as a bidirectional axis line). Default enforces a forward ray.")

    # Outputs
    ap.add_argument("--out_points", default="lineaments_points.csv")
    ap.add_argument("--out_segments", default="lineaments_segments.csv")
    ap.add_argument("--out_summary", default="lineaments_summary.csv")
    ap.add_argument("--out_debug_candidates", default="", help="If set, writes a CSV of candidate evaluations for tuning.")

    args = ap.parse_args()

    # Read points
    if args.format == "csv":
        pts = read_points_csv(
            args.input_path,
            col_line=args.col_line,
            col_e=args.col_easting,
            col_n=args.col_northing,
            col_strike=args.col_strike
        )
    else:
        pts = read_points_xyz(
            args.input_path,
            idx_line=args.idx_line,
            idx_e=args.idx_easting,
            idx_n=args.idx_northing,
            idx_strike=args.idx_strike,
            has_header=args.xyz_has_header
        )

    debug_rows: Optional[List[Dict[str, Any]]] = [] if args.out_debug_candidates else None

    comps, segs = build_lineaments(
        pts,
        cone_half_angle_deg=args.cone_half_angle_deg,
        strike_match_tol_deg=args.strike_match_tol_deg,
        max_link_distance=args.max_link_distance,
        min_points=args.min_points,
        cell_size=args.cell_size,
        require_reciprocal=args.require_reciprocal,
        skip_lines=args.skip_lines,
        skip_distance_factor=args.skip_distance_factor,
        skip_angle_bonus_deg=args.skip_angle_bonus_deg,
        debug_rows=debug_rows,
        enforce_forward_ray=not args.allow_backward
    )

    write_outputs(
        pts, comps, segs,
        out_points_csv=args.out_points,
        out_segments_csv=args.out_segments,
        out_summary_csv=args.out_summary,
        max_link_distance=args.max_link_distance,
        min_points=args.min_points,
        strike_match_tol_deg=args.strike_match_tol_deg,
        cone_half_angle_deg=args.cone_half_angle_deg
    )

    if args.out_debug_candidates and debug_rows is not None:
        write_debug_candidates(debug_rows, args.out_debug_candidates)

    print(f"Input points: {len(pts)}")
    print(f"Saved lineaments (>= {args.min_points} pts): {len(comps)}")
    print(f"Total kept segments: {len(segs)}")
    print(f"Wrote: {args.out_points}")
    print(f"Wrote: {args.out_segments}")
    print(f"Wrote: {args.out_summary}")
    if args.out_debug_candidates:
        print(f"Wrote: {args.out_debug_candidates}")


if __name__ == "__main__":
    main()
