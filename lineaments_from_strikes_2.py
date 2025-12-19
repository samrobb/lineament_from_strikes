import math
import csv
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import List, Dict, Tuple, Set, Optional


@dataclass(frozen=True)
class Point:
    idx: int
    line: int
    e: float
    n: float
    strike: float  # degrees, expected 0-180 (but we normalize)


def axis_azimuth_deg_from_to(p: Point, q: Point) -> float:
    """Axis azimuth [0,180) of vector p->q."""
    return vector_axis_azimuth_deg(p.e, p.n, q.e, q.n)


def ray_perp_distance(p: Point, strike_axis_deg: float, q: Point) -> Tuple[float, float]:
    """
    Distance metrics relative to a ray from p in strike direction (axis):
      - along: signed distance along the ray direction (must be >0 to be 'forward')
      - perp: perpendicular distance to the ray line

    Uses a unit direction vector for strike_axis_deg interpreted as azimuth-from-North.
    """
    a = normalize_axis_deg(strike_axis_deg)
    # azimuth is from North clockwise; direction unit vector in (E,N):
    rad = math.radians(a)
    ux = math.sin(rad)  # East component
    uy = math.cos(rad)  # North component

    dx = q.e - p.e
    dy = q.n - p.n

    along = dx * ux + dy * uy
    perp = abs(dx * (-uy) + dy * (ux))  # magnitude of cross-track
    return along, perp

def order_by_line(points: List[Point], idxs: List[int]) -> List[int]:
    """Order a set of point indices by (LineNumber, then along-line position)."""
    return sorted(idxs, key=lambda i: (points[i].line, points[i].e, points[i].n))


def polyline_length(points: List[Point], ordered_idxs: List[int]) -> float:
    """Sum of segment lengths along an ordered polyline."""
    if len(ordered_idxs) < 2:
        return 0.0
    total = 0.0
    for a, b in zip(ordered_idxs[:-1], ordered_idxs[1:]):
        pa, pb = points[a], points[b]
        total += dist(pa.e, pa.n, pb.e, pb.n)
    return total


def best_match_on_line(points: List[Point],
                       grid: Dict[Tuple[int, int], List[int]],
                       p_idx: int,
                       target_line: int,
                       max_link_distance: float,
                       cell_size: float,
                       cone_half_angle_deg: float,
                       strike_match_tol_deg: float) -> Optional[int]:
    """
    From point p, find the single best matching point q on target_line within:
      - distance <= max_link_distance
      - azimuth(P->Q) within cone_half_angle of P.strike
      - Q.strike within strike_match_tol of P.strike

    Tie-break: smallest angular misfit, then smallest perpendicular distance to ray,
              then smallest Euclidean distance.
    """
    p = points[p_idx]
    cand_idxs = nearby_candidates(grid, points, p.e, p.n, max_link_distance, cell_size)

    best = None  # (score_tuple, q_idx)
    for j in cand_idxs:
        q = points[j]
        if q.line != target_line:
            continue

        d = dist(p.e, p.n, q.e, q.n)
        if d == 0 or d > max_link_distance:
            continue

        az_axis = vector_axis_azimuth_deg(p.e, p.n, q.e, q.n)

        ang_misfit = axis_angle_diff_deg(az_axis, p.strike)
        if ang_misfit > cone_half_angle_deg:
            continue

        strike_misfit = axis_angle_diff_deg(q.strike, p.strike)
        if strike_misfit > strike_match_tol_deg:
            continue

        along, perp = ray_perp_distance(p, p.strike, q)

        score = (ang_misfit, perp, d, strike_misfit)
        if best is None or score < best[0]:
            best = (score, j)

    return None if best is None else best[1]


def normalize_axis_deg(a: float) -> float:
    """
    Normalize an axis direction to [0, 180).
    Strike is an axis: 0 and 180 are the same.
    """
    a = a % 180.0
    if a < 0:
        a += 180.0
    return a


def axis_angle_diff_deg(a: float, b: float) -> float:
    """
    Smallest difference between two axis directions (both in [0,180)).
    """
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


def vector_axis_azimuth_deg(e1: float, n1: float, e2: float, n2: float) -> float:
    """
    Azimuth of the vector from p1->p2 as an axis direction [0,180).
    (i.e., 0-180, where 0 ~ North, 90 ~ East)
    """
    de = e2 - e1
    dn = n2 - n1
    # azimuth from North, clockwise:
    az = math.degrees(math.atan2(de, dn))  # atan2(x=east, y=north)
    az = az % 360.0
    # convert to axis (0-180)
    if az >= 180.0:
        az -= 180.0
    return az


def dist(e1: float, n1: float, e2: float, n2: float) -> float:
    de = e2 - e1
    dn = n2 - n1
    return math.hypot(de, dn)


def parse_line_number(s: str) -> int:
    """
    Extract the integer from line labels like 'L12', 'L0012', 'l12', 'Line L12'.
    Falls back to plain int if the whole string is numeric.
    """
    s = str(s).strip()
    m = re.search(r'(\d+)', s)
    if m:
        return int(m.group(1))
    # if nothing numeric is found, this is a hard error
    raise ValueError(f"Could not parse LineNumber from: {s!r}")


def read_points_flexible(path: str, args) -> List[Point]:
    points: List[Point] = []

    # delimiter handling
    if args.delimiter.lower() == "space":
        delimiter = None  # split on whitespace
    else:
        delimiter = args.delimiter

    with open(path, "r", encoding="utf-8-sig") as f:
        if args.no_header:
            # Headerless (XYZ-style)
            if None in (args.idx_easting, args.idx_northing,
                        args.idx_strike, args.idx_line):
                raise ValueError("Headerless file requires --idx_* arguments")

            for i, line in enumerate(f):
                if not line.strip():
                    continue
                parts = line.split() if delimiter is None else line.split(delimiter)

                e = float(parts[args.idx_easting])
                n = float(parts[args.idx_northing])
                strike = normalize_axis_deg(float(parts[args.idx_strike]))
                line_num = parse_line_number(parts[args.idx_line])

                points.append(Point(i, line_num, e, n, strike))

        else:
            # Header-based (CSV or delimited text)
            import csv
            reader = csv.DictReader(f, delimiter=delimiter)

            for i, row in enumerate(reader):
                e = float(row[args.col_easting])
                n = float(row[args.col_northing])
                strike = normalize_axis_deg(float(row[args.col_strike]))
                line_num = parse_line_number(row[args.col_line])

                points.append(Point(i, line_num, e, n, strike))

    return points



def build_spatial_hash(points: List[Point], cell_size: float) -> Dict[Tuple[int, int], List[int]]:
    """
    Simple grid index: maps (cx,cy) -> list of point indices.
    """
    grid: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for p in points:
        cx = int(math.floor(p.e / cell_size))
        cy = int(math.floor(p.n / cell_size))
        grid[(cx, cy)].append(p.idx)
    return grid


def nearby_candidates(grid: Dict[Tuple[int, int], List[int]],
                      points: List[Point],
                      e: float, n: float,
                      radius: float,
                      cell_size: float) -> List[int]:
    """
    Return point indices that fall in grid cells within radius (coarse),
    caller should still distance-check.
    """
    cx = int(math.floor(e / cell_size))
    cy = int(math.floor(n / cell_size))
    r_cells = int(math.ceil(radius / cell_size))
    out: List[int] = []
    for dx in range(-r_cells, r_cells + 1):
        for dy in range(-r_cells, r_cells + 1):
            out.extend(grid.get((cx + dx, cy + dy), []))
    return out


def build_lineaments(points: List[Point],
                     cone_half_angle_deg: float = 10.0,
                     strike_match_tol_deg: float = 10.0,
                     max_link_distance: float = 200.0,
                     min_points: int = 5,
                     cell_size: float = 200.0,
                     require_reciprocal: bool = True) -> Tuple[List[List[int]], List[Tuple[int, int]]]:
    """
    Build lineaments as stepwise chains across adjacent lines (no connected-components merging).
    This prevents multiple points from the same line entering one lineament and reduces zig-zag.

    Returns:
      - list of lineaments (each is an ordered list of point indices)
      - list of segments (u,v) for all kept lineaments
    """
    # group point indices by numeric LineNumber
    by_line: Dict[int, List[int]] = defaultdict(list)
    for p in points:
        by_line[p.line].append(p.idx)

    sorted_lines = sorted(by_line.keys())
    if len(sorted_lines) < 2:
        return [], []

    # next/prev line lookup by acquisition order
    next_line: Dict[int, Optional[int]] = {}
    prev_line: Dict[int, Optional[int]] = {}
    for i, ln in enumerate(sorted_lines):
        prev_line[ln] = sorted_lines[i - 1] if i > 0 else None
        next_line[ln] = sorted_lines[i + 1] if i < len(sorted_lines) - 1 else None

    # spatial grids per line
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
    kept_edges: List[Tuple[int, int]] = []

    # helper: attempt to link p -> q on next line with optional reciprocal requirement
    def link_forward(p_idx: int) -> Optional[int]:
        p = points[p_idx]
        nl = next_line[p.line]
        if nl is None:
            return None

        q_idx = best_match_on_line(
            points, grid_by_line[nl], p_idx, nl,
            max_link_distance, cell_size,
            cone_half_angle_deg, strike_match_tol_deg
        )
        if q_idx is None or q_idx in used:
            return None

        if require_reciprocal:
            # q must also pick p as its best match on q's previous line
            pl = prev_line[points[q_idx].line]
            if pl is None:
                return None
            back = best_match_on_line(
                points, grid_by_line[pl], q_idx, pl,
                max_link_distance, cell_size,
                cone_half_angle_deg, strike_match_tol_deg
            )
            if back != p_idx:
                return None

        return q_idx

    # Build chains starting from every unused point, but only accept chains with >= min_points
    # We also extend backward (using prev_line) by building reverse links, then stitch.
    def link_backward(p_idx: int) -> Optional[int]:
        p = points[p_idx]
        pl = prev_line[p.line]
        if pl is None:
            return None

        q_idx = best_match_on_line(
            points, grid_by_line[pl], p_idx, pl,
            max_link_distance, cell_size,
            cone_half_angle_deg, strike_match_tol_deg
        )
        if q_idx is None or q_idx in used:
            return None

        if require_reciprocal:
            # q must also pick p as its best match on q's next line
            nl = next_line[points[q_idx].line]
            if nl is None:
                return None
            fwd = best_match_on_line(
                points, grid_by_line[nl], q_idx, nl,
                max_link_distance, cell_size,
                cone_half_angle_deg, strike_match_tol_deg
            )
            if fwd != p_idx:
                return None

        return q_idx

    # iterate lines in order to seed sensibly
    all_idxs = sorted(range(len(points)), key=lambda i: (points[i].line, points[i].e, points[i].n))

    for seed in all_idxs:
        if seed in used:
            continue

        # grow backward
        backward_chain = [seed]
        cur = seed
        while True:
            prev_idx = link_backward(cur)
            if prev_idx is None or prev_idx in backward_chain:
                break
            backward_chain.append(prev_idx)
            cur = prev_idx

        backward_chain = list(reversed(backward_chain))  # now earliest -> seed

        # grow forward from seed (but avoid duplicating seed)
        forward_chain = []
        cur = seed
        while True:
            nxt = link_forward(cur)
            if nxt is None or nxt in forward_chain:
                break
            forward_chain.append(nxt)
            cur = nxt

        chain = backward_chain + forward_chain

        # enforce one point per line (should already hold, but keep it bulletproof)
        seen_lines = set()
        dedup_chain = []
        for idx in chain:
            ln = points[idx].line
            if ln in seen_lines:
                continue
            seen_lines.add(ln)
            dedup_chain.append(idx)

        if len(dedup_chain) >= min_points:
            # mark used and save
            for idx in dedup_chain:
                used.add(idx)
            lineaments.append(dedup_chain)

            # edges for outputs
            for a, b in zip(dedup_chain[:-1], dedup_chain[1:]):
                u, v = (a, b) if a < b else (b, a)
                kept_edges.append((u, v))

    return lineaments, kept_edges



def order_component(points: List[Point], comp: List[int]) -> List[int]:
    """
    Give the lineament a stable point order for export.
    Approach: principal-axis projection (PCA-ish) using covariance eigenvector.
    """
    xs = [points[i].e for i in comp]
    ys = [points[i].n for i in comp]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)

    # covariance
    sxx = sum((x - mx) ** 2 for x in xs) / len(xs)
    syy = sum((y - my) ** 2 for y in ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / len(xs)

    # eigenvector of largest eigenvalue for 2x2 covariance
    # handle near-degenerate cases
    if abs(sxy) < 1e-12 and abs(sxx - syy) < 1e-12:
        # basically a blob: sort by Easting then Northing
        return sorted(comp, key=lambda i: (points[i].e, points[i].n))

    # compute principal direction angle
    theta = 0.5 * math.atan2(2 * sxy, (sxx - syy))
    ux = math.cos(theta)
    uy = math.sin(theta)

    # project onto principal axis and sort
    def proj(i: int) -> float:
        return (points[i].e - mx) * ux + (points[i].n - my) * uy

    return sorted(comp, key=proj)

def mean_strike_axis_deg(strikes_deg: List[float]) -> float:
    """
    Mean of strike values treated as an AXIS (0-180).
    Uses the 2-theta method:
      mean_axis = 0.5 * atan2(sum(sin 2θ), sum(cos 2θ))
    Returns in [0, 180).
    """
    if not strikes_deg:
        return float("nan")

    s2 = 0.0
    c2 = 0.0
    for a in strikes_deg:
        a = normalize_axis_deg(a)
        r = math.radians(2.0 * a)
        s2 += math.sin(r)
        c2 += math.cos(r)

    if abs(s2) < 1e-12 and abs(c2) < 1e-12:
        # Degenerate case (e.g., perfectly symmetric distribution)
        return float("nan")

    mean2 = math.degrees(math.atan2(s2, c2))  # in degrees, (-180,180]
    mean = 0.5 * mean2
    return normalize_axis_deg(mean)


def write_outputs(points: List[Point],
                  components: List[List[int]],
                  edges: List[Tuple[int, int]],
                  out_points_csv: str = "lineaments_points.csv",
                  out_segments_csv: str = "lineaments_segments.csv",
                  out_summary_csv: str = "lineaments_summary.csv") -> None:

    # map point idx -> lineament id
    idx_to_lid: Dict[int, int] = {}
    ordered_components: List[List[int]] = []
    for lid, comp in enumerate(components, start=1):
        ordered = order_by_line(points, comp)
        ordered_components.append(ordered)
        for i in ordered:
            idx_to_lid[i] = lid

    # Points export
    with open(out_points_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["LineamentID", "Seq", "PointIndex", "LineNumber", "Easting", "Northing", "Strike"])
        for lid, comp in enumerate(ordered_components, start=1):
            for seq, idx in enumerate(comp, start=1):
                p = points[idx]
                w.writerow([lid, seq, p.idx, p.line, p.e, p.n, p.strike])

    # Segments export (only edges where both endpoints belong to same saved lineament)
    with open(out_segments_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["LineamentID", "U_PointIndex", "V_PointIndex",
                    "U_Easting", "U_Northing", "V_Easting", "V_Northing", "VectorAxisAzimuth"])
        for u, v in edges:
            lu = idx_to_lid.get(u)
            lv = idx_to_lid.get(v)
            if lu is None or lv is None or lu != lv:
                continue
            pu, pv = points[u], points[v]
            az = vector_axis_azimuth_deg(pu.e, pu.n, pv.e, pv.n)
            w.writerow([lu, u, v, pu.e, pu.n, pv.e, pv.n, az])

    # Summary export (one point + one strike per lineament)
        # Summary export (one point + one strike + length per lineament)
    with open(out_summary_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["LineamentID", "NumPoints", "CenterEasting", "CenterNorthing", "MeanStrike", "Length"])
        for lid, comp in enumerate(ordered_components, start=1):
            # comp is already ordered in this new approach, but we’ll be safe:
            ordered = order_by_line(points, comp)

            es = [points[i].e for i in ordered]
            ns = [points[i].n for i in ordered]
            ss = [points[i].strike for i in ordered]

            center_e = sum(es) / len(es)
            center_n = sum(ns) / len(ns)
            mean_strike = mean_strike_axis_deg(ss)
            length = polyline_length(points, ordered)

            w.writerow([lid, len(ordered), center_e, center_n, mean_strike, length])


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Build lineaments by linking strike-consistent points on neighbouring lines.")
    ap.add_argument("input_path", help="Input CSV or XYZ file")

    # File formats
    ap.add_argument("--delimiter", default=",",
                    help="Field delimiter (',' for CSV, ' ' for XYZ). Use 'space' for whitespace.")
    ap.add_argument("--no_header", action="store_true",
                    help="Set if the file has no header row")

    # Column headers
    ap.add_argument("--col_easting", default="Easting")
    ap.add_argument("--col_northing", default="Northing")
    ap.add_argument("--col_strike", default="Strike")
    ap.add_argument("--col_line", default="LineNumber")

    # Column indicies (for XYZ)
    ap.add_argument("--idx_easting", type=int)
    ap.add_argument("--idx_northing", type=int)
    ap.add_argument("--idx_strike", type=int)
    ap.add_argument("--idx_line", type=int)

    # Linenament parameters
    ap.add_argument("--cone_half_angle", type=float, default=10.0)
    ap.add_argument("--strike_match_tol", type=float, default=10.0)
    ap.add_argument("--max_link_distance", type=float, default=200.0, help="Max distance to link points (map units)")
    ap.add_argument("--min_points", type=int, default=5, help="Minimum points for saving a lineament (>= this)")
    ap.add_argument("--cell_size", type=float, default=200.0, help="Spatial hash cell size (map units)")
    ap.add_argument("--out_points", default="lineaments_points.csv", help="Output points CSV")
    ap.add_argument("--out_segments", default="lineaments_segments.csv", help="Output segments CSV")
    args = ap.parse_args()

    pts = read_points_flexible(args.input_path, args)
    comps, eds = build_lineaments(
        pts,
        cone_half_angle_deg=args.cone_half_angle,
        strike_match_tol_deg=args.strike_match_tol,
        max_link_distance=args.max_link_distance,
        min_points=args.min_points,
        cell_size=args.cell_size,
        require_reciprocal=args.require_reciprocal
    )

    print("Lineaments found:", len(comps))
    print("Edges kept:", len(eds))
    write_outputs(pts, comps, eds, args.out_points, args.out_segments, "lineaments_summary.csv")

    print(f"Input points: {len(pts)}")
    print(f"Saved lineaments (>= {args.min_points} pts): {len(comps)}")
    print(f"Total edges: {len(eds)}")
    print(f"Wrote: {args.out_points}")
    print(f"Wrote: {args.out_segments}")


if __name__ == "__main__":
    main()
