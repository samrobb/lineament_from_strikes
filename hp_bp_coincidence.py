#!/usr/bin/env python3
"""
HP/BP coincidence + intersections helper.

Inputs:
- HP lineament outputs (summary + segments)
- BP lineament outputs (summary + segments)

Outputs:
- hp_bp_matches.csv : best BP match for each HP lineament (by distance + strike similarity)
- hp_bp_classified.csv : HP lineaments labelled as "HP-only", "HP+BP coincident"
- hp_bp_intersections.csv : intersection points between HP and BP segment sets (with ids and angle info)

Typical use:
python hp_bp_coincidence.py \
  --hp_summary hp_lineaments_summary.csv --hp_segments hp_lineaments_segments.csv \
  --bp_summary bp_lineaments_summary.csv --bp_segments bp_lineaments_segments.csv \
  --match_dist 1500 --strike_tol 15 --cell 2000
"""
import csv
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable

# ---- helpers

def normalize_axis_deg(a: float) -> float:
    a = a % 180.0
    if a < 0:
        a += 180.0
    return a

def axis_angle_diff_deg(a: float, b: float) -> float:
    d = abs(normalize_axis_deg(a) - normalize_axis_deg(b)) % 180.0
    return min(d, 180.0 - d)

def dist(e1: float, n1: float, e2: float, n2: float) -> float:
    return math.hypot(e2 - e1, n2 - n1)

def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))

# ---- data

@dataclass
class Lineament:
    lid: int
    e: float
    n: float
    strike: float
    length: float
    conf: Optional[float] = None

@dataclass
class Segment:
    lid: int
    e1: float
    n1: float
    e2: float
    n2: float
    az: float

# ---- readers

def read_summary(path: str) -> Dict[int, Lineament]:
    out: Dict[int, Lineament] = {}
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        for row in r:
            lid = int(row["LineamentID"])
            out[lid] = Lineament(
                lid=lid,
                e=float(row["CenterEasting"]),
                n=float(row["CenterNorthing"]),
                strike=float(row["MeanStrike"]),
                length=float(row.get("Length", "0") or 0),
                conf=float(row.get("Confidence")) if row.get("Confidence") not in (None, "", "nan", "NaN") else None
            )
    return out

def read_segments(path: str) -> List[Segment]:
    segs: List[Segment] = []
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        for row in r:
            segs.append(Segment(
                lid=int(row["LineamentID"]),
                e1=float(row["U_Easting"]),
                n1=float(row["U_Northing"]),
                e2=float(row["V_Easting"]),
                n2=float(row["V_Northing"]),
                az=float(row.get("VectorAxisAzimuth", "nan"))
            ))
    return segs

# ---- matching

def best_match(hp: Lineament, bp_list: List[Lineament], *, match_dist: float, strike_tol: float) -> Optional[Tuple[Lineament, float, float]]:
    best = None  # (score, bp, d, sdiff)
    for bp in bp_list:
        d = dist(hp.e, hp.n, bp.e, bp.n)
        if d > match_dist:
            continue
        sdiff = axis_angle_diff_deg(hp.strike, bp.strike)
        if sdiff > strike_tol:
            continue
        # score: weighted distance + strike diff
        score = (d / match_dist) + 0.7 * (sdiff / strike_tol)
        if best is None or score < best[0]:
            best = (score, bp, d, sdiff)
    if best is None:
        return None
    _, bp, d, sdiff = best
    return bp, d, sdiff

# ---- intersections (segment vs segment)

def seg_bbox(s: Segment) -> Tuple[float, float, float, float]:
    return min(s.e1, s.e2), min(s.n1, s.n2), max(s.e1, s.e2), max(s.n1, s.n2)

def ccw(ax, ay, bx, by, cx, cy) -> float:
    return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax)

def segment_intersection(a: Segment, b: Segment) -> Optional[Tuple[float, float]]:
    # Fast reject bbox
    ax0, ay0, ax1, ay1 = seg_bbox(a)
    bx0, by0, bx1, by1 = seg_bbox(b)
    if ax1 < bx0 or bx1 < ax0 or ay1 < by0 or by1 < ay0:
        return None

    x1,y1,x2,y2 = a.e1,a.n1,a.e2,a.n2
    x3,y3,x4,y4 = b.e1,b.n1,b.e2,b.n2

    d1 = ccw(x1,y1,x2,y2,x3,y3)
    d2 = ccw(x1,y1,x2,y2,x4,y4)
    d3 = ccw(x3,y3,x4,y4,x1,y1)
    d4 = ccw(x3,y3,x4,y4,x2,y2)

    # Proper intersection test (including colinear simplification later)
    if (d1 == 0 and d2 == 0 and d3 == 0 and d4 == 0):
        # Colinear: skip (ambiguous "intersection point")
        return None
    if (d1 * d2 > 0) or (d3 * d4 > 0):
        return None

    # Compute intersection of lines (not segments) then ensure within segment bounds
    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    if abs(denom) < 1e-12:
        return None
    px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)) / denom

    # Bounds check with small epsilon
    eps = 1e-9
    if (min(x1,x2)-eps <= px <= max(x1,x2)+eps and min(y1,y2)-eps <= py <= max(y1,y2)+eps and
        min(x3,x4)-eps <= px <= max(x3,x4)+eps and min(y3,y4)-eps <= py <= max(y3,y4)+eps):
        return px, py
    return None

# Spatial hash for segments
def seg_cell_coords(x: float, y: float, cell: float) -> Tuple[int,int]:
    return int(math.floor(x / cell)), int(math.floor(y / cell))

def build_seg_grid(segs: List[Segment], cell: float) -> Dict[Tuple[int,int], List[int]]:
    grid: Dict[Tuple[int,int], List[int]] = {}
    for i,s in enumerate(segs):
        x0,y0,x1,y1 = seg_bbox(s)
        cx0,cy0 = seg_cell_coords(x0,y0,cell)
        cx1,cy1 = seg_cell_coords(x1,y1,cell)
        for cx in range(cx0, cx1+1):
            for cy in range(cy0, cy1+1):
                grid.setdefault((cx,cy), []).append(i)
    return grid

def iter_nearby_seg_indices(grid: Dict[Tuple[int,int], List[int]], s: Segment, cell: float) -> Iterable[int]:
    x0,y0,x1,y1 = seg_bbox(s)
    cx0,cy0 = seg_cell_coords(x0,y0,cell)
    cx1,cy1 = seg_cell_coords(x1,y1,cell)
    seen = set()
    for cx in range(cx0, cx1+1):
        for cy in range(cy0, cy1+1):
            for idx in grid.get((cx,cy), []):
                if idx not in seen:
                    seen.add(idx)
                    yield idx

def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="HP/BP coincidence matching and intersection extraction.")

    ap.add_argument("--hp_summary", required=True)
    ap.add_argument("--hp_segments", required=True)
    ap.add_argument("--bp_summary", required=True)
    ap.add_argument("--bp_segments", required=True)

    ap.add_argument("--match_dist", type=float, default=1500.0, help="Center-to-center distance tolerance (map units)")
    ap.add_argument("--strike_tol", type=float, default=15.0, help="Strike tolerance for coincidence (deg)")
    ap.add_argument("--cell", type=float, default=2000.0, help="Spatial hash cell size for intersections")
    ap.add_argument("--out_matches", default="hp_bp_matches.csv")
    ap.add_argument("--out_classified", default="hp_bp_classified.csv")
    ap.add_argument("--out_intersections", default="hp_bp_intersections.csv")

    args = ap.parse_args()

    hp = read_summary(args.hp_summary)
    bp = read_summary(args.bp_summary)
    hp_list = list(hp.values())
    bp_list = list(bp.values())

    # Match
    matches: List[Dict[str, object]] = []
    classified: List[Dict[str, object]] = []
    coincident_pairs = set()

    for h in hp_list:
        m = best_match(h, bp_list, match_dist=args.match_dist, strike_tol=args.strike_tol)
        if m is None:
            matches.append({"HP_ID": h.lid, "BP_ID": "", "CenterDist": "", "StrikeDiff": ""})
            classified.append({"HP_ID": h.lid, "Class": "HP-only", "HP_Strike": h.strike, "HP_Length": h.length, "HP_Confidence": h.conf})
        else:
            b, d, sdiff = m
            coincident_pairs.add((h.lid, b.lid))
            matches.append({"HP_ID": h.lid, "BP_ID": b.lid, "CenterDist": d, "StrikeDiff": sdiff})
            classified.append({"HP_ID": h.lid, "Class": "HP+BP coincident", "HP_Strike": h.strike, "HP_Length": h.length, "HP_Confidence": h.conf,
                               "BP_ID": b.lid, "BP_Strike": b.strike, "BP_Length": b.length, "BP_Confidence": b.conf, "CenterDist": d, "StrikeDiff": sdiff})

    with open(args.out_matches, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["HP_ID","BP_ID","CenterDist","StrikeDiff"])
        w.writeheader()
        w.writerows(matches)

    # classified headers: union
    cls_keys = sorted({k for r in classified for k in r.keys()})
    with open(args.out_classified, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cls_keys)
        w.writeheader()
        w.writerows(classified)

    # Intersections
    hp_segs = read_segments(args.hp_segments)
    bp_segs = read_segments(args.bp_segments)
    bp_grid = build_seg_grid(bp_segs, args.cell)

    intersections: List[Dict[str, object]] = []
    for s in hp_segs:
        for j in iter_nearby_seg_indices(bp_grid, s, args.cell):
            t = bp_segs[j]
            ip = segment_intersection(s, t)
            if ip is None:
                continue
            x, y = ip
            adiff = axis_angle_diff_deg(s.az, t.az) if (not math.isnan(s.az) and not math.isnan(t.az)) else ""
            intersections.append({
                "Easting": x, "Northing": y,
                "HP_ID": s.lid, "BP_ID": t.lid,
                "HP_SegAz": s.az, "BP_SegAz": t.az,
                "SegAzDiff": adiff
            })

    with open(args.out_intersections, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Easting","Northing","HP_ID","BP_ID","HP_SegAz","BP_SegAz","SegAzDiff"])
        w.writeheader()
        w.writerows(intersections)

    print(f"HP lineaments: {len(hp_list)}")
    print(f"BP lineaments: {len(bp_list)}")
    print(f"Coincident HP: {sum(1 for r in classified if r.get('Class')=='HP+BP coincident')}")
    print(f"Intersections: {len(intersections)}")
    print(f"Wrote: {args.out_matches}")
    print(f"Wrote: {args.out_classified}")
    print(f"Wrote: {args.out_intersections}")

if __name__ == "__main__":
    main()
