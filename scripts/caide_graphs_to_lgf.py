#!/usr/bin/env python3
import sys
from collections import defaultdict

def rel_cap(rel: int) -> float:
    return 1.0

def main(inp: str, outp: str):
    asn_to_id = {}
    def get_id(asn: int) -> int:
        if asn not in asn_to_id:
            asn_to_id[asn] = len(asn_to_id)
        return asn_to_id[asn]

    # undirected edge aggregation (u < v)
    cap = defaultdict(float)

    with open(inp, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            a = int(parts[0])
            b = int(parts[1])
            r = int(parts[2])

            u = get_id(a)
            v = get_id(b)
            if u == v:
                continue

            x, y = (u, v) if u < v else (v, u)
            # Dedup with max-cap; or sum-cap if you prefer:
            cap[(x, y)] = max(cap[(x, y)], rel_cap(r))
            # cap[(x, y)] += rel_cap(r)  # alternative: aggregate

    n = len(asn_to_id)
    m = len(cap)

    with open(outp, "w", encoding="utf-8") as g:
        g.write("@nodes\n")
        g.write("label\n")
        for i in range(n):
            g.write(f"{i}\n")

        g.write("@arcs\n")
        ctr = 0
        g.write("  label     cost\n")
        # if you don’t use len, keep it 1.0 or omit depending on your parser
        for (u, v), c in cap.items():
            g.write(f"{u} {v} {ctr} 1.0\n")
            ctr+= 1

    print(f"Wrote {outp}: n={n}, m={m}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: caida_to_lgf.py as-caida20040105.txt out.lgf")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
