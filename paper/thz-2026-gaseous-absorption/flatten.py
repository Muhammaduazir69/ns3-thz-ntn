#!/usr/bin/env python3
"""Expand \\input{...} recursively (macros.tex excluded) for screening tools."""
import re, sys, os
def expand(path):
    out = []
    for line in open(path):
        m = re.match(r"\s*\\input\{([^}]+)\}", line)
        if m and "macros" not in m.group(1):
            p = m.group(1) if m.group(1).endswith(".tex") else m.group(1) + ".tex"
            out.append(expand(p))
        else:
            out.append(line)
    return "".join(out)
open(sys.argv[2], "w").write(expand(sys.argv[1]))
