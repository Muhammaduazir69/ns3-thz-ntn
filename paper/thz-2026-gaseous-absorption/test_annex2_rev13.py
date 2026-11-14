"""thzatm.p676_13_annex2 against the ITU-R SG3 validation workbook (P.676-13 Annex 2)."""
import math
import openpyxl
import thzatm as A

wb = openpyxl.load_workbook("data/CG-3M3J-13-ValEx-Rev8.3.0.xlsx", read_only=True, data_only=True)
worst = 0.0
n = 0
for r in list(wb["P.676-13 A_Gas_A2_INST"].iter_rows(values_only=True))[23:]:
    if r is None or len(r) < 33 or not isinstance(r[5], (int, float)):
        continue
    Ps, Ts, f, rho, el, ref = r[2], r[4], r[5], r[8], r[9], r[32]
    got = float(A.p676_13_annex2(f, el, rho, Ps, Ts)[0])
    worst = max(worst, abs(got - ref) / ref)
    n += 1
rows = list(wb["P.676-13 A_Gas_A2_STAT"].iter_rows(values_only=True))
hdr = list(rows[21])
w0 = [i for i, h in enumerate(hdr) if isinstance(h, str) and h.strip().startswith("f")][1]
m = 0
for r in rows[23:]:
    if not r or not isinstance(r[w0], (int, float)):
        continue
    f, Ps, Ts, rho, p, V, el = r[w0:w0 + 7]
    Awv = r[w0 + 13]
    # water-vapour part of the statistical method (eq. 40): K_V(mean surface) * V(p)
    Aw = float(A.p676_13_annex2(f, el, rho, Ps, Ts, V_t=V)[0]) - float(A.p676_13_annex2(f, el, rho, Ps, Ts, V_t=0.0)[0])
    worst = max(worst, abs(Aw - Awv) / Awv)
    m += 1
print(f"{n} instantaneous + {m} statistical cases, worst relative error {worst:.2e}")
assert worst < 1e-6
