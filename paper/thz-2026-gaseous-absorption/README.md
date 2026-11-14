# Code and data: simplified gaseous absorption models against ITU-R P.676-13, 340 to 850 GHz

Companion material for M. Uzair, "Deviation of Simplified Gaseous Absorption Models From Recommendation P.676-13 on Low Earth Orbit Slant Paths From 340 to 850 GHz," 2026.

The study verifies the module's molecular absorption engine against Recommendation ITU-R P.676-13 and the am radiative-transfer code, then measures how far simplified absorption models deviate on LEO slant paths, using radiosonde profiles and the standard atmospheres.

| Path | Contents |
|---|---|
| `run_all.sh` | Regenerates every number, table and figure; options in its header (`SKIP_CPP`, `SKIP_AM`, `PAR`) |
| `thzatm.py`, `run_models*.py`, `analysis.py`, `heff_model.py`, `budget.py` | Absorption models, slant integration, statistics and link budgets |
| `cpp_sweeps.sh` | Sweeps of the module's C++ engine (`thz-ntn-absorption-sweep` example) |
| `fetch_sondes.sh`, `prep_sondes.py`, `fetch_arctic.sh`, `arctic.py`, `data/` | Radiosonde inputs and line tables, with sources |
| `am_check/` | am 14.0 configurations and outputs used as the independent reference |
| `out.tar.xz.part-*` | All result files |
| `check_macros.py`, `test_annex2_rev13.py`, `p835v7_check.py` | Consistency and reference checks |
| `requirements.txt` | Pinned Python packages |

Large result folders are stored as xz archives. Where an archive is split into `*.part-NN` files, join and check it first: `cat out.tar.xz.part-* > out.tar.xz && sha256sum -c out.tar.xz.sha256 && tar -xJf out.tar.xz`.

Contact: Muhammad Uzair, muhammaduzairr69@gmail.com, ORCID 0009-0002-4104-2680. License: as the repository (GPL-2.0-only for code).
