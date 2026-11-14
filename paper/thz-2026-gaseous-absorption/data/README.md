# Data files

Every input the analysis reads is in this directory, so the scripts in `sim/`
do not depend on any path outside the paper folder.

| File | Source | SHA-256 |
|---|---|---|
| p676_13_annex1_lines.json | ITU-R P.676-13 (08/2022) Annex 1, Table 1 (44 oxygen lines, f0, a1..a6) and Table 2 (35 water-vapour lines including the 1780 GHz pseudo-line, f0, b1..b6). Transcribed in the thz-ntn module (`contrib/thz-ntn/tools/p676-lut-gen.py`, the same table the C++ model compiles in) and checked equal, coefficient for coefficient (max difference 0.0), to the ITU-Rpy 0.4.0 files `v12_lines_oxygen.txt` and `v12_lines_water_vapour.txt`; revision 13 left Tables 1 and 2 unchanged from revision 12. | 5d4ad76fccc4283b75b78f131fa03df7eec5f7134cd2a0a4c4f7701cc709c015 |
| p676_13_annex2_part1.txt | ITU-R P.676-13 Annex 2 data Part 1 (h_o coefficients), https://www.itu.int/dms_pub/itu-r/oth/11/01/R11010000020001TXTM.txt | 9c38d278274d156c2c4dcacadff9d5a9d591e467acb3b04fd330c8b03694e7be |
| p676_13_annex2_part2.txt | ITU-R P.676-13 Annex 2 data Part 2 (K_V coefficients), https://www.itu.int/dms_pub/itu-r/oth/11/01/R11010000020002TXTM.txt | a57cb11235018c77a6ccea0bc985df6a1b2c6a841c6dfe9e65154e962ff2f95d |
| CG-3M3J-13-ValEx-Rev8.3.0.xlsx | ITU-R SG3 validation examples, https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev8.3.0.xlsx (downloaded 2026-10-19) | 5d2f4b6789de111de50b2be0cc91d9518c873d0e97299e4328abce55088af768 |
| sondes_nyalesund.npz | Ny-Alesund (WMO 01004) radiosondes, 12 UTC, 1 Dec 2023 to 29 Feb 2024 and July 2024, University of Wyoming upper-air archive, https://weather.uwyo.edu/wsgi/sounding (TEXT:CSV, BUFR source), downloaded 2026-10-19 with `../fetch_sondes.sh` and resampled by `../prep_sondes.py` (soundings that reach 300 hPa); kept only to mark which days the closed form had seen when it was written | ebc715390a077bdf10f80515d6b704182ce6570e3c736dfa3640b72119f8ae0a |
| sondes_nyalesund_year.npz | The same archive and station, 12 UTC, every day from 1 Dec 2023 to 30 Nov 2024, downloaded 2026-11-03 with `../fetch_sondes.sh` and resampled by `../prep_sondes.py --year` with the same filter; drives `analysis.py sonde`, `heff_model.py` and `sonde_rate.py` (356 soundings; the 109 of the set above are identical in it, and the archive now also returns 10 further December and January days) | bffb8bec27180b20abd054254b3404b13a07273d1d4ef4e2fa2315d9cf489a76 |
| sondes_arctic.npz | The same archive, 12 UTC, stations 01028, 04320, 70026, 71082, 71917, 02836, every second day from Dec 2023 to Feb 2024 and every sixth day from Mar to Nov 2024, downloaded 2026-11-03 with `../fetch_arctic.sh` and resampled by `../arctic.py --prep`, same filter | 1e943e6160ebd6ca7918157f095165998109230aa334b49c65e1cb64da93706b |
| raw_soundings.tar.xz | The raw TEXT:CSV files behind both sets above, as downloaded on 2026-11-03; the first-record question cannot be answered from the resampled files, and the archive can change | 7af287de2d8555574bf1c57820af214cc722b8fc7de80f8fe96b284bba710fcc |

Both Annex 2 coefficient files were re-downloaded on 2026-10-19 and are
byte-identical to the copies used.

Pinned environment: `../requirements.txt` (Python packages) and am 14.0
(https://lweb.cfa.harvard.edu/~spaine/am/, Zenodo DOI 10.5281/zenodo.13748403).
The code, these data and the cached outputs in `../out/` are to be archived
with a DOI on acceptance (open item for the author).

Published measured values (round 3, review D13) are transcribed, with the
table or equation each comes from, in `../am_check/meas_data.json` and read
by `../am_check/meas_check2.py`:
Archibald et al. 2002 (MNRAS 336, 1, Table 1; DOI 10.1046/j.1365-8711.2002.05582.x),
Dempsey et al. 2013 (MNRAS 430, 2534, Sec. 3.2; DOI 10.1093/mnras/stt090) with the
SCUBA-2 filter half-power bands of Holland et al. 2013 (MNRAS 430, 2513, Sec. 2;
DOI 10.1093/mnras/sts612), Radford and Peterson 2016 (PASP 128, 075001,
eqs. 8 and 9; DOI 10.1088/1538-3873/128/965/075001), Chamberlin, Lane and
Stark 1997 (ApJ 476, 428, Table 2; DOI 10.1086/303621), Yang, Mandehgar and
Grischkowsky 2014 (Opt. Express 22, 4388, pp. 4390 and 4399;
DOI 10.1364/OE.22.004388) and Yang, Shutler and Grischkowsky 2011
(Opt. Express 19, 8830, Table 1; DOI 10.1364/OE.19.008830). Every value was
read from the full text and every DOI checked against Crossref on 2026-10-20.
