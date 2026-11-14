#!/bin/sh
# Ny-Alesund (WMO 01004) 12 UTC soundings, 1 Dec 2023 to 30 Nov 2024, University
# of Wyoming archive. The December to February and July subset is the held-out
# test of the absorbing-height analysis; the whole year drives the sounding
# rule of the availability study.
# Usage: fetch_sondes.sh <dir>. Files already holding a sounding are kept, and
# a "server busy" reply is retried.
cd "$1" || exit 1
for d in $(python3 -c "
import datetime as D
s=D.date(2023,12,1)
while s<=D.date(2024,11,30):
    print(s.isoformat()); s+=D.timedelta(days=1)"); do
  f="nya_${d}_12.csv"
  grep -q pressure_hPa "$f" 2>/dev/null && continue
  for try in 1 2 3 4 5 6 7 8; do
    curl -s -m 60 "https://weather.uwyo.edu/wsgi/sounding?datetime=${d}%2012:00:00&id=01004&type=TEXT:CSV&src=BUFR" -o "$f"
    grep -q "too busy" "$f" 2>/dev/null || break
    sleep 30
  done
  sleep 4
done
echo FETCH_DONE
