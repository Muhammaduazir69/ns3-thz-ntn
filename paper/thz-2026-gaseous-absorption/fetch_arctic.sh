#!/bin/sh
# Six further Arctic stations, 12 UTC, University of Wyoming archive: every second day from
# December 2023 to February 2024 and every sixth day from March to November 2024.
# Usage: fetch_arctic.sh <dir>. Files already holding a sounding are kept; a "busy" reply is retried.
cd "$1" || exit 1
DATES=$(python3 -c "
import datetime as D
s=D.date(2023,12,1); i=0
while s<=D.date(2024,2,29):
    print(s.isoformat()); s+=D.timedelta(days=2)
s=D.date(2024,3,1)
while s<=D.date(2024,11,30):
    print(s.isoformat()); s+=D.timedelta(days=6)")
for d in $DATES; do
  for st in 70026 71917 71082 04320 01028 02836; do
    f="s${st}_${d}_12.csv"
    grep -q pressure_hPa "$f" 2>/dev/null && continue
    for try in 1 2 3 4 5; do
      curl -s -m 60 "https://weather.uwyo.edu/wsgi/sounding?datetime=${d}%2012:00:00&id=${st}&type=TEXT:CSV&src=BUFR" -o "$f"
      grep -q "too busy" "$f" 2>/dev/null || break
      sleep 30
    done
    sleep 3
  done
done
echo FETCH_DONE
