#!/bin/sh
# C++ sweeps of the thz-ntn absorption engine used by the paper.
set -e
O=$(cd "$(dirname "$0")" && pwd)/out
cd "$(dirname "$0")/../../../ns-3-dev"
B=./build/contrib/thz-ntn/examples/ns3.43-thz-ntn-absorption-sweep-optimized
$B --mode=inproc --pwvs=1,2,5,10,15 --out=$O/cpp_inproc.csv
$B --mode=inproc --profile=low-latitude --pwvs=20,30,45,50 --out=$O/cpp_inproc_lowlat.csv
$B --mode=lut --pwvs=none --out=$O/cpp_lut.csv
$B --mode=inproc --pwvs=none --out=$O/cpp_inproc_pub.csv
$B --mode=lutfile --lutFile=$O/lut_legacy_shipped.csv --pwvs=none --fStop=500 --out=$O/cpp_lut_legacy.csv
$B --mode=inproc --pwvs=none --steps=200 --elevs=10,30,90 --fStep=1 --out=$O/cpp_inproc_steps200.csv
for p in low-latitude midlatitude-summer midlatitude-winter high-latitude-winter; do
  $B --mode=inproc --pwvs=none --profile=$p --elevs=30 --fStep=1 --out=$O/cpp_profile_$p.csv
done
