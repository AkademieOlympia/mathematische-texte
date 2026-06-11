#!/usr/bin/env zsh
# Gegenstück zur bash-Version: in zsh *muss* $P in zwei Felder zerlegt werden.
#   set -- $P     → oft EIN Wort  → fehlerhafte ./heeger-Argumente
#   set -- ${=P}  → 11, 41 getrennt (SH_WORD_SPLIT für diese Expansion)

emulate -L zsh
set -e
cd -- "$(dirname "$0")"

N="${1:-5000000000}"
OUT="${2:-axes_5e9.csv}"
H="${3:-./heeger}"
[[ -x $H ]] || { print -u2 "fehlt/ausfuehrbar: $H"; exit 1 }

print -r "N,R1,R2,dist,dominant_mode,twin_count,sector_count,delta_ons,p_ons,d16,p16,d25,p25,d34,p34,purity,d_eff,purity_over_iso" > "$OUT"

for P in "11 41" "17 47" "59 29" "11 29" "59 17" "59 41" "11 17" "41 47"; do
  set -- ${=P}
  "$H" "$N" "$1" "$2" | grep -A2 CSV_SUMMARY | tail -n 1
done >> "$OUT"

print "OK -> $OUT ($(wc -l < "$OUT" | tr -d ' ') Zeilen, inkl. Header)"
