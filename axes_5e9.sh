#!/usr/bin/env bash
# Sammelt eine CSV (Header + 8 Sektormuster) per ./heeger
#
# Warum nicht “for P in …; do set -- $P; ./heeger … $1 $2” in zsh?
#   In zsh ist SH_WORD_SPLIT standardmäßig AUS: $P bleibt EIN Wort.
#   Dann fehlt $2 oder es landen falsche argv → u.a. “R1==R2”-Fehler.
# Lösungen:
#   (1) Dieses Skript mit bash:     bash ./axes_5e9.sh
#   (1b) Oder:                        zsh ./axes_5e9.zsh
#   (2) In zsh explizit splitten:   set -- ${=P}    # oder: set -- ${(s: :)P}
#   (3) Ohne set:   read r1 r2 <<< "$P"   (mit Anführungszeichen testen)
#
# Vollständiger 5e9-Lauf dauert pro Aufruf viele Minuten (Sieb).

set -euo pipefail
cd "$(dirname "$0")"

N="${1:-5000000000}"
OUT="${2:-axes_5e9.csv}"
HEEGER="${3:-./heeger}"

if [[ ! -x "$HEEGER" ]]; then
  echo "Nicht ausführbar: $HEEGER (in diesem Ordner bauen: g++ -o heeger Heeger.cpp)" >&2
  exit 1
fi

echo "N,R1,R2,dist,dominant_mode,twin_count,sector_count,delta_ons,p_ons,d16,p16,d25,p25,d34,p34,purity,d_eff,purity_over_iso" > "$OUT"

while read -r r1 r2; do
  [[ -z "${r1:-}" ]] && continue
  "$HEEGER" "$N" "$r1" "$r2" | grep -A2 CSV_SUMMARY | tail -n 1
done >> "$OUT" <<'PAIRS'
11 41
17 47
59 29
11 29
59 17
59 41
11 17
41 47
PAIRS

echo "Wrote $OUT ($(wc -l < "$OUT") lines)"
