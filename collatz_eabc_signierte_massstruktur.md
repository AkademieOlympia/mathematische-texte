# EABC: Signierte Maßstruktur auf dem Zustandsgraphen (jenseits Wigner)

**Status:** Modellabbildung + diskrete Numerik  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-euklidische-hebung`  
**Tao-Labels:** Definition | Analogie | Hypothese | Modellabbildung | offene Frage

**Primäre Geometrie:** `collatz_eabc_uebergangsraum.md` — Übergangsraum, $S^1$, Fluss/Wilson-Neulesung.

**Querverweise:**
- `collatz_eabc_wigner_analog.md` — historische 7-Abschnitte-Wigner-Analogie
- `collatz_eabc_hodge_eabc.py` — `signed_measure_graph`, `orientation_information_test`
- `collatz_eabc_wigner_field.py` — $W_E(i,j;N)$, Informationsüberschuss-Test

---

## 0. Boxed Leitfrage

$$\boxed{\;\text{Signierte Maßstruktur auf } G_E=(V,E,w)\text{ — Information in Orientierung, nicht in }p_n.\;}$$

---

## 1. Lokales Zustandsfeld $W_E(N)$

$$W_E(N) = \#\mathrm{ABCE} - \#\mathrm{CEAB},\qquad
S_W(N) = \frac{W_E(N)}{N_+^{(4)}+N_-^{(4)}} \in [-1,1].$$

**Neulesung (Übergangsraum):** $S_W \approx \tanh\Theta_E$ — normierte Flussdichte, nicht Quasi-Wahrscheinlichkeit (`collatz_eabc_uebergangsraum.md` §4).

Vorzeichendomänen von $W_E^{\mathrm{cum}}(n)$: **chirale Grenze** bei $W_E\approx 0$ (nicht Maxima).

---

## 2. Gerichteter gewichteter Graph $G_E=(V,E,w)$

$V=\{E,A,B,C\}$, $w_{ij}$ aus $W$-Matrix — **Netzwerktheorie**, nicht klassische Zahlentheorie.

**Numerik:** `signed_measure_graph` in `collatz_eabc_hodge_eabc.py`.

---

## 3. Feynman-/Propagator-Lesart

Übergänge $E\!\to\!A\!\to\!\cdots$, nicht Primzahlen $p$ — Pfadintegral auf $G_E$.

**Label:** Propagator-Bild = **Analogie**.

---

## 4. Near-zero: chirale Grenze

$W_E\approx 0$: Grenze zwischen Orientierungsdomänen (wie $\lambda_k\approx 0$ in Spektralgeometrie).

---

## 5. $L = D - W$, Spektrum $\lambda_1\le\cdots\le\lambda_4$

Near-zero-Moden = Interferenz-/Orientierungsregionen.

**Numerik:** `laplacian_from_W`.

---

## 6. Holonomie $\neq$ Wahrscheinlichkeit

$C_E = \oint\omega$; Wigner = signierter Phasenraum; EABC = signierte Zyklen im Restklassenraum.

---

## 7. OFFEN: Arithmetische Wigner-Negativität

$$\boxed{\;\text{Trägt die Vorzeichenstruktur von }W_E(i,j;N)\text{ Information jenseits der Marginalen }(N_++N_-,\,S_E)?\;}$$

**Test:** `orientation_information_test` — $(N_+,N_-)$ algebraisch aus $(N_{\mathrm{tot}},S_E)$ rekonstruierbar; **Kantenfeld** $W_E(i,j)$ aus Marginalen **nicht** (`information_excess_test` in `collatz_eabc_wigner_field.py`).

**Label:** Arithmetische Wigner-Negativität = **Hypothese** + **Experiment**.

---

## 8. Artefakte

```bash
python3 collatz_eabc_wigner_field.py --max-p 100000
python3 collatz_eabc_hodge_eabc.py --max-p 1000000
pytest tests/test_eabc_wigner_field.py tests/test_eabc_hodge_eabc.py -q
```
