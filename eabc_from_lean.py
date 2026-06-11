from __future__ import annotations

import json
from enum import Enum
from math import isqrt
from typing import Any


class EClass(Enum):
    E = "E"
    A = "A"
    B = "B"
    C = "C"


def residue(eclass: EClass) -> int:
    mapping = {
        EClass.E: 1,
        EClass.A: 5,
        EClass.B: 7,
        EClass.C: 11,
    }
    return mapping[eclass]


def class_of(n: int) -> EClass | None:
    mod = n % 12
    mapping = {
        1: EClass.E,
        5: EClass.A,
        7: EClass.B,
        11: EClass.C,
    }
    return mapping.get(mod)


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    limit = isqrt(n)
    for d in range(3, limit + 1, 2):
        if n % d == 0:
            return False
    return True


def is_prime_quadruplet(p: int) -> bool:
    return _is_prime(p) and _is_prime(p + 2) and _is_prime(p + 6) and _is_prime(p + 8)


def q(p: int) -> list[int]:
    return [p, p + 2, p + 6, p + 8]


def center(p: int) -> int:
    return p + 4


GAP_SIGNATURE = [2, 4, 2]


class Chirality(Enum):
    ABCE = "ABCE"
    CEAB = "CEAB"


def chirality_order(chirality: Chirality) -> list[EClass]:
    if chirality is Chirality.ABCE:
        return [EClass.A, EClass.B, EClass.C, EClass.E]
    return [EClass.C, EClass.E, EClass.A, EClass.B]


def t(x: EClass) -> EClass:
    mapping = {
        EClass.E: EClass.A,
        EClass.A: EClass.B,
        EClass.B: EClass.C,
        EClass.C: EClass.E,
    }
    return mapping[x]


def t4(x: EClass) -> EClass:
    return t(t(t(t(x))))


def definitions_overview() -> dict[str, Any]:
    return {
        "enums": {
            "EClass": {
                "description": "Vier EABC-Familien als Klassen.",
                "values": [member.value for member in EClass],
            },
            "Chirality": {
                "description": "Die zwei chiralen Reihenfolgen.",
                "values": [member.value for member in Chirality],
            },
        },
        "constants": {
            "GAP_SIGNATURE": {
                "description": "Interne Gap-Signatur des Vierlings Q(p).",
                "value": GAP_SIGNATURE,
            }
        },
        "functions": {
            "residue": "Ordnet einer EClass den Modulo-12-Rest zu.",
            "class_of": "Liefert die EClass einer Zahl modulo 12 oder None.",
            "is_prime_quadruplet": "Prueft, ob p, p+2, p+6, p+8 alle prim sind.",
            "q": "Gibt den Vierling [p, p+2, p+6, p+8] zurueck.",
            "center": "Gibt den Mittelpunkt p+4 des Vierlings zurueck.",
            "chirality_order": "Liefert die Klassenfolge fuer eine Chiralitaet.",
            "t": "Rotiert E->A->B->C->E.",
            "t4": "Vierfache Anwendung von t; ergibt wieder das Ausgangselement.",
        },
    }


if __name__ == "__main__":
    print(json.dumps(definitions_overview(), indent=2, ensure_ascii=False))
