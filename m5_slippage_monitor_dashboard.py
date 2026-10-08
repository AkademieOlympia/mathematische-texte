import zlib

import numpy as np
import pandas as pd
import streamlit as st

REGIMES = {
    "Normal (wenig Jitter)": {
        "slippage_bps": 0.08,
        "legacy_multiple": 11.2,
        "p99_ns": 350,
        "fill_rate": 99.99,
        "legacy_fill_rate": 89.01,
        "jitter_ns": 12.0,
        "slip_sigma": 0.015,
    },
    "Hohe Volatilität (L2-Lücken)": {
        "slippage_bps": 0.15,
        "legacy_multiple": 11.2,
        "p99_ns": 420,
        "fill_rate": 99.95,
        "legacy_fill_rate": 86.40,
        "jitter_ns": 22.0,
        "slip_sigma": 0.03,
    },
    "Extremer Stress (Flash Crash)": {
        "slippage_bps": 0.28,
        "legacy_multiple": 11.2,
        "p99_ns": 580,
        "fill_rate": 99.85,
        "legacy_fill_rate": 78.20,
        "jitter_ns": 40.0,
        "slip_sigma": 0.06,
    },
}

TRANCHE_COUNT = 240


def market_slippage_bps(base_bps: float, volume_m: float) -> float:
    """Square-root impact: larger clips walk the book."""
    return base_bps * (1.0 + 0.35 * np.sqrt(volume_m / 50.0))


def execution_stats(volume_m: float, regime: str, tolerance_bps: float) -> dict:
    spec = REGIMES[regime]
    raw_bps = market_slippage_bps(spec["slippage_bps"], volume_m)
    overshoot_bps = max(0.0, raw_bps - tolerance_bps)
    clipped = overshoot_bps > 0.0
    effective_bps = min(raw_bps, tolerance_bps)
    # Capping the price costs a little fill, not the whole order.
    fill_rate = max(95.0, spec["fill_rate"] - overshoot_bps * 8.0)
    legacy_bps = market_slippage_bps(spec["slippage_bps"] * spec["legacy_multiple"], volume_m)
    notional = volume_m * 1_000_000
    cost = notional * (effective_bps / 10_000.0)
    legacy_cost = notional * (legacy_bps / 10_000.0)
    return {
        "raw_bps": raw_bps,
        "effective_bps": effective_bps,
        "legacy_bps": legacy_bps,
        "clipped": clipped,
        "p99_ns": spec["p99_ns"],
        "fill_rate": fill_rate,
        "legacy_fill_rate": spec["legacy_fill_rate"],
        "cost": cost,
        "legacy_cost": legacy_cost,
        "savings": legacy_cost - cost,
    }


def _trail(end: float, noise: float, seed: int, n: int = 18) -> list[float]:
    rng = np.random.default_rng(seed)
    walk = np.cumsum(rng.normal(0.0, noise, n))
    walk = walk - walk[-1] + end
    return walk.tolist()


def _seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode()
    return zlib.adler32(payload)


@st.cache_data
def build_monitor(volume_m: float, regime: str, tolerance_bps: float) -> dict:
    stats = execution_stats(volume_m, regime, tolerance_bps)
    spec = REGIMES[regime]
    rng = np.random.default_rng(_seed(volume_m, regime, tolerance_bps))
    index = np.arange(1, TRANCHE_COUNT + 1)
    slippage = np.clip(rng.normal(stats["raw_bps"], spec["slip_sigma"], TRANCHE_COUNT), 0.01, None)
    latency = np.clip(rng.normal(spec["p99_ns"], spec["jitter_ns"], TRANCHE_COUNT), 80.0, None)
    capped = slippage > tolerance_bps
    effective = np.minimum(slippage, tolerance_bps)
    notional = volume_m * 1_000_000 / TRANCHE_COUNT
    profile = pd.DataFrame(
        {
            "Tranche": index,
            "Latenz (ns)": latency,
            "p99-Ziel (ns)": np.full(TRANCHE_COUNT, spec["p99_ns"]),
            "Markt (bps)": slippage,
            "Nach Guardrail (bps)": effective,
        }
    )
    tranches = pd.DataFrame(
        {
            "Tranche": index,
            "Latenz (ns)": np.rint(latency).astype(int),
            "Markt-Slippage (bps)": np.round(slippage, 3),
            "Effektive Slippage (bps)": np.round(effective, 3),
            "Kosten ($)": np.round(notional * (effective / 10_000.0), 2),
            "Status": np.where(capped, "Gekappt", "Im Limit"),
        }
    )
    seed = _seed("trail", volume_m, regime, tolerance_bps)
    trails = {
        "slippage": _trail(stats["effective_bps"], max(stats["effective_bps"] * 0.04, 0.005), seed),
        "latency": _trail(float(stats["p99_ns"]), 8.0, seed + 1),
        "fill": _trail(stats["fill_rate"], 0.08, seed + 2),
        "savings": _trail(stats["savings"], max(stats["savings"] * 0.03, 1.0), seed + 3),
    }
    return {"stats": stats, "profile": profile, "tranches": tranches, "trails": trails}


st.set_page_config(
    page_title="Order-Fill und Slippage",
    page_icon=":material/monitoring:",
    layout="wide",
)

st.title("Order-Fill und Slippage")
st.caption(
    "Ausführungsmonitor für die M5-Topos-Guardrail. "
    "Liegt die Markt-Slippage über der Toleranz, wird der Preis gekappt. Die Fill-Rate gibt nur leicht nach."
)

with st.sidebar:
    st.header("Ausführung")
    order_volume_m = st.slider(
        "Handelsvolumen (Mio. $)",
        min_value=1.0,
        max_value=500.0,
        value=50.0,
        step=5.0,
    )
    volatility_regime = st.selectbox("Marktregime", list(REGIMES))
    max_slippage_tolerance_bps = st.slider(
        "Maximale Slippage-Toleranz (bps)",
        min_value=0.05,
        max_value=1.0,
        value=0.20,
        step=0.01,
    )
    st.caption("Die Toleranz ist die Kappungsgrenze der Guardrail, nicht nur eine Anzeige.")

view = build_monitor(order_volume_m, volatility_regime, max_slippage_tolerance_bps)
stats = view["stats"]

with st.container(horizontal=True):
    if stats["clipped"]:
        st.badge("Guardrail greift", icon=":material/shield:", color="orange")
    else:
        st.badge("Innerhalb der Toleranz", icon=":material/check_circle:", color="green")
    st.badge(volatility_regime, icon=":material/show_chart:", color="blue")

saved_share = stats["savings"] / stats["legacy_cost"] if stats["legacy_cost"] else 0.0

top_left, top_right = st.columns(2, border=True)
with top_left:
    st.metric(
        "Effektive Slippage (bps)",
        stats["effective_bps"],
        delta=stats["effective_bps"] - stats["legacy_bps"],
        delta_color="inverse",
        delta_description="gegen Legacy",
        format="%.2f",
        chart_data=view["trails"]["slippage"],
        chart_type="area",
        help="Markt-Slippage nach Square-Root-Impact, gekappt bei der Toleranz.",
    )
with top_right:
    st.metric(
        "p99-Latenz (ns)",
        stats["p99_ns"],
        delta=stats["p99_ns"] - 1000,
        delta_color="inverse",
        delta_description="gegen 1 µs",
        format="%d",
        chart_data=view["trails"]["latency"],
        chart_type="line",
    )

bottom_left, bottom_right = st.columns(2, border=True)
with bottom_left:
    st.metric(
        "Fill-Rate (%)",
        stats["fill_rate"],
        delta=stats["fill_rate"] - stats["legacy_fill_rate"],
        delta_description="Prozentpunkte gegen Legacy",
        format="%.2f",
        chart_data=view["trails"]["fill"],
        chart_type="line",
        help="Pro Basispunkt über der Toleranz sinkt die Fill-Rate um 8 Prozentpunkte, mindestens auf 95 %.",
    )
with bottom_right:
    st.metric(
        "Eingespartes Kapital",
        stats["savings"],
        delta=f"{saved_share:.0%}",
        format="dollar",
        delta_description="der Legacy-Kosten",
        chart_data=view["trails"]["savings"],
        chart_type="area",
    )

left, right = st.columns(2)
with left:
    with st.container(border=True):
        st.subheader("Fill-Latenz")
        st.line_chart(
            view["profile"],
            x="Tranche",
            y=["Latenz (ns)", "p99-Ziel (ns)"],
            color=["green", "orange"],
            height=320,
            alt="Latenz je Tranche gegen das p99-Ziel",
        )
with right:
    with st.container(border=True):
        st.subheader("Slippage je Tranche")
        st.line_chart(
            view["profile"],
            x="Tranche",
            y=["Markt (bps)", "Nach Guardrail (bps)"],
            color=["orange", "green"],
            height=320,
            alt="Markt-Slippage und gekappte Slippage je Tranche",
        )
        st.caption(
            f"Markt {stats['raw_bps']:.2f} bps, effektiv {stats['effective_bps']:.2f} bps, "
            f"Legacy {stats['legacy_bps']:.2f} bps."
        )

with st.container(border=True):
    st.subheader("Tranchen")
    capped_count = int((view["tranches"]["Status"] == "Gekappt").sum())
    st.caption(
        f"{capped_count} von {TRANCHE_COUNT} Tranchen über {max_slippage_tolerance_bps:.2f} bps, "
        "ausgeführt zum Limit."
    )
    st.dataframe(
        view["tranches"],
        hide_index=True,
        height=360,
        column_config={
            "Latenz (ns)": st.column_config.NumberColumn(format="%d"),
            "Markt-Slippage (bps)": st.column_config.NumberColumn(format="%.3f"),
            "Effektive Slippage (bps)": st.column_config.NumberColumn(format="%.3f"),
            "Kosten ($)": st.column_config.NumberColumn(format="$%.2f"),
        },
    )
