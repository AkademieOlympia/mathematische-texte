import React, { useEffect, useMemo, useState } from "react";
import { motion } from "framer-motion";
import { Play, Pause, RotateCcw } from "lucide-react";

const quadruplets = [
  { p: 5, values: [5, 7, 11, 13], phase: "A-B-C-E", angle: 0 },
  { p: 11, values: [11, 13, 17, 19], phase: "C-E-A-B", angle: 42 },
  { p: 101, values: [101, 103, 107, 109], phase: "A-B-C-E", angle: 112 },
  { p: 191, values: [191, 193, 197, 199], phase: "C-E-A-B", angle: 156 },
  { p: 821, values: [821, 823, 827, 829], phase: "A-B-C-E", angle: 218 },
  { p: 1481, values: [1481, 1483, 1487, 1489], phase: "A-B-C-E", angle: 272 },
];

function Card({ children, className = "" }) {
  return <div className={`card ${className}`}>{children}</div>;
}

function CardContent({ children, className = "" }) {
  return <div className={`card-content ${className}`}>{children}</div>;
}

function Button({ children, onClick }) {
  return (
    <button type="button" className="btn" onClick={onClick}>
      {children}
    </button>
  );
}

function xMap(n) {
  return 70 + (Math.log(n - 1) / Math.log(1490)) * 680;
}

function yCurvature(n, t) {
  const centers = [9, 15, 105, 195, 825, 1485];
  let y = 250 + 18 * Math.sin(0.018 * n + t * 0.8);
  for (const c of centers) {
    const d = Math.log(Math.abs(n - c) + 1);
    y += 48 * Math.exp((-d * d) / 1.2);
  }
  return y;
}

function orbitPath(cx, cy, rx, ry, angleDeg) {
  const a = (angleDeg * Math.PI) / 180;
  const pts = [];
  for (let i = 0; i <= 160; i++) {
    const u = (i / 160) * Math.PI * 2;
    const x = rx * Math.cos(u);
    const y = ry * Math.sin(u);
    const xr = cx + x * Math.cos(a) - y * Math.sin(a);
    const yr = cy + x * Math.sin(a) + y * Math.cos(a);
    pts.push(`${xr.toFixed(1)},${yr.toFixed(1)}`);
  }
  return `M ${pts[0]} L ${pts.slice(1).join(" L ")} Z`;
}

export default function App() {
  const [running, setRunning] = useState(true);
  const [t, setT] = useState(0);

  useEffect(() => {
    if (!running) return;
    const id = setInterval(() => setT((v) => (v + 0.014) % 6), 32);
    return () => clearInterval(id);
  }, [running]);

  const activeIndex = Math.min(5, Math.floor(t));
  const localT = t - activeIndex;

  const curve = useMemo(() => {
    const pts = [];
    for (let n = 2; n <= 1500; n += n < 220 ? 2 : 10) {
      pts.push(`${xMap(n).toFixed(1)},${yCurvature(n, t).toFixed(1)}`);
    }
    return `M ${pts.join(" L ")}`;
  }, [t]);

  return (
    <div className="page">
      <div className="container">
        <div className="row">
          <div>
            <h1>
              Primzahlvierlinge als gekruemmte arithmetische Raumzeit
            </h1>
            <p>
              Einsteinische Lesart: Die Vierlinge entstehen nicht als gestoerte
              perfekte Ellipsen, sondern als lokale Geodaeten in einer durch
              Siebung und Restklassen gekruemmten Zahlraumzeit.
            </p>
          </div>
          <div className="buttons">
            <Button onClick={() => setRunning(!running)}>
              {running ? <Pause size={16} /> : <Play size={16} />}
              {running ? "Pause" : "Start"}
            </Button>
            <Button onClick={() => setT(0)}>
              <RotateCcw size={16} />
              Reset
            </Button>
          </div>
        </div>

        <Card>
          <CardContent>
            <svg
              viewBox="0 0 840 520"
              style={{
                width: "100%",
                height: "560px",
                borderRadius: "0.75rem",
                background:
                  "linear-gradient(to bottom, rgb(2,6,23), rgb(15,23,42))",
              }}
            >
              <defs>
                <radialGradient id="well" cx="50%" cy="50%" r="50%">
                  <stop offset="0%" stopColor="white" stopOpacity="0.16" />
                  <stop offset="70%" stopColor="white" stopOpacity="0.04" />
                  <stop offset="100%" stopColor="white" stopOpacity="0" />
                </radialGradient>
                <filter id="glow">
                  <feGaussianBlur stdDeviation="4" result="coloredBlur" />
                  <feMerge>
                    <feMergeNode in="coloredBlur" />
                    <feMergeNode in="SourceGraphic" />
                  </feMerge>
                </filter>
              </defs>

              {Array.from({ length: 12 }).map((_, i) => (
                <line
                  key={`v-${i}`}
                  x1={60 + i * 65}
                  y1="65"
                  x2={60 + i * 65}
                  y2="455"
                  stroke="white"
                  strokeOpacity="0.06"
                />
              ))}
              {Array.from({ length: 8 }).map((_, i) => (
                <line
                  key={`h-${i}`}
                  x1="50"
                  y1={90 + i * 50}
                  x2="790"
                  y2={90 + i * 50}
                  stroke="white"
                  strokeOpacity="0.06"
                />
              ))}

              <path
                d={curve}
                fill="none"
                stroke="white"
                strokeOpacity="0.34"
                strokeWidth="3"
              />
              <path
                d={curve}
                fill="none"
                stroke="white"
                strokeOpacity="0.08"
                strokeWidth="14"
                filter="url(#glow)"
              />

              {quadruplets.map((q, i) => {
                const cx = xMap(q.values[0] + 4);
                const cy = yCurvature(q.values[0] + 4, t);
                const visible = i <= activeIndex;
                const pulse =
                  i === activeIndex ? 1 + 0.1 * Math.sin(localT * Math.PI * 2) : 1;
                return (
                  <g key={q.p} opacity={visible ? 1 : 0.18}>
                    <ellipse
                      cx={cx}
                      cy={cy + 8}
                      rx={72 * pulse}
                      ry={48 * pulse}
                      fill="url(#well)"
                    />
                    <path
                      d={orbitPath(
                        cx,
                        cy - 8,
                        58,
                        27,
                        q.angle + 22 * Math.sin(t * 2 + i),
                      )}
                      fill="none"
                      stroke="white"
                      strokeOpacity={visible ? 0.75 : 0.18}
                      strokeWidth="2.5"
                    />
                    <path
                      d={orbitPath(
                        cx,
                        cy - 8,
                        58,
                        27,
                        q.angle + 22 * Math.sin(t * 2 + i),
                      )}
                      fill="none"
                      stroke="white"
                      strokeOpacity={visible ? 0.18 : 0.04}
                      strokeWidth="10"
                      filter="url(#glow)"
                    />
                    <line
                      x1={cx}
                      y1={cy - 8}
                      x2={
                        cx +
                        58 *
                          Math.cos(
                            ((q.angle + 22 * Math.sin(t * 2 + i)) * Math.PI) /
                              180,
                          )
                      }
                      y2={
                        cy -
                        8 +
                        58 *
                          Math.sin(
                            ((q.angle + 22 * Math.sin(t * 2 + i)) * Math.PI) /
                              180,
                          )
                      }
                      stroke="white"
                      strokeOpacity={visible ? 0.8 : 0.16}
                      strokeWidth="2"
                    />
                    <text
                      x={cx - 58}
                      y={cy - 62}
                      fill="white"
                      opacity={visible ? 0.9 : 0.25}
                      fontSize="14"
                    >
                      ({q.values.join(", ")})
                    </text>
                    <text
                      x={cx - 42}
                      y={cy - 43}
                      fill="white"
                      opacity={visible ? 0.55 : 0.18}
                      fontSize="12"
                    >
                      Phase: {q.phase}
                    </text>
                    {q.values.map((n) => (
                      <g key={n}>
                        <circle
                          cx={xMap(n)}
                          cy={yCurvature(n, t)}
                          r={visible ? 5.2 : 3}
                          fill="white"
                          opacity={visible ? 0.95 : 0.24}
                        />
                        <text
                          x={xMap(n) - 7}
                          y={yCurvature(n, t) + 24}
                          fill="white"
                          opacity={visible ? 0.72 : 0.2}
                          fontSize="11"
                        >
                          {n}
                        </text>
                      </g>
                    ))}
                  </g>
                );
              })}

              <motion.line
                x1={xMap(5 + t * 250)}
                x2={xMap(5 + t * 250)}
                y1="60"
                y2="460"
                stroke="white"
                strokeOpacity="0.25"
                strokeDasharray="6 7"
                animate={{ opacity: [0.15, 0.45, 0.15] }}
                transition={{ duration: 1.8, repeat: Infinity }}
              />

              <text x="55" y="488" fill="white" opacity="0.65" fontSize="13">
                Zahlraum, logarithmisch komprimiert
              </text>
              <text x="500" y="488" fill="white" opacity="0.65" fontSize="13">
                Zeit / Entstehungsfolge der sechs Vierlings-Geodaeten
              </text>
            </svg>
          </CardContent>
        </Card>

        <div className="grid">
          {quadruplets.map((q, i) => (
            <Card key={q.p}>
              <CardContent>
                <div style={{ fontSize: "0.9rem", color: "#94a3b8" }}>
                  Vierling {i + 1}
                </div>
                <div style={{ fontSize: "1.2rem", fontWeight: 600, marginTop: "0.2rem" }}>
                  ({q.values.join(", ")})
                </div>
                <div style={{ fontSize: "0.9rem", marginTop: "0.5rem" }}>
                  EABC-Phase: {q.phase}
                </div>
                <div style={{ fontSize: "0.85rem", marginTop: "0.45rem", color: "#94a3b8" }}>
                  Perihelwinkel: ca.{" "}
                  {Math.round(q.angle + 22 * Math.sin(t * 2 + i))}°
                </div>
              </CardContent>
            </Card>
          ))}
        </div>

        <Card>
          <CardContent>
            <p>
              Modellidee: Jeder Primzahlvierling erzeugt eine lokale
              Kruemmungsmulde. Die elliptische Bahn ist nicht primaer eine
              gestoerte Newton-Bahn, sondern eine Geodaete in dieser Mulde.
              Die sichtbare Periheldrehung ist daher ein Messsignal der
              arithmetischen Raumzeitkruemmung: vom Urvierling (5,7,11,13)
              ueber den ueberlappenden Vierling (11,13,17,19) bis zu den
              weiteren Kruemmungsinseln (101,103,107,109), (191,193,197,199),
              (821,823,827,829) und (1481,1483,1487,1489).
            </p>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
