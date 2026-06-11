import React, { useEffect, useMemo, useState } from "react";
import { motion } from "framer-motion";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Play, Pause, RotateCcw } from "lucide-react";

const quadruplets = [
  { p: 5, values: [5, 7, 11, 13], phase: "A-B-C-E", color: "border-blue-400", angle: 0 },
  { p: 11, values: [11, 13, 17, 19], phase: "C-E-A-B", color: "border-purple-400", angle: 42 },
  { p: 101, values: [101, 103, 107, 109], phase: "A-B-C-E", color: "border-amber-400", angle: 112 },
  { p: 191, values: [191, 193, 197, 199], phase: "C-E-A-B", color: "border-emerald-400", angle: 156 },
  { p: 821, values: [821, 823, 827, 829], phase: "A-B-C-E", color: "border-rose-400", angle: 218 },
  { p: 1481, values: [1481, 1483, 1487, 1489], phase: "A-B-C-E", color: "border-cyan-400", angle: 272 },
];

function xMap(n) {
  // Logarithmische Kompression: 5..1489 passen in eine einzige Raumzeit-Szene.
  return 70 + Math.log(n - 1) / Math.log(1490) * 680;
}

function yCurvature(n, t) {
  // Arithmetische Raumzeitkrümmung: Senken bei den Vierlingszentren.
  const centers = [9, 15, 105, 195, 825, 1485];
  let y = 250 + 18 * Math.sin(0.018 * n + t * 0.8);
  for (const c of centers) {
    const d = Math.log(Math.abs(n - c) + 1);
    y += 48 * Math.exp(-d * d / 1.2);
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

export default function PrimeQuadrupletSpacetimeAnimation() {
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
    <div className="min-h-screen w-full bg-slate-950 text-slate-100 p-6 flex flex-col gap-4">
      <div className="max-w-5xl mx-auto w-full">
        <div className="flex items-center justify-between gap-4 mb-4">
          <div>
            <h1 className="text-2xl md:text-4xl font-semibold tracking-tight">Primzahlvierlinge als gekrümmte arithmetische Raumzeit</h1>
            <p className="text-slate-300 mt-2 max-w-3xl">
              Einsteinische Lesart: Die Vierlinge entstehen nicht als gestörte perfekte Ellipsen, sondern als lokale Geodäten in einer durch Siebung und Restklassen gekrümmten Zahlraumzeit.
            </p>
          </div>
          <div className="flex gap-2">
            <Button variant="secondary" onClick={() => setRunning(!running)} className="rounded-2xl">
              {running ? <Pause className="w-4 h-4 mr-2" /> : <Play className="w-4 h-4 mr-2" />}
              {running ? "Pause" : "Start"}
            </Button>
            <Button variant="secondary" onClick={() => setT(0)} className="rounded-2xl">
              <RotateCcw className="w-4 h-4 mr-2" /> Reset
            </Button>
          </div>
        </div>

        <Card className="bg-slate-900/70 border-slate-700 rounded-2xl shadow-2xl">
          <CardContent className="p-4">
            <svg viewBox="0 0 840 520" className="w-full h-[560px] rounded-xl bg-gradient-to-b from-slate-950 to-slate-900">
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

              {/* grid */}
              {Array.from({ length: 12 }).map((_, i) => (
                <line key={`v-${i}`} x1={60 + i * 65} y1="65" x2={60 + i * 65} y2="455" stroke="white" strokeOpacity="0.06" />
              ))}
              {Array.from({ length: 8 }).map((_, i) => (
                <line key={`h-${i}`} x1="50" y1={90 + i * 50} x2="790" y2={90 + i * 50} stroke="white" strokeOpacity="0.06" />
              ))}

              {/* curved arithmetic spacetime */}
              <path d={curve} fill="none" stroke="white" strokeOpacity="0.34" strokeWidth="3" />
              <path d={curve} fill="none" stroke="white" strokeOpacity="0.08" strokeWidth="14" filter="url(#glow)" />

              {/* curvature wells */}
              {quadruplets.map((q, i) => {
                const cx = xMap(q.values[0] + 4);
                const cy = yCurvature(q.values[0] + 4, t);
                const visible = i <= activeIndex;
                const pulse = i === activeIndex ? 1 + 0.1 * Math.sin(localT * Math.PI * 2) : 1;
                return (
                  <g key={q.p} opacity={visible ? 1 : 0.18}>
                    <ellipse cx={cx} cy={cy + 8} rx={72 * pulse} ry={48 * pulse} fill="url(#well)" />
                    <path
                      d={orbitPath(cx, cy - 8, 58, 27, q.angle + 22 * Math.sin(t * 2 + i))}
                      fill="none"
                      stroke="white"
                      strokeOpacity={visible ? 0.75 : 0.18}
                      strokeWidth="2.5"
                    />
                    <path
                      d={orbitPath(cx, cy - 8, 58, 27, q.angle + 22 * Math.sin(t * 2 + i))}
                      fill="none"
                      stroke="white"
                      strokeOpacity={visible ? 0.18 : 0.04}
                      strokeWidth="10"
                      filter="url(#glow)"
                    />
                    <line
                      x1={cx}
                      y1={cy - 8}
                      x2={cx + 58 * Math.cos(((q.angle + 22 * Math.sin(t * 2 + i)) * Math.PI) / 180)}
                      y2={cy - 8 + 58 * Math.sin(((q.angle + 22 * Math.sin(t * 2 + i)) * Math.PI) / 180)}
                      stroke="white"
                      strokeOpacity={visible ? 0.8 : 0.16}
                      strokeWidth="2"
                    />
                    <text x={cx - 58} y={cy - 62} fill="white" opacity={visible ? 0.9 : 0.25} fontSize="14">
                      ({q.values.join(", ")})
                    </text>
                    <text x={cx - 42} y={cy - 43} fill="white" opacity={visible ? 0.55 : 0.18} fontSize="12">
                      Phase: {q.phase}
                    </text>
                    {q.values.map((n, j) => (
                      <g key={n}>
                        <circle cx={xMap(n)} cy={yCurvature(n, t)} r={visible ? 5.2 : 3} fill="white" opacity={visible ? 0.95 : 0.24} />
                        <text x={xMap(n) - 7} y={yCurvature(n, t) + 24} fill="white" opacity={visible ? 0.72 : 0.2} fontSize="11">
                          {n}
                        </text>
                      </g>
                    ))}
                  </g>
                );
              })}

              {/* propagation light cone / emergence front */}
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

              <text x="55" y="488" fill="white" opacity="0.65" fontSize="13">Zahlraum, logarithmisch komprimiert</text>
              <text x="500" y="488" fill="white" opacity="0.65" fontSize="13">Zeit / Entstehungsfolge der sechs Vierlings-Geodäten</text>
            </svg>
          </CardContent>
        </Card>

        <div className="grid md:grid-cols-3 gap-3 mt-4">
          {quadruplets.map((q, i) => (
            <Card key={q.p} className={`bg-slate-900/70 border ${i === activeIndex ? "border-white" : "border-slate-700"} rounded-2xl`}>
              <CardContent className="p-4">
                <div className="text-sm text-slate-400">Vierling {i + 1}</div>
                <div className="text-xl font-semibold mt-1">({q.values.join(", ")})</div>
                <div className="text-sm text-slate-300 mt-2">EABC-Phase: {q.phase}</div>
                <div className="text-sm text-slate-400 mt-2">
                  Perihelwinkel: ca. {Math.round(q.angle + 22 * Math.sin(t * 2 + i))}°
                </div>
              </CardContent>
            </Card>
          ))}
        </div>

        <Card className="bg-slate-900/70 border-slate-700 rounded-2xl mt-4">
          <CardContent className="p-4 text-slate-300 leading-relaxed">
            <p>
              Modellidee: Jeder Primzahlvierling erzeugt eine lokale Krümmungsmulde. Die elliptische Bahn ist nicht primär eine gestörte Newton-Bahn, sondern eine Geodäte in dieser Mulde. Die sichtbare Periheldrehung ist daher ein Messsignal der arithmetischen Raumzeitkrümmung: vom Urvierling (5,7,11,13) über den überlappenden Vierling (11,13,17,19) bis zu den weiteren Krümmungsinseln (101,103,107,109), (191,193,197,199), (821,823,827,829) und (1481,1483,1487,1489).
            </p>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
