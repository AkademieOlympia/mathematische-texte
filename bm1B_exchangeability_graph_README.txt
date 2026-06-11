Exchangeability graphs (within-family substitution):
- A-A edges: two A-primes that co-occur under the same (B,C) context (enough support + alternatives).
- B-B edges: two B-primes that co-occur under the same (A,C) context.
- C-C edges: two C-primes that co-occur under the same (A,B) context.

Edge weight = sum over contexts of (entropy_norm(context) * total(context)) / C(k,2),
where k = number of distinct primes in that context.
Edge label in DOT: "weight/contexts" (contexts = how many qualifying contexts contributed).

Render:
  dot -Tsvg bm1B_exchangeability_graph_AA.dot -o AA.svg
  dot -Tsvg bm1B_exchangeability_graph_BB.dot -o BB.svg
  dot -Tsvg bm1B_exchangeability_graph_CC.dot -o CC.svg
