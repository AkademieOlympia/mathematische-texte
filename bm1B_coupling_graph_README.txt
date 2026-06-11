Graphviz Rendering:
  dot -Tpng bm1B_coupling_graph_top50.dot -o bm1B_coupling_graph_top50.png
  dot -Tsvg bm1B_coupling_graph_top50.dot -o bm1B_coupling_graph_top50.svg

Edge label format: "lift_1B / count_1B"
Edge style:
  solid   = in_both (stable across 100M and 1B)
  dashed  = only_1B (new at 1B scale)
  dotted  = only_100M (was top at 100M only)

Node colors:
  lightblue         = A_pair
  lightgreen        = B_pair
  lightgoldenrod1   = C_pair
