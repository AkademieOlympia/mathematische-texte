import Mathlib.Geometry.Euclidean.Sphere.Ptolemy

open EuclideanGeometry

section

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {P : Type*} [MetricSpace P] [NormedAddTorsor V P]
variable {a b c d p : P}

example
    (h : Cospherical ({a, b, c, d} : Set P))
    (hapc : angle a p c = Real.pi)
    (hbpd : angle b p d = Real.pi) :
    dist a b * dist c d + dist b c * dist d a = dist a c * dist b d := by
  simpa using
    EuclideanGeometry.mul_dist_add_mul_dist_eq_mul_dist_of_cospherical
      h hapc hbpd

end
