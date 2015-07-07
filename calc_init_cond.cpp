#include "calc_init_cond.hpp"

vector<ComputationalCell> calc_init_cond
(const Tessellation& tess)
{
  ComputationalCell c;
  c.density = 1;
  c.pressure = 1;
  c.velocity = Vector2D(0,0);
  return vector<ComputationalCell>
    (static_cast<size_t>(tess.GetPointNo()),
     c);
}
