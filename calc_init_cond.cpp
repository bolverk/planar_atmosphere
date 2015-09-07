#include <cmath>
#include "calc_init_cond.hpp"

vector<ComputationalCell> calc_init_cond
(const Tessellation& tess,
 const double /*g*/, // adiabatic index
 const double a, // acceleration
 const double /*s*/, // entropy
 const double /*h*/) // Atmosphere height
{
  vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    const double z = tess.GetCellCM(static_cast<int>(i)).y;
    //    const double x = tess.GetCellCM(static_cast<int>(i)).x;
    res[i].density = 1;
    res[i].pressure = 2-a*z;
    res[i].velocity = Vector2D(0,0);
  }
  return res;
}
