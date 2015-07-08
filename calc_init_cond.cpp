#include <cmath>
#include "calc_init_cond.hpp"

vector<ComputationalCell> calc_init_cond
(const Tessellation& tess,
 const double g, // adiabatic index
 const double a, // acceleration
 const double s, // entropy
 const double h) // Atmosphere height
{
  vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    const double z = tess.GetCellCM(static_cast<int>(i)).y;
    res[i].density = pow(a*(g-1)*s*h*(1.-z/h)/g,1./(g-1));
    res[i].pressure = s*pow(res[i].density,g);
    res[i].velocity = Vector2D(0,0);
  }
  return res;
}
