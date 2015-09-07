#include "sim_data.hpp"
#include "source/misc/mesh_generator.hpp"
#include "calc_init_cond.hpp"

SimData::SimData(double adiabatic_index,
		 const Vector2D& acceleration,
		 double entropy,
		 double atmosphere_height):
  pg_(),
  outer_(Vector2D(0,0),
	 Vector2D(1,0.2)),
  tess_(cartesian_mesh(200*2,20*2,
		       outer_.getBoundary().first,
		       outer_.getBoundary().second),
	outer_),
  eos_(adiabatic_index),
  rs_(),
  pm_(),
  grav_acc_(acceleration),
  force_(grav_acc_),
  tsf_(0.3),
  fc_(tess_,
      pg_,
      acceleration,
      rs_),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_,
		      adiabatic_index,
		      abs(acceleration),
		      entropy,
		      atmosphere_height),
       eos_,
       pm_,
       force_,
       tsf_,
       fc_,
       eu_,
       cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}
