#include "sim_data.hpp"
#include "source/misc/mesh_generator.hpp"
#include "calc_init_cond.hpp"

SimData::SimData(void):
  pg_(),
  outer_(Vector2D(-0.05,0),
	 Vector2D(0.05,1)),
  tess_(cartesian_mesh(5,100,
		       outer_.getBoundary().first,
		       outer_.getBoundary().second),
	outer_),
  eos_(5./3.),
  rs_(),
  pm_(),
  grav_acc_(Vector2D(0,-10)),
  force_(grav_acc_),
  tsf_(0.03),
  fc_(tess_,
      pg_,
      rs_),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_),
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
