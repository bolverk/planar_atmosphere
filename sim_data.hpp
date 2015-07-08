#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/source_terms/ConstantGravity.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "gravity_support.hpp"

class SimData
{
public:

  SimData(double adiabatic_index,
	  double acceleration,
	  double entropy,
	  double atmosphere_height);

  hdsim& getSim(void);

private:
  const SlabSymmetry pg_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  const Hllc rs_;
  const Eulerian pm_;
  const ConstantGravity grav_acc_;
  const ConservativeForce force_;
  const SimpleCFL tsf_;
  const GravitySupport fc_;
  const SimpleExtensiveUpdater eu_;
  const SimpleCellUpdater cu_;
  hdsim sim_;
};

#endif // SIM_DATA_HPP
