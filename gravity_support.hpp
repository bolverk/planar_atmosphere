#ifndef GRAVITY_SUPPORT_HPP
#define GRAVITY_SUPPORT_HPP 1

#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"

class GravitySupport: public FluxCalculator
{
public:

  GravitySupport(const Tessellation& tess,
		 const PhysicalGeometry& pg,
		 const Vector2D& acceleration,
		 const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double time,
   const double dt) const;

private:
  const double bottom_area_;
  const Vector2D acceleration_;
  const RiemannSolver& rs_;

  Conserved calcHydroFlux(const Tessellation& tess,
			  const vector<Vector2D>& point_velocities,
			  const vector<ComputationalCell>& cells,
			  const EquationOfState& eos,
			  const size_t i,
			  const double support) const;
};

#endif // GRAVITY_SUPPORT_HPP
