#include "gravity_support.hpp"
#include "bracket.hpp"

namespace {

  double calc_tracer_flux(size_t i,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const std::string& name,
			  const Conserved& hf)
  {
    const Edge& edge = tess.GetEdge(static_cast<int>(i));
    if(hf.Mass>0 && 
       edge.neighbors.first>0 &&
       edge.neighbors.first < tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.first)].tracers.find(name)->second;
    if(hf.Mass<0 &&
       edge.neighbors.second>0 &&
       edge.neighbors.second < tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.second)].tracers.find(name)->second;
    return 0;
  }

  pair<bool,pair<size_t,bool> > check_boundary_edge
  (const Edge& edge, const Bracket& b)
  {
    if(!b(edge.neighbors.first)){
      assert(b(edge.neighbors.second));
      return pair<bool,pair<size_t,bool> >
	(true,pair<size_t,bool>
	 (static_cast<size_t>(edge.neighbors.second),
	  false));
    }
    if(!b(edge.neighbors.second))
      return pair<bool,pair<size_t,bool> >
	(true,pair<size_t,bool>
	 (static_cast<size_t>(edge.neighbors.first),
	  true));
    return pair<bool,pair<size_t,bool> >
      (false,pair<size_t,bool>(0,false));
  }

  bool check_bottom_edge
  (const Edge& edge,
   const Bracket& b,
   const Tessellation& tess)
  {
    const pair<bool,pair<size_t,bool> > temp =
      check_boundary_edge(edge,b);
    if(!temp.first)
      return false;
    const Vector2D mp =
      tess.GetMeshPoint(static_cast<int>(temp.second.first));
    return (mp.y>edge.vertices.first.y &&
	    mp.y>edge.vertices.second.y);
  }

  double calc_bottom_area(const Tessellation& tess,
			  const PhysicalGeometry& pg)
  {
    const vector<Edge>& edges = tess.getAllEdges();
    const Bracket b(0,tess.GetPointNo());
    double res = 0;
    for(size_t i=0;i<edges.size();++i){
      const Edge& edge = edges[i];
      if(check_bottom_edge(edge,b,tess))
	res += pg.calcArea(edge);
    }
    return res;
  }

  Primitive boost(const Primitive& origin,
		  const Vector2D& v)
  {
    Primitive res = origin;
    res.Velocity += v;
    return res;
  }
}

GravitySupport::GravitySupport(const Tessellation& tess,
			       const PhysicalGeometry& pg,
			       const Vector2D& acceleration,
			       const RiemannSolver& rs):
  bottom_area_(calc_bottom_area(tess,pg)),
  acceleration_(acceleration),
  rs_(rs) {}

vector<Extensive> GravitySupport::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*extensives*/,
   const CacheData& /*cd*/,
   const EquationOfState& eos,
   const double /*time*/,
   const double dt) const
{
  const double support = dt*abs(acceleration_);
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<tess.getAllEdges().size();++i){
    const Conserved hydro_flux = calcHydroFlux
      (tess, point_velocities, cells, eos, i, support);
    res[i].mass = hydro_flux.Mass;
    res[i].momentum = hydro_flux.Momentum;
    res[i].energy = hydro_flux.Energy;
    for(std::map<std::string,double>::const_iterator it = 
	  cells.front().tracers.begin();
	it!=cells.front().tracers.end();++it)
      res[i].tracers[it->first] =
	calc_tracer_flux(i,tess,cells,it->first,hydro_flux);
  }
  return res;
}

namespace {
  Conserved support_riemann(const RiemannSolver& rs,
			    const Tessellation& tess,
			    const Edge& edge,
			    const ComputationalCell& cc,
			    const EquationOfState& eos,
			    double support,
			    bool left_real)
  {
    const Primitive cell = convert_to_primitive(cc,eos);
    const Vector2D& p = Parallel(edge);
    Conserved res = left_real ?
      rotate_solve_rotate_back
      (rs,
       cell,
       boost(reflect(cell,p),Vector2D(0,support)),
       0,
       remove_parallel_component
       (edge.vertices.second -
	tess.GetMeshPoint(edge.neighbors.first),p),
       p) :
      rotate_solve_rotate_back
      (rs,
       boost(reflect(cell,p),Vector2D(0,support)),
       cell,
       0,
       remove_parallel_component
       (tess.GetMeshPoint(edge.neighbors.second)-
	edge.vertices.second,p),
       p);
    res.Mass = 0;
    res.Energy = 0;
    return res;
  }

  ComputationalCell gravinterpolate
  (const ComputationalCell& source,
   const Vector2D& cm,
   const Vector2D& centroid,
   const Vector2D& gravity)
  {
    ComputationalCell res = source;
    res.pressure += res.density*ScalarProd(gravity,centroid-cm);
    return res;
  }

  Conserved bulk_riemann
  (const RiemannSolver& rs,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const Edge& edge,
   const Vector2D& gravity)
  {
    const Vector2D pos_left =
      tess.GetCellCM(edge.neighbors.first);
    const Vector2D pos_right =
      tess.GetCellCM(edge.neighbors.second);
    const size_t left_index =
      static_cast<size_t>(edge.neighbors.first);
    const size_t right_index =
      static_cast<size_t>(edge.neighbors.second);
    const Vector2D centroid =
      0.5*(edge.vertices.first+edge.vertices.second);
    const Primitive left =
      convert_to_primitive
      (gravinterpolate
       (cells[left_index],
	pos_left,
	centroid,
	gravity),
       eos);
    const Primitive right =
      convert_to_primitive
      (gravinterpolate
       (cells[right_index],
	pos_right,
	centroid,
	gravity),
       eos);
    const Vector2D p = Parallel(edge);
    const Vector2D n =
      tess.GetMeshPoint(edge.neighbors.second) -
      tess.GetMeshPoint(edge.neighbors.first);
    const double velocity = Projection
      (tess.CalcFaceVelocity
       (point_velocities.at(left_index),
	point_velocities.at(right_index),
	tess.GetCellCM(edge.neighbors.first),
	tess.GetCellCM(edge.neighbors.second),
	calc_centroid(edge)),n);
    return rotate_solve_rotate_back
      (rs,left,right,velocity,n,p);
  }
}

Conserved GravitySupport::calcHydroFlux
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const size_t i,
 double support) const
{
  const Edge& edge = tess.GetEdge(static_cast<int>(i));
  const Bracket b(0,tess.GetPointNo());
  if(!b(edge.neighbors.first)){
    assert(b(edge.neighbors.second));
    return support_riemann
      (rs_,
       tess,
       edge,
       cells.at(static_cast<size_t>(edge.neighbors.second)),
       eos,
       support,
       false);
  }
  if(!b(edge.neighbors.second))
    return support_riemann
      (rs_,
       tess,
       edge,
       cells.at(static_cast<size_t>(edge.neighbors.first)),
       eos,
       support,
       true);
  return bulk_riemann(rs_,
		      tess,
		      point_velocities,
		      cells,
		      eos,
		      edge,
		      acceleration_);
}
