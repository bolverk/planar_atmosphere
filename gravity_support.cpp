#include "gravity_support.hpp"

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

  class Bracket
  {
  public:

    Bracket(const double low, const double high):
      low_high_(pair<double,double>(low,high)) {}

    bool operator()(int arg) const
    {
      return ((arg>=low_high_.first)&&
	      (arg<=low_high_.second));
    }

  private:
    const pair<int,int> low_high_;
  };

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
	 (static_cast<size_t>(edge.neighbors.second),
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
}

GravitySupport::GravitySupport(const Tessellation& tess,
			       const PhysicalGeometry& pg,
			       const RiemannSolver& rs):
  bottom_area_(calc_bottom_area(tess,pg)),
  rs_(rs) {}

namespace {
  double calc_total_downward_momentum(const vector<Extensive>& extensives)
  {
    const Vector2D dir(0,-1);
    double res = 0;
    for(size_t i=0;i<extensives.size();++i)
      res += ScalarProd(extensives[i].momentum,dir);
    return res;
  }
}

vector<Extensive> GravitySupport::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const EquationOfState& eos,
   const double /*time*/,
   const double dt) const
{
  const double support =
    fmax(0,calc_total_downward_momentum(extensives)/bottom_area_/dt);
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
  Conserved reflect_riemann(const RiemannSolver& rs,
			    const Tessellation& tess,
			    const Edge& edge,
			    const ComputationalCell& cc,
			    const EquationOfState& eos,
			    bool left_real)
  {
    const Primitive cell = convert_to_primitive(cc,eos);
    const Vector2D& p = Parallel(edge);
    return left_real ?
      rotate_solve_rotate_back
      (rs,
       cell,
       reflect(cell,p),
       0,
       remove_parallel_component
       (edge.vertices.second -
	tess.GetMeshPoint(edge.neighbors.first),p),
       p) :
      rotate_solve_rotate_back
      (rs,
       reflect(cell,p),
       cell,
       0,
       remove_parallel_component
       (tess.GetMeshPoint(edge.neighbors.first)-
	edge.vertices.second,p),
       p);
  }

  Conserved support_riemann(double support,
			    bool left_real)
  {
    return left_real ?
      Conserved(0,Vector2D(0,-support),0) :
      Conserved(0,Vector2D(0,support),0);
  }

  Conserved regular_riemann
  (const RiemannSolver& rs,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const Edge& edge)
  {
    const size_t left_index =
      static_cast<size_t>(edge.neighbors.first);
    const size_t right_index =
      static_cast<size_t>(edge.neighbors.second);
    const Primitive left =
      convert_to_primitive(cells[left_index], eos);
    const Primitive right =
      convert_to_primitive(cells[right_index], eos);
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
    const Vector2D r = tess.GetMeshPoint(edge.neighbors.second);
    if(r.y>edge.vertices.first.y && r.y>edge.vertices.second.y)
      return support_riemann(support,false);
    return reflect_riemann
      (rs_,
       tess,
       edge,
       cells.at(static_cast<size_t>(edge.neighbors.second)),
       eos,
       false);
  }
  if(!b(edge.neighbors.second)){
    const Vector2D r = tess.GetMeshPoint(edge.neighbors.first);
    if(r.y>edge.vertices.first.y && r.y>edge.vertices.second.y)
      return support_riemann(support,true);
    return reflect_riemann
      (rs_,
       tess,
       edge,
       cells.at(static_cast<size_t>(edge.neighbors.first)),
       eos,
       true);
  }
  return regular_riemann(rs_,
			 tess,
			 point_velocities,
			 cells,
			 eos,
			 edge);
}
