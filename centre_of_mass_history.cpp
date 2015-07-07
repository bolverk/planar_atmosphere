#include <fstream>
#include "centre_of_mass_history.hpp"

namespace {
  Vector2D calc_cm(const hdsim& sim)
  {
    double mass = 0;
    Vector2D res = Vector2D(0,0);
    const Tessellation& tess = sim.getTessellation();
    const vector<ComputationalCell>& cells = 
      sim.getAllCells();
    const vector<Extensive>& extensives = 
      sim.getAllExtensives();
    for(size_t i=0;i<cells.size();++i){
      mass += extensives[i].mass;
      res += extensives[i].mass*
	tess.GetCellCM(static_cast<int>(i));
    }
    return res/mass;
  }
}

CentreOfMassHistory::CentreOfMassHistory(const string& fname):
  fname_(fname), time_cm_list_() {}

void CentreOfMassHistory::operator()(const hdsim& sim)
{

  time_cm_list_.push_back
    (pair<double,Vector2D>
     (sim.getTime(),calc_cm(sim)));
}

CentreOfMassHistory::~CentreOfMassHistory(void)
{
  std::ofstream f(fname_.c_str());
  for(size_t i=0;i<time_cm_list_.size();++i)
    f << time_cm_list_[i].first << " "
      << time_cm_list_[i].second.x << " " 
      << time_cm_list_[i].second.y << std::endl;
  f.close();
}

