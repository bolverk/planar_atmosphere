#ifndef CENTRE_OF_MASS_HISTORY_HPP
#define CENTRE_OF_MASS_HISTORY_HPP 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"

class CentreOfMassHistory: public DiagnosticFunction
{
public:

  CentreOfMassHistory(const string& fname);

  void operator()(const hdsim& sim);

  ~CentreOfMassHistory(void);

private:
  const string fname_;
  mutable vector<pair<double,Vector2D> > time_cm_list_;
};

#endif // CENTRE_OF_MASS_HISTORY_HPP
