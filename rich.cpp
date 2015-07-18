#include <iostream>
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  SimData sim_data(5./3.,
		   Vector2D(0,-10),
		   1,
		   1.1);
  hdsim& sim = sim_data.getSim();
  my_main_loop(sim);
  return 0;
}
