#include "my_main_loop.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "write_cycle.hpp"
#include "centre_of_mass_history.hpp"

using namespace simulation2d;

void my_main_loop(hdsim& sim)
{
  write_snapshot_to_hdf5(sim,"initial.h5");
  const double tf = 10;
  SafeTimeTermination term_cond(tf,1e6);
  MultipleDiagnostics diag
    (VectorInitialiser<DiagnosticFunction*>
     (new ConsecutiveSnapshots(new ConstantTimeInterval(tf/1000),
			       new Rubric("snapshot_",".h5")))
     (new WriteTime("time.txt"))
     (new WriteCycle("cycle.txt"))
     (new CentreOfMassHistory("cm_history.txt"))
     ());
  main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance,
	    &diag);
  write_snapshot_to_hdf5(sim,"final.h5");
}
