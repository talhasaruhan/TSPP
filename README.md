# TSPP
Travelling Salesman Path Problem Solver

Currently implemented:

* MST 2-approximation
* Iterative 2-opt, 3-opt methods, (both first gain and best gain)
* Branch and Bound algorithm (using 2-Approx + 3-OPT as initial solution)
  * FixedStartVert and ArbitraryStartVert options
  
Future:
* In the Branch&Bound code, nodes used to only keep the "diff"s for constraints, but due to a bug and the time constraints, I had to roll that back to a naive way. Fix that.
* Since I was focusing on being able to find exact solutions for small-medium Ns (~100), B&B dominated the complexity. So the approximation code isn't really optimized for large N values.
