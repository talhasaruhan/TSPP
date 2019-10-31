#include "BranchAndBound.h"
#include "Approx.h"
#include <cassert>
#include <chrono>
#include <limits>

namespace branch_and_bound
{
inline void ApplyConstraint(const int n, Constraints* constraints, const ConstraintDiff& diff)
{
	constraints->Set(EdgeId(n, diff.row, diff.col), diff.constraint);
	constraints->Set(EdgeId(n, diff.col, diff.row), diff.constraint);
}
inline void RevertConstraint(const int n, Constraints* constraints, const ConstraintDiff& diff)
{
	constraints->Set(EdgeId(n, diff.row, diff.col), ConstraintTypes::Allowed);
	constraints->Set(EdgeId(n, diff.col, diff.row), ConstraintTypes::Allowed);
}
VertexPair CandidateEdge(const int n, const Constraints& constraints, const Matrix& distance_matrix)
{
	assert(n == distance_matrix.num_rows && n == distance_matrix.num_cols &&
	       n * n == constraints.size);
	int min_edge = -1;
	float min_edge_weight = std::numeric_limits<float>::max();
	for (int row = 0; row < n; ++row)
	{
		for (int col = row + 1; col < n; ++col)
		{
			edge_t edge_id = EdgeId(n, row, col);
			if (constraints.Query(edge_id) == ConstraintTypes::Allowed)
			{
				float weight = distance_matrix[row][col];
				if (weight < min_edge_weight)
				{
					min_edge = edge_id;
					min_edge_weight = weight;
				}
			}
		}
	}
	return EdgeId(n, min_edge);
}

DBGExpConstraint DBGExpandConstraints(const int n, const Constraints& constraints)
{
	DBGExpConstraint exp;
	for (int row = 0; row < n; ++row)
	{
		for (int col = 0; col < n; ++col)
		{
			exp.push_back((ConstraintTypes)constraints.Query(EdgeId(n, row, col)));
		}
	}
	return exp;
}

/* Tighten constraints to ensure the validity of the partial solution.
 * Returns whether given constraints describe a leaf node or not. */
template <TSPPType tspp_type>
TightenConstraintsReturnCode TightenConstraints(const int n,
                                                const Constraints& constraints,
                                                ConstraintTighteningPatch* patch,
                                                const int start_vert)
{
	/* For a vertex:
	 * If two neighboring edges are required all other edges must be disallowed.
	 * If all but 2 edges are disallowed those 2 edge must be required.
	 * We check for the cycle condition during the search,
	 * so we don't need to impose that constraint here. */
	bool is_leaf = true;
	struct Stats
	{
		int num_allowed = 0;
		int num_disallowed = 0;
		int num_required = 0;
	};
	std::vector<Stats> stats(n);
	int num_endpoints = 0;
	for (int row = 0; row < n; ++row)
	{
		Stats& stat = stats[row];
		for (int col = 0; col < n; ++col)
		{
			switch (constraints.Query(EdgeId(n, row, col)))
			{
			case ConstraintTypes::Allowed:
				stat.num_allowed++;
				break;
			case ConstraintTypes::Disallowed:
				stat.num_disallowed++;
				break;
			case ConstraintTypes::Required:
				stat.num_required++;
				break;
			default:
				assert(false);
				break;
			}
		}
		if (stat.num_required > 2 || stat.num_disallowed > n - 1)
		{
			return TightenConstraintsReturnCode::NonValid;
		}
		if (stat.num_required == 1 && stat.num_allowed == 0)
		{
			num_endpoints++;
		}
		if (row == start_vert && stat.num_required > 1)
		{
			return TightenConstraintsReturnCode::NonValid;
		}
	}

	if (num_endpoints > 2)
	{
		return TightenConstraintsReturnCode::NonValid;
	}

	for (int row = 0; row < n; ++row)
	{
		const Stats& stat = stats[row];

		bool tighten = true;
		ConstraintTypes set;
		if (!stat.num_allowed)
		{
			tighten = false;
		}
		else if (stat.num_required == 2)
		{
			set = ConstraintTypes::Disallowed;
		}
		else if (stat.num_disallowed == n - 1)
		{
			set = ConstraintTypes::Required;
		}
		else
		{
			tighten = false;
			is_leaf = false;
		}

		if (tighten)
		{
			for (int col = 0; col < n; ++col)
			{
				ConstraintTypes constraint =
				    (ConstraintTypes)constraints.Query(EdgeId(n, row, col));
				if (constraint == ConstraintTypes::Allowed)
				{
					patch->emplace_back(ConstraintDiff{(int8_t)row, (int8_t)col, set});
				}
			}
		}
	}

	std::sort(patch->begin(), patch->end());
	patch->erase(std::unique(patch->begin(), patch->end()), patch->end());
	return (TightenConstraintsReturnCode)is_leaf;
}

template <TSPPType tspp_type>
TightenConstraintsReturnCode IterativelyTightenConstraints(const int n,
                                                           Constraints* constraints,
                                                           ConstraintTighteningPatch* patch,
                                                           const int start_vert)
{
	TightenConstraintsReturnCode ret;
	int prev_size = 0;
	while (true)
	{
		ret = TightenConstraints<tspp_type>(n, *constraints, patch, start_vert);
		ApplyConstraintTighteningPatch(n, constraints, *patch);
		if (patch->size() == prev_size)
		{
			break;
		}
		prev_size = patch->size();
	}
	return ret;
}

void ApplyConstraintTighteningPatch(const int n,
                                    Constraints* constraints,
                                    const ConstraintTighteningPatch& patch)
{
	Constraints& _constraints = *constraints;
	for (const ConstraintDiff& diff : patch)
	{
		ApplyConstraint(n, constraints, diff);
	}
}
void RevertConstraintTighteningPatch(const int n,
                                     Constraints* constraints,
                                     const ConstraintTighteningPatch& patch)
{
	Constraints& _constraints = *constraints;
	for (const ConstraintDiff& diff : patch)
	{
		RevertConstraint(n, constraints, diff);
	}
}

float CalculateSolutionLength(const int n,
                              const Constraints& constraints,
                              const Matrix& distance_matrix)
{
	float len = 0.0f;
	for (int row = 0; row < n; ++row)
	{
		for (int col = row + 1; col < n; ++col)
		{
			if (constraints.Query(EdgeId(n, row, col)) == ConstraintTypes::Required)
			{
				len += distance_matrix[row][col];
			}
		}
	}
	return len;
}

/* Calculate simple lower bound for the given edge set
 * by summing two minimum weighted adjacent edges (excluding
 * disallowed edges and prioritizing required edges)
 * for every vertex in the graph. */
template <TSPPType tspp_type>
float LowerBound(const Constraints& constraints,
                 const Matrix& distance_matrix,
                 const int start_vert)
{
	const int n = distance_matrix.num_rows;
	/* Must be a square matrix. */
	assert(n == distance_matrix.num_cols);

	float cost = 0.0f;
	float max_min_edge_weight = 0.0f;
	float start_vert_second_min;
	for (int row = 0; row < n; ++row)
	{
		int num_required = 0;
		/* Allocate 3 edge objects on the stack, at the end,
		 * if num_required == 2: use first 2,
		 * if num_required == 1: use first 2,
		 * if num_required == 0: use last 2 */
		float edge_weights[3];
		float& _req_edge_weight = edge_weights[0];
		float& min_edge_weight = edge_weights[1];
		float& second_min_edge_weight = edge_weights[2];
		min_edge_weight = std::numeric_limits<float>::max();
		second_min_edge_weight = std::numeric_limits<float>::max();

		for (int col = 0; (col < n) && (num_required < 2); ++col)
		{
			const float len = distance_matrix.get(row, col);
			edge_t edge_id = EdgeId(n, row, col);
			const ConstraintTypes& edge_constraint = (ConstraintTypes)constraints.Query(edge_id);
			switch (edge_constraint)
			{
			/* Prioritize required edges. */
			case ConstraintTypes::Required:
				edge_weights[num_required++] = len;
				break;
			case ConstraintTypes::Allowed:
				if (len < min_edge_weight)
				{
					second_min_edge_weight = min_edge_weight;
					min_edge_weight = len;
				}
				else if (len < second_min_edge_weight)
				{
					second_min_edge_weight = len;
				}
			/* Exlude disallowed edges. */
			case ConstraintTypes::Disallowed:
				break;
			default:
				assert(false);
				break;
			}
		}

		/* Calculate the sum of the selected two adjacent edges. */
		if (num_required > 0)
		{
			cost += edge_weights[0] + edge_weights[1];
			if constexpr (tspp_type == TSPPType::ArbitraryStartVert)
			{
				max_min_edge_weight =
				    std::max(max_min_edge_weight, std::max(edge_weights[0], edge_weights[1]));
			}
			else if constexpr (tspp_type == TSPPType::FixedStartVert)
			{
				if (row == start_vert)
				{
					start_vert_second_min = std::max(edge_weights[0], edge_weights[1]);
				}
			}
		}
		else
		{
			cost += edge_weights[1] + edge_weights[2];
			if constexpr (tspp_type == TSPPType::ArbitraryStartVert)
			{
				max_min_edge_weight = std::max(max_min_edge_weight, edge_weights[2]);
			}
			else if constexpr (tspp_type == TSPPType::FixedStartVert)
			{
				if (row == start_vert)
				{
					start_vert_second_min = std::max(edge_weights[1], edge_weights[2]);
				}
			}
		}
	}

	if constexpr (tspp_type == TSPPType::ArbitraryStartVert)
	{
		return cost * 0.5 - max_min_edge_weight;
	}
	else if constexpr (tspp_type == TSPPType::FixedStartVert)
	{
		return cost * 0.5 - start_vert_second_min;
	}
}

/* Do a bidirectional search from u to v using only allowed nodes.
 * If a valid path is found, it implies that insertion of edge (u, v)
 * will create a subtour. */
bool YieldsCycle(const int n, const Constraints& constraints, int u, int v)
{
	typedef Bitset<1> Flags;
	typedef std::queue<int> Queue;

	Flags visited_u(n);
	Flags visited_v(n);
	Queue queue_u;
	Queue queue_v;

	queue_u.push(u);
	queue_v.push(v);

	auto Step = [](const int n, const Constraints& constraints, Queue& queue, Flags& visited) {
		int vert = queue.front();
		queue.pop();
		visited.Set(vert);
		assert(visited.Query(vert));
		/* Push unvisited neighbors into the queue only if the edge is required. */
		for (int i = 0; i < n; ++i)
		{
			if (i == vert || visited.Query(i) ||
			    constraints.Query(EdgeId(n, vert, i)) != ConstraintTypes::Required)
			{
				continue;
			}
			queue.push(i);
		}
	};

	while (true)
	{
		bool qu_empty = queue_u.empty();
		bool qv_empty = queue_v.empty();
		if (qu_empty && qv_empty)
		{
			break;
		}
		if (!qu_empty)
		{
			Step(n, constraints, queue_u, visited_u);
		}
		if (!qv_empty)
		{
			Step(n, constraints, queue_v, visited_v);
		}
		bool intersect = BitsetIntersect(n, visited_u, visited_v);
		if (intersect)
		{
			return true;
		}
	}

	return false;
}

void PrintConstraints(const int n, const Constraints& constraints)
{
	printf("XX |");
	for (int col = 0; col < n; ++col)
	{
		printf("%02d |", col);
	}
	printf("\n");
	for (int row = 0; row < n; ++row)
	{
		printf("%02d |", row);
		for (int col = 0; col < n; ++col)
		{
			ConstraintTypes type = (ConstraintTypes)constraints.Query(EdgeId(n, row, col));
			switch (type)
			{
			case ConstraintTypes::Allowed:
				printf(" - |");
				break;
			case ConstraintTypes::Disallowed:
				printf(" D |");
				break;
			case ConstraintTypes::Required:
				printf(" R |");
				break;
			}
		}
		printf("\n");
	}
	printf("\n");
}

template <TSPPType tspp_type>
Path SolveMetricTSPP(const int n,
                     const Matrix& distance_matrix,
                     const int start_vert,
                     const int exec_time_limit)
{
	const auto exec_start_time = std::chrono::high_resolution_clock::now();
	assert(n == distance_matrix.num_rows && n == distance_matrix.num_cols);

	/* Create an initial approximation using MST method (2-approximator). */
	Path initial_approx = mst_approx::Approximate(n, distance_matrix, start_vert);
	iterative_improv::improv_3opt(n,
	                              distance_matrix,
	                              iterative_improv::OptStrategy::FirstGain,
	                              &initial_approx);

	/* Stacks hold the nodes to be visited and the path from root to node. */
	std::stack<Node> tobevisited;
	struct BacktrackInfo
	{
		Node parent;
		ConstraintTighteningPatch constraint_patch;
	};
	std::deque<BacktrackInfo> backtrackinfo;
	/* Current set of constraints. */
	// Constraints cur_constraints(n * n, ConstraintTypes::Allowed);
	// for (int i = 0; i < n; ++i)
	//{
	//	cur_constraints.Set(EdgeId(n, i, i), ConstraintTypes::Disallowed);
	//}
	float cur_best_length;
	Constraints cur_best;

	/* Initialize the root node: empty set, all edges are allowed. */
	tobevisited.emplace(
	    Node{0, 0, ConstraintTypes::Disallowed, 0, Constraints(n * n, ConstraintTypes::Allowed)});
	for (int i = 0; i < n; ++i)
	{
		tobevisited.top().constraints.Set(EdgeId(n, i, i), ConstraintTypes::Disallowed);
	}
	/* Initialize the current best length. */
	bool found_better_than_initial_approx = false;
	cur_best_length = geometry::CalculatePathLength(initial_approx, distance_matrix);
	printf("Initial approx. cost: %f\n", cur_best_length);

	int step = 0;
	enum BacktrackCode
	{
		NoBacktrack,
		Leaf,
		Pruned,
		Nonvalid
	};
	//int backtrack = 0;
	while (!tobevisited.empty())
	{
		const auto t = std::chrono::high_resolution_clock::now();
		const int elapsed_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(t - exec_start_time).count();
		if (elapsed_ms > exec_time_limit)
		{
			printf("Exec. time limit is surpassed.\n%d steps.\nBest result found so far: %f\n",
			       step,
			       cur_best_length);
			if (found_better_than_initial_approx)
			{
				printf("This result was achieved with branch and bound.\n");
			}
			else
			{
				printf("No better results than initial approximation has been found.\n");
			}
			if (found_better_than_initial_approx)
			{
				AssertConstraintsValidity<tspp_type>(n, cur_best, start_vert);
				Path path = ConstraintsToPath(n, cur_best, start_vert);
				return path;
			}
			else
			{
				return initial_approx;
			}
		}

		step++;
		Node node = tobevisited.top();
		tobevisited.pop();

		Constraints& cur_constraints = node.constraints;

		if (step > 1)
		{
			ApplyConstraint(n, &cur_constraints, {node.row, node.col, node.diff_constraint});
		}

		ConstraintTighteningPatch constraint_tightening_patch;
		TightenConstraintsReturnCode tightenretcode =
		    IterativelyTightenConstraints<tspp_type>(n,
		                                             &cur_constraints,
		                                             &constraint_tightening_patch,
		                                             start_vert);

		//backtrack = BacktrackCode::NoBacktrack;

		if (tightenretcode == TightenConstraintsReturnCode::NonValid)
		{
			/* Discard the node. */
			//backtrack = BacktrackCode::Nonvalid;
		}
		else if (tightenretcode == TightenConstraintsReturnCode::NonValid)
		{
			/* Discard the node. */
			//backtrack = BacktrackCode::Nonvalid;
		}
		else if (tightenretcode == TightenConstraintsReturnCode::ValidAndLeaf)
		{
			/* Calculate the length of the solution node corresponds to. */
			float length = CalculateSolutionLength(n, cur_constraints, distance_matrix);
			if (length < cur_best_length)
			{
				found_better_than_initial_approx = true;
				cur_best_length = length;
				cur_best = node.constraints;
			}
			//backtrack = BacktrackCode::Leaf;
		}
		else
		{
			float node_lower_bound =
			    LowerBound<tspp_type>(cur_constraints, distance_matrix, start_vert);

			if (node_lower_bound > cur_best_length)
			{
				/* Discard the node. */
				//backtrack = BacktrackCode::Pruned;
			}
			else
			{
				backtrackinfo.push_back(
				    {{node.row, node.col, node.diff_constraint, 2}, constraint_tightening_patch});
				/* Find the candidate edge (shortest allowed edge) e */
				int candidate_row;
				int candidate_col;
				/* Cycle check is done lazily. */
				while (true)
				{
					const auto candidate = CandidateEdge(n, cur_constraints, distance_matrix);
					candidate_row = candidate.row;
					candidate_col = candidate.col;
					if (EdgeId(n, candidate_row, candidate_col) < 0)
					{
						/* No viable edge exists. */
						candidate_row = -1;
						candidate_col = -1;
						break;
					}
					if (YieldsCycle(n, cur_constraints, candidate_row, candidate_col))
					{
						ApplyConstraint(n,
						                &cur_constraints,
						                {(int8_t)candidate_row,
						                 (int8_t)candidate_col,
						                 ConstraintTypes::Disallowed});
					}
					else
					{
						break;
					}
				}
				if (candidate_row >= 0 && candidate_col >= 0)
				{
					/* Insert left children (requires e) */
					tobevisited.emplace(Node{(int8_t)candidate_row,
					                         (int8_t)candidate_col,
					                         ConstraintTypes::Disallowed,
					                         0,
					                         cur_constraints});
					/* Insert left children (disallows e) */
					tobevisited.emplace(Node{(int8_t)candidate_row,
					                         (int8_t)candidate_col,
					                         ConstraintTypes::Required,
					                         0,
					                         cur_constraints});
					//backtrack = BacktrackCode::NoBacktrack;
				}
				else
				{
					//backtrack = BacktrackCode::Nonvalid;
				}
			}
		}

		/* Backtrack to parent until parent has a viable children. */
		// if (backtrack)
		//{
		//	// RevertConstraint(n, &cur_constraints, {node.row, node.col});
		//	// RevertConstraintTighteningPatch(n, &cur_constraints, constraint_tightening_patch);
		//	BacktrackInfo* bck = &backtrackinfo.back();
		//	int backtrack_steps = 1;
		//	while (bck->parent.counter == 1)
		//	{
		//		backtrack_steps++;
		//		/* node was the last viable child of the parent, backtrack. */
		//		backtrackinfo.pop_back();
		//		// RevertConstraint(n, &cur_constraints, {bck->parent.row, bck->parent.col});
		//		// RevertConstraintTighteningPatch(n, &cur_constraints, bck->constraint_patch);
		//		if (backtrackinfo.size())
		//		{
		//			bck = &backtrackinfo.back();
		//		}
		//		else
		//		{
		//			break;
		//		}
		//	}
		//	bck->parent.counter--;
		//}
	}
	printf("steps: %d, min length: %f\n", step, cur_best_length);
	AssertConstraintsValidity<tspp_type>(n, cur_best, start_vert);
	Path path = ConstraintsToPath(n, cur_best, start_vert);
	return path;
}

/* Check the validity of the constraint matrix, (only in debug mode) */
template <TSPPType tspp_type>
void AssertConstraintsValidity(const int n, const Constraints& constraints, const int start_vert)
{
#ifndef NDEBUG
	/* Vertices with only 1 edges as neighbors. */
	int n_endpoints = 0;
	for (int row = 0; row < n; ++row)
	{
		int num_required = 0;
		int num_disallowed = 0;
		for (int col = 0; col < n; ++col)
		{
			edge_t edge_id = EdgeId(n, row, col);
			ConstraintTypes constraint = (ConstraintTypes)constraints.Query(edge_id);
			assert(constraint == ConstraintTypes::Required ||
			       constraint == ConstraintTypes::Disallowed);
			if (constraint == ConstraintTypes::Required)
			{
				num_required++;
			}
			else if (constraint == ConstraintTypes::Disallowed)
			{
				num_disallowed++;
			}
		}
		if constexpr (tspp_type == TSPPType::FixedStartVert)
		{
			assert(row != start_vert || num_required == 1);
		}
		assert(num_required >= 0);
		assert(num_required <= 2);
		if (num_required == 1)
		{
			n_endpoints++;
		}
		assert(n_endpoints <= 2);
	}
	assert(n_endpoints == 2);
	printf("Asserted correctness of the solution.\n");
#else
	printf("Not asserted correctness of the solution. Build in DEBUG mode.\n");
#endif
}

Path ConstraintsToPath(const int n, const Constraints& constraints, const int start_vert)
{
	int prev_vert = -1;
	int vert = start_vert;
	Path path;
	path.reserve(n);
	while (true)
	{
		/* Aside from start and end, every vertex must have 2 adj edges. */
		bool endpoint = true;
		path.push_back(vert);
		for (int other = 0; other < n; ++other)
		{
			if (other == prev_vert)
			{
				continue;
			}
			ConstraintTypes constraint = (ConstraintTypes)constraints.Query(EdgeId(n, vert, other));
			if (constraint == ConstraintTypes::Required)
			{
				prev_vert = vert;
				vert = other;
				endpoint = false;
				break;
			}
		}
		if (vert != start_vert && endpoint)
		{
			break;
		}
	}
	#ifndef NDEBUG
	/* If we're in debug mode, make sure the path has every vertex exactly once. */
	Path _path_cpy(path.begin(), path.end());
	std::sort(_path_cpy.begin(), _path_cpy.end());
	assert(std::unique(_path_cpy.begin(), _path_cpy.end()) == _path_cpy.end());
	#endif
	return path;
}
}; // namespace branch_and_bound
