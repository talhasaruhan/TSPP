#include "Approx.h"
#include "FlagArray.h"
#include "MST.h"
#include <iostream>
#include <stack>
#include <vector>

namespace mst_approx
{
using namespace data;

Path Approximate(const int num_verts, const Matrix& distance_matrix, const int start_vert)
{
	MSTParentArray mst_parents = MSTPrim(num_verts, distance_matrix, start_vert);
	MSTEdgeArr edge_arr = BuildMSTEdgeArr(num_verts, mst_parents);
	const int num_edges = 2 * num_verts - 2;
	assert(edge_arr.size() == num_edges);

	FlagArr is_visited(num_verts);
	std::stack<int> stack;
	stack.push(start_vert);
	Path tspath;
	tspath.reserve(num_verts);

	while (!stack.empty())
	{
		const int u = stack.top();
		stack.pop();
		is_visited.Set(u);
		tspath.push_back(u);

		int adj = GetFirstAdjVert(edge_arr, u);
		for (int e = adj; e < num_edges; ++e)
		{
			MSTEdge& edge = edge_arr[e];
			if (edge.u != u)
			{
				break;
			}
			int v = edge.v;
			if (!is_visited.Query(v))
			{
				stack.push(v);
			}
		}
	}

	return tspath;
}
}; // namespace mst_approx

namespace iterative_improv
{
void improv_2opt(const int num_verts,
                 const Matrix& distance_matrix,
                 const OptStrategy strategy,
                 Path* path)
{
	const int num_edges = num_verts - 1;

	bool improv_found = true;
	float total_gain = 0.0f;

	printf("2-OPT w/ MST 2-approx as input\n");
	float initial_length = geometry::CalculatePathLength(*path, distance_matrix);
	printf("Initial length: %f\n", initial_length);

	float max_gain = 0.0f;
	struct Opt2Segment
	{
		int i;
		int j;
	} max_gain_segment;
	const float epsilon = 1e-4;

	/* Iterate until no local improvement is found. */
	int step = 0;
	while (improv_found)
	{
		improv_found = false;
		max_gain = epsilon;
		for (int i = 0; i < num_edges; ++i)
		{
			const int e1u = (*path)[i];
			const int e1v = (*path)[i + 1];
			float edge1_len = distance_matrix.get(e1u, e1v);
			for (int j = i + 2; j < num_edges; ++j)
			{
				const int e2u = (*path)[j];
				const int e2v = (*path)[j + 1];
				float edge2_len = distance_matrix.get(e2u, e2v);

				float gain = edge1_len + edge2_len - distance_matrix.get(e1u, e2u) -
				             distance_matrix.get(e1v, e2v);
				/* Do the first move that has positive gain. */
				if (strategy == OptStrategy::FirstGain)
				{
					if (gain > epsilon)
					{
						improv_found = true;
						float _l0 = geometry::CalculatePathLength(*path, distance_matrix);
						std::reverse(path->begin() + i + 1, path->begin() + j + 1);
						float _l1 = geometry::CalculatePathLength(*path, distance_matrix);
						assert(fabs(_l0 - gain - _l1) < epsilon);
						total_gain += gain;
						break;
					}
				}
				else if (strategy == OptStrategy::BestGain)
				{
					if (gain > max_gain)
					{
						max_gain = gain;
						max_gain_segment = {i, j};
					}
				}
			}
			if (strategy == OptStrategy::FirstGain && improv_found)
			{
				break;
			}
		}
		if (strategy == OptStrategy::BestGain)
		{
			/* Apply the best opt2 move. */
			std::reverse(path->begin() + max_gain_segment.i + 1,
			             path->begin() + max_gain_segment.j + 1);
			total_gain += max_gain;
		}
		step++;
	}

	float _l = geometry::CalculatePathLength(*path, distance_matrix);
	printf("2OPT, Final length: %f | %f\n", initial_length - total_gain, _l);
}

void improv_3opt(const int num_verts,
                 const Matrix& distance_matrix,
                 const OptStrategy strategy,
                 Path* path)
{
	const int num_edges = num_verts - 1;

	bool improv_found = true;
	float total_gain = 0.0f;

	printf("3-OPT w/ MST 2-approx as input\n");
	float initial_length = geometry::CalculatePathLength(*path, distance_matrix);
	printf("Initial length: %f\n", initial_length);

	float max_gain = 0.0f;
	enum Opt3MoveCase
	{
		case_0, // No change
		case_1, // One possible 2-opt move, one segment left unchanged
		case_2, // One possible 2-opt move, one segment left unchanged
		case_3, // One possible 2-opt move, one segment left unchanged
		case_4  // Equivalent to three or more 2-opt moves. ad, be, cf
	};
	struct Opt3Segment
	{
		int i;
		int j;
		int k;
		Opt3MoveCase move_case;
	} max_gain_segment;
	const float epsilon = 1e-4;

	auto ApplyOpt3Move = [](Path* _path, Opt3Segment opt3move) {
		if (opt3move.move_case == Opt3MoveCase::case_0)
		{
		}
		else if (opt3move.move_case == Opt3MoveCase::case_1)
		{
			std::reverse(_path->begin() + opt3move.i + 1, _path->begin() + opt3move.j + 1);
		}
		else if (opt3move.move_case == Opt3MoveCase::case_2)
		{
			std::reverse(_path->begin() + opt3move.j + 1, _path->begin() + opt3move.k + 1);
		}
		else if (opt3move.move_case == Opt3MoveCase::case_3)
		{
			std::reverse(_path->begin() + opt3move.i + 1, _path->begin() + opt3move.k + 1);
		}
		else if (opt3move.move_case == Opt3MoveCase::case_4)
		{
			/* No need to copy the whole thing or do 3 seperate reverse operations. */
			Path _cpy(opt3move.k - opt3move.i);
			for (int x = 0, _end = opt3move.k - opt3move.j; x < _end; ++x)
			{
				_cpy[x] = (*_path)[opt3move.j + 1 + x];
			}
			for (int x = 0, _end = opt3move.j - opt3move.i; x < _end; ++x)
			{
				_cpy[opt3move.k-opt3move.j+x] = (*_path)[opt3move.i + 1 + x];
			}
			for (int x = 0, _end = opt3move.k - opt3move.i; x < _end; ++x)
			{
				(*_path)[opt3move.i + 1 + x] = _cpy[x];
			}
		}
	};

	/* Iterate until no local improvement is found. */
	int step = 0;
	while (improv_found)
	{
		improv_found = false;
		max_gain = epsilon;
		for (int i = 0; i < num_edges; ++i)
		{
			const int A = (*path)[i];
			const int B = (*path)[i + 1];
			float dAB = distance_matrix.get(A, B);
			for (int j = i + 2; j < num_edges; ++j)
			{
				const int C = (*path)[j];
				const int D = (*path)[j + 1];
				float dCD = distance_matrix.get(C, D);
				float dAC = distance_matrix.get(A, C);
				float dAD = distance_matrix.get(A, D);
				float dBD = distance_matrix.get(B, D);
				for (int k = j + 2; k < num_edges; ++k)
				{
					const int E = (*path)[k];
					const int F = (*path)[k + 1];
					float dEF = distance_matrix.get(E, F);
					float dEB = distance_matrix.get(E, B);
					float dFB = distance_matrix.get(F, B);
					float dCE = distance_matrix.get(C, E);
					float dDF = distance_matrix.get(D, F);
					float dEA = distance_matrix.get(E, A);
					float dCF = distance_matrix.get(C, F);

					float d0 = dAB + dCD + dEF;
					float d1 = dAC + dBD + dEF;
					float d2 = dAB + dCE + dDF;
					float d3 = dAD + dEB + dCF;
					float d4 = dFB + dCD + dEA;

					Opt3MoveCase c = Opt3MoveCase::case_0;
					float gain = 0.0f;
					if (d0 > d1)
					{
						c = Opt3MoveCase::case_1;
						gain = d0 - d1;
					}
					else if (d0 > d2)
					{
						c = Opt3MoveCase::case_2;
						gain = d0 - d2;
					}
					else if (d0 > d4)
					{
						c = Opt3MoveCase::case_3;
						gain = d0 - d4;
					}
					else if (d0 > d3)
					{
						c = Opt3MoveCase::case_4;
						gain = d0 - d3;
					}

					/* Do the first move that has positive gain. */
					if (strategy == OptStrategy::FirstGain)
					{
						if (gain > epsilon)
						{
							improv_found = true;
							float _l0 = geometry::CalculatePathLength(*path, distance_matrix);
							Path _path_cpy(path->begin(), path->end());
							ApplyOpt3Move(path, {i, j, k, c});
							float _l1 = geometry::CalculatePathLength(*path, distance_matrix);
							assert(fabs(_l0 - gain - _l1) < epsilon);
							total_gain += gain;
							break;
						}
					}
					else if (strategy == OptStrategy::BestGain)
					{
						if (gain > max_gain)
						{
							max_gain = gain;
							max_gain_segment = {i, j, k, c};
						}
					}
				}
				if (strategy == OptStrategy::FirstGain && improv_found)
				{
					break;
				}
			}
			if (strategy == OptStrategy::FirstGain && improv_found)
			{
				break;
			}
		}
		if (strategy == OptStrategy::BestGain)
		{
			/* Apply the best opt3 move. */
			ApplyOpt3Move(path, max_gain_segment);
			total_gain += max_gain;
		}
		step++;
	}

	float _l = geometry::CalculatePathLength(*path, distance_matrix);
	printf("3OPT, Final length: %f | %f\n", initial_length - total_gain, _l);
}
}; // namespace iterative_improv
