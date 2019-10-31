#pragma once
#include "FlagArray.h"
#include "MST.h"
#include "Matrix.h"
#include <algorithm>
#include <queue>
#include <stack>
#include <vector>

namespace branch_and_bound
{
/* Search tree starts with an empty set, and at each step, each node
 * branches into two nodes, one that includes 'e' and one that doesn't.
 * 'e' is selected in the increasing length order. Node is disacarded
 * if the lower bound for the given node is greater than the current
 * best solution or approximation found. */
enum ConstraintTypes : uint8_t
{
	/* Let the default initialized value be Allowed. */
	Allowed,
	Disallowed,
	Required,
	_NUM_CONSTRAINTS
};
/* Reserve 2 bits for each constraint. */
static_assert(_NUM_CONSTRAINTS < 4);

typedef int edge_t;
inline edge_t EdgeId(int n, int row, int col)
{
	return row * n + col;
}
struct VertexPair
{
	int row;
	int col;
};
inline VertexPair EdgeId(int n, edge_t edge_id)
{
	std::div_t res = std::div(edge_id, n);
	return {res.quot, res.rem};
};

struct Edge
{
	edge_t edge_id;
	float weight;
	bool operator<(const Edge& other) const
	{
		return weight > other.weight;
	}
};

enum TSPPType
{
	FixedStartVert,
	ArbitraryStartVert
};

/* Holds resolved data for the constraints for all edges (n^2). */
typedef Bitset<2> Constraints;

/* Instead of carrying n^2 number of constraints,
 * carry the diff from the parent. Tightening is applied later on. */
struct Node
{
	/* Instead of holding all the constraints, hold the diff from parent. */
	int8_t row;
	int8_t col;
	int8_t diff_constraint;
	int8_t counter;
	Constraints constraints;
};
static_assert((1 << 8) > MAX_CITIES);
// static_assert(sizeof(Node) == sizeof(uint32_t));

/* Tighten constraints given these three rules:
 * When all but two edges adjacent to a vertex are excluded, those two edges
have to be included as otherwise it would be impossible for a tour to exist.
 * When two edges adjacent to a vertex are included, all other edges adjacent to
this vertex have to be excluded.
 * When including an edge would complete a subtour with other included edges,
this edge has to be excluded.
 * Note that a ConstraintDiff takes an 'Allowed' edge and marks it
 * whatever constraint member is. Reversing a patch is as simple as
 * marking edge as 'Allowed' again. */
struct ConstraintDiff
{
	int8_t row;
	int8_t col;
	int16_t constraint;
	bool operator<(const ConstraintDiff& other) const
	{
		const auto [_row, _col] = std::minmax(row, col);
		const auto [_orow, _ocol] = std::minmax(other.row, other.col);
		if (_row == _orow)
		{
			return _col < _ocol;
		}
		return _row < _orow;
	}
	bool operator==(const ConstraintDiff& other) const
	{
		return (row == other.row && col == other.col) || (row == other.col && col == other.row);
	}
};

void PrintConstraints(const int n, const Constraints& constraints);
inline void ApplyConstraint(const int n, Constraints* constraints, const ConstraintDiff& diff);
inline void RevertConstraint(const int n, Constraints* constraints, const ConstraintDiff& diff);

typedef std::vector<ConstraintDiff> ConstraintTighteningPatch;
void ApplyConstraintTighteningPatch(const int n,
                                    Constraints* constraints,
                                    const ConstraintTighteningPatch& patch);
void RevertConstraintTighteningPatch(const int n,
                                     Constraints* constraints,
                                     const ConstraintTighteningPatch& patch);
enum TightenConstraintsReturnCode
{
	Valid,
	ValidAndLeaf,
	NonValid
};
/* Cycle checks are done lazily during the traversal phase. */
template <TSPPType tspp_type>
TightenConstraintsReturnCode IterativelyTightenConstraints(const int n,
                                                           Constraints* constraints,
                                                           ConstraintTighteningPatch* patch,
                                                           const int start_vert);
template <TSPPType tspp_type>
TightenConstraintsReturnCode TightenConstraints(const int n,
                                                const Constraints& constraints,
                                                ConstraintTighteningPatch* patch,
                                                const int start_vert);
template TightenConstraintsReturnCode
TightenConstraints<TSPPType::ArbitraryStartVert>(const int n,
                                                 const Constraints& constraints,
                                                 ConstraintTighteningPatch* patch,
                                                 const int start_vert);
template TightenConstraintsReturnCode
TightenConstraints<TSPPType::FixedStartVert>(const int n,
                                             const Constraints& constraints,
                                             ConstraintTighteningPatch* patch,
                                             const int start_vert);

typedef std::vector<ConstraintTypes> DBGExpConstraint;
DBGExpConstraint DBGExpandConstraints(const int n, const Constraints& constraints);

float CalculateSolutionLength(const int n,
                              const Constraints& constraints,
                              const Matrix& distance_matrix);

/* Find the shortest length edge e that is allowed. */
VertexPair
CandidateEdge(const int n, const Constraints& constraints, const Matrix& distance_matrix);

template <TSPPType tspp_type>
float LowerBound(const Constraints& node, const Matrix& distance_matrix, const int start_vert);
template float LowerBound<TSPPType::FixedStartVert>(const Constraints& node,
                                                    const Matrix& distance_matrix,
                                                    const int start_vert);
template float LowerBound<TSPPType::ArbitraryStartVert>(const Constraints& node,
                                                        const Matrix& distance_matrix,
                                                        const int start_vert);

/* Solve the TSPP problem for a complete graph with given distance matrix. */
template <TSPPType tspp_type>
Path SolveMetricTSPP(const int n,
                     const Matrix& distance_matrix,
                     const int start_vert,
                     const int exec_time_limit);
template Path SolveMetricTSPP<TSPPType::ArbitraryStartVert>(const int n,
                                                            const Matrix& distance_matrix,
                                                            const int start_vert,
                                                            const int exec_time_limit);
template Path SolveMetricTSPP<TSPPType::FixedStartVert>(const int n,
                                                        const Matrix& distance_matrix,
                                                        const int start_vert,
                                                        const int exec_time_limit);

/* Bidirectional search, if a path is found over required edges,
 * insertion of the new edge (u, v) will create a sub-tour. */
bool YieldsCycle(const int n, const Constraints& constraints, int u, int v);

template <TSPPType tspp_type>
void AssertConstraintsValidity(const int n, const Constraints& constraints, const int start_vert);
template void
AssertConstraintsValidity<TSPPType::ArbitraryStartVert>(const int n,
                                                        const Constraints& constraints,
                                                        const int start_vert);
template void AssertConstraintsValidity<TSPPType::FixedStartVert>(const int n,
                                                                  const Constraints& constraints,
                                                                  const int start_vert);
Path ConstraintsToPath(const int n, const Constraints& constraints, const int start_vert);
}; // namespace branch_and_bound
