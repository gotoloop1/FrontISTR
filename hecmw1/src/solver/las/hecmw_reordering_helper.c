#include "lib.h"
#include "reordering.h"

void single_sort_wrap(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    single_sort(*l, x, *max_val);
}
void single_sort_wrap_(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    single_sort(*l, x, *max_val);
}
void single_sort_wrap__(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    single_sort(*l, x, *max_val);
}
void SINGLE_SORT_WRAP(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    single_sort(*l, x, *max_val);
}

void first_element_sort_wrap(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    first_element_sort(*l, x, *max_val);
}
void first_element_sort_wrap_(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    first_element_sort(*l, x, *max_val);
}
void first_element_sort_wrap__(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    first_element_sort(*l, x, *max_val);
}
void FIRST_ELEMENT_SORT_WRAP(
    const int *l,
    int *restrict x,
    const int *max_val
) {
    first_element_sort(*l, x, *max_val);
}

void reverse_cuthill_mckee_wrap(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    reverse_cuthill_mckee(*n, edge_separator, edge, order);
}
void reverse_cuthill_mckee_wrap_(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    reverse_cuthill_mckee(*n, edge_separator, edge, order);
}
void reverse_cuthill_mckee_wrap__(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    reverse_cuthill_mckee(*n, edge_separator, edge, order);
}
void REVERSE_CUTHILL_MCKEE_WRAP(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    reverse_cuthill_mckee(*n, edge_separator, edge, order);
}

void edge_based_cost_minimization_wrap(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    edge_based_cost_minimization(*n, edge_separator, edge, order);
}
void edge_based_cost_minimization_wrap_(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    edge_based_cost_minimization(*n, edge_separator, edge, order);
}
void edge_based_cost_minimization_wrap__(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    edge_based_cost_minimization(*n, edge_separator, edge, order);
}
void EDGE_BASED_COST_MINIMIZATION_WRAP(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    edge_based_cost_minimization(*n, edge_separator, edge, order);
}

void node_based_cost_minimization_wrap(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    node_based_cost_minimization(*n, edge_separator, edge, order);
}
void node_based_cost_minimization_wrap_(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    node_based_cost_minimization(*n, edge_separator, edge, order);
}
void node_based_cost_minimization_wrap__(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    node_based_cost_minimization(*n, edge_separator, edge, order);
}
void NODE_BASED_COST_MINIMIZATION(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    node_based_cost_minimization(*n, edge_separator, edge, order);
}

void hybrid_cost_minimization_wrap(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void hybrid_cost_minimization_wrap_(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void hybrid_cost_minimization_wrap__(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void HYBRID_COST_MINIMIZATION_WRAP(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    hybrid_cost_minimization(*n, edge_separator, edge, order);
}

void coarsen_refine_hybrid_cost_minimization_wrap(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    coarsen_refine_hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void coarsen_refine_hybrid_cost_minimization_wrap_(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    coarsen_refine_hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void coarsen_refine_hybrid_cost_minimization_wrap__(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    coarsen_refine_hybrid_cost_minimization(*n, edge_separator, edge, order);
}
void COARSEN_REFINE_HYBRID_COST_MINIMIZATION_WRAP(
    const int *n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict edge_weight,
    int *restrict order
) {
    coarsen_refine_hybrid_cost_minimization(*n, edge_separator, edge, order);
}
