"""
Algorithms
==========

``pdppy.algorithms``

Implementations of Pickup Delivery Problem algorithms and solvers for instances.

Algorithms and solvers take in graphical instances consistent with the
description of H found in instances.request_graph and return a graphical
representation of the tour which contains the edges in the final solution
produced by the algorithm or solver. The implementations represent a series
of developed algorithms, heuristics, integer programming and linear
programming solvers for instances of the PDP.
"""

# All functions follow the same basic structure:
#     Parameters
#     ----------
#     H: graph
#         The graph of the PDP instance.
#
#     Returns
#     -------
#     P: graph
#         PDP tour solution.
#         Contains all attributes from H.
#         Contains only the edges representative of the final tour.
#         Edges have additional attribute 'value' set to 1 if the
#         edge is contained in the final tour. (All edges in P should
#         have attribute 'value' == 1.)

# TODO: Parameters and Returns not being output well in docs, restructure.

import networkx as nx
import heapq
import gurobipy as gp
from pdppy import helper

__author__ = 'Adrian Hernandez <ah695@cornell.edu>'

__all__ = ['path_build_alg', 'two_traversal_tree_alg',
           'cheapest_feasible_insertion', 'four_traversal_mst_alg',
           'five_traversal_alg', 'double_tour_alg', 'linear_prog',
           'branch_cut_integer_prog']


def path_build_alg(H):
    """
    **1-TT**:  Path Build Algorithm: Heuristic.

    Greedy path building algorithm adds feasible edges to a path starting
    with the source node.

    Algorithm starts at source node, and in a single path:

        1. Adds minimum cost edge ensuring each origin nodes is visited
        before its corresponding destination node, and

        2. Once all requests are satisfied, adds edge from *final node* to
        **t** and from **t** to **s**.
    """

    P = nx.Graph()
    P.graph.update(H.graph)
    P.graph['type'] = '1-TT'

    s = P.graph['s']
    t = P.graph['t']
    requests = P.graph['requests']

    # Set starting node to begin algorithm.
    u = s
    # Transfers node attributes from H to T.
    P.add_nodes_from([(u, H.nodes[u])])
    visited_requests = set()

    # Continue to add edges while nodes in T are not connected
    # by popping the min cost edge and attempting to add it to T.
    while P.number_of_nodes() < H.number_of_nodes():
        S = []
        # Consider edges on boundary of connected nodes in T.
        for v in nx.neighbors(H, u):
            if v not in P:
                boundary_edge = (u, v)
                wgt = H.edges[boundary_edge]['weight']
                heapq.heappush(S, (wgt, v))
        wgt, v = heapq.heappop(S)

        # Add v to P so long as it is:
        #   a destination whose origin is already in P,
        #   or t once all other nodes are in P.
        while (requests[v][1] == 'd'
               and requests[v][0] not in visited_requests) \
                or (v == t and P.number_of_nodes() < H.number_of_nodes() - 1):
            wgt, v = heapq.heappop(S)
        # Add any v that is an origin node.
        if requests[v][1] == 'o':
            visited_requests.add(requests[v][0])  # this 'o' is now visited

        # Add the feasible node and edge.
        P.add_nodes_from([(v, H.nodes[v])])
        P.add_edges_from([(u, v, H.edges[u, v])])
        # Value attribute indicates value of edge variable in solution.
        P.edges[u, v]['value'] = 1
        # Update the branching node.
        u = v

    P.add_edges_from([(t, s, H.edges[t, s])])
    P.edges[t, s]['value'] = 1

    P.graph['dist'] = P.size(weight='weight')

    return P


def two_traversal_tree_alg(H):
    """
    **2-TT**: Two Traversal Tree Algorithm: Heuristic.

    Adaptation of Prim's MST algorithm for PDP instances.
    Composes a tour from spanning tree through pre-order traversal (see
    :ref:`helper.preorder\_st\_traversal<Helper>`).

    1. Spanning tree construction:

        a. Adds the branch **(s, t)** and

        b. Adds subsequent branches greedily to ensure final tree is
        traversed a maximum of *2* times when using pre-order traversal.

    2. Traverses spanning tree using pre-order traversal.
    """

    s = H.graph['s']
    t = H.graph['t']
    requests = H.graph['requests']
    k = (H.number_of_nodes() - 2) // 2

    # Dict to store which components of each request minor {s,t,o_i,d_i}
    # are in T.
    minors = {}
    # key = request number, value = {s, t, o_i, d_i} [set] when full
    for i in range(1, (k + 1)):
        minors[i] = {s, t}

    # Initialization of spanning tree T. Contains edge (s,t).
    T = nx.Graph()
    T.graph.update(H.graph)
    T.graph['type'] = '2-TT'
    T.add_nodes_from([(s, H.nodes[s])])
    T.add_nodes_from([(t, H.nodes[t])])
    T.add_edges_from([(s, t, H.edges[s, t])])

    # Priority queue of edges leaving T to all nodes not in T.
    S = []
    # Go through all edges radiating from nodes in T to nodes not in T.
    for w in list(nx.neighbors(H, s)):
        if w not in T:
            boundary_edge = (s, w)
            wgt = H.edges[boundary_edge]['weight']
            heapq.heappush(S, (wgt, boundary_edge))
    for w in list(nx.neighbors(H, t)):
        if w not in T:
            boundary_edge = (t, w)
            wgt = H.edges[boundary_edge]['weight']
            heapq.heappush(S, (wgt, boundary_edge))

    # While not all nodes in T, pop new min. cost edge and add if feasible.
    while T.number_of_nodes() < H.number_of_nodes():
        wgt, e = heapq.heappop(S)

        if e[0] not in T or e[1] not in T:
            if e[0] not in T:
                v = e[0]
            else:
                v = e[1]

            v_req_num = requests[v][0]
            # Add tentative node v and tentative edge e with all information.
            T.add_nodes_from([(v, H.nodes[v])])
            T.add_edges_from([(*e, H.edges[e])])

            # If minor will be completed, check if addition will be feasible.
            if len(minors[v_req_num]) == 3:
                current_set = minors[v_req_num].copy()
                current_set.add(v)
                if not helper.valid_minor(T, current_set):
                    T.remove_edge(*e)
                    T.remove_node(v)
                    continue

            minors[v_req_num].add(v)
            for w in list(nx.neighbors(H, v)):
                if w not in T:
                    boundary_edge = (v, w)
                    wgt = H.edges[boundary_edge]['weight']
                    heapq.heappush(S, (wgt, boundary_edge))

    # Traverse T in pre-order s-t oriented traversal.
    P = helper.preorder_st_traversal(H, T)

    return P


def cheapest_feasible_insertion(H):
    """
    **CFI**: Cheapest Feasible Insertion Algorithm: Heuristic.

    Tour production:

        1. Finds a sub-tour *T_o* with nodes *{* **s** *} U {o_1,...,o_n}*
        (*1.5*-factor optimal tour),

        2. Iteratively adds in nodes in *N_d: {d_1,...,d_n}* through cheapest
        feasible edge insertion, and

        3. Adds edge **{s,t}** to path to produce final tour **P**.
    """

    s = H.graph['s']
    t = H.graph['t']
    requests = H.graph['requests']

    # Construct origin (G_o) sub-graph.
    G_o = nx.Graph(H)
    G_o.remove_node(t)
    d_nodes = []
    for node, req_info in requests.items():
        if req_info[1] == 'd':
            G_o.remove_node(node)
            d_nodes.append(node)

    # Perform Christofides' algorithm on graph and have s as source node.
    total_tour = helper.christofides_approx(G_o, s)

    # For each destination node, find minimum cost feasible insertion point.
    while not len(d_nodes) == 0:
        d = d_nodes.pop(0)
        visited = set()
        u_min = total_tour[len(total_tour)-2]
        v_min = total_tour[len(total_tour)-1]
        cut = -1
        c_min = (H.edges[u_min, d]['weight'] + H.edges[d, v_min]['weight']
                 - H.edges[u_min, v_min]['weight'])
        # Iterate through current tour nodes and find feasible insertion point.
        for i in range(1, (len(total_tour)-1)):
            visited.add(requests[total_tour[i]][0])
            if requests[d][0] in visited:
                u_new = total_tour[i]
                v_new = total_tour[i+1]
                c_new = (H.edges[u_new, d]['weight']
                         + H.edges[d, v_new]['weight']
                         - H.edges[u_new, v_new]['weight'])
                if c_new < c_min:
                    c_min = c_new
                    cut = i + 1
        # Insert at found point.
        total_tour = total_tour[0:cut] + [d] + total_tour[cut:]

    # Insert t between last node before s and s
    total_tour = total_tour[0:-1] + [t] + total_tour[-1:]

    # Compose final tour P from total_tour.
    P = nx.Graph()
    P.graph.update(H.graph)
    P.graph['type'] = 'CFI'
    nx.add_path(P, total_tour)

    # Include all the attributes from H in P.
    for u, v in P.edges():
        P.nodes[u].update(H.nodes[u])
        P.nodes[v].update(H.nodes[v])
        P.edges[u, v].update(H.edges[u, v])
        # Value attribute indicates value of edge variable in solution.
        P.edges[u, v]['value'] = 1

    P.graph['dist'] = P.size(weight='weight')

    return P


def four_traversal_mst_alg(H):
    """
    **4-TT**: Four Traversal Minimum Spanning Tree Algorithm: *4*-factor approx.

    Adaptation of Prim's MST algorithm for PDP instances, composes a tour
    from this MST through pre-order traversal (see
    :ref:`helper.preorder\_st\_traversal<Helper>`).
    """

    # Initialization of MST T from input graph attributes.
    T = nx.Graph()
    T.graph.update(H.graph)
    T.graph['type'] = '4-TT'

    s = T.graph['s']

    T.add_nodes_from([(s, H.nodes[s])])
    S = []
    # Go through all edges radiating from nodes in T to nodes not in T.
    for w in list(H.adj[s]):
        boundary_edge = (s, w)
        wgt = H.edges[boundary_edge]['weight']
        heapq.heappush(S, (wgt, boundary_edge))
    # While not all nodes in T, pop new min. cost edge and add if feasible.
    while T.number_of_nodes() < H.number_of_nodes():
        wgt, e = heapq.heappop(S)

        if e[0] not in T or e[1] not in T:
            if e[0] not in T:
                v = e[0]
            else:
                v = e[1]

            T.add_nodes_from([(v, H.nodes[v])])
            T.add_edges_from([(*e, H.edges[e])])

            # Add in new boundary edges by addition of v to T.
            for w in list(H.adj[v]):
                boundary_edge = (v, w)
                wgt = H.edges[boundary_edge]['weight']
                heapq.heappush(S, (wgt, boundary_edge))

    # Traverse T in pre-order s-t oriented traversal.
    P = helper.preorder_st_traversal(H, T)

    return P


def five_traversal_alg(H):
    """
    **5-TT**: Five Traversal Spanning Tree Algorithm: *5*-factor approx.

    1. Produces two sub-trees:

        a. *M_o* of nodes *{* **s** *} U {o_1,...,o_n}*

        b. *M_d* of nodes *{* **t** *} U {d_1,...,d_n}*

    2) Traverses each sub-tree in pre-order traversal, and

    3) Links the path traversals by *e = {* **s**, **t** *}*
    """

    s = H.graph['s']
    t = H.graph['t']
    requests = H.graph['requests']

    # Construct origin (G_o) and destination (G_d) sub-graphs.
    G_o = nx.Graph(H)
    G_d = nx.Graph(H)
    for node, req_info in requests.items():
        if node == s or req_info[1] == 'o':
            G_d.remove_node(node)
        else:
            G_o.remove_node(node)

    # Find MST and tour traversal of each tree.
    T_o = nx.minimum_spanning_tree(G_o)
    tour_o = list(nx.dfs_preorder_nodes(T_o, s))
    T_d = nx.minimum_spanning_tree(G_d)
    tour_d = list(nx.dfs_preorder_nodes(T_d, t))[1:]
    total_tour = tour_o + tour_d + [t, s]

    # Compose final tour P from total_tour.
    P = nx.Graph()
    P.graph.update(H.graph)
    nx.add_path(P, total_tour)

    # Include all the attributes from H in P.
    for u, v in P.edges():
        P.nodes[u].update(H.nodes[u])
        P.nodes[v].update(H.nodes[v])
        P.edges[u, v].update(H.edges[u, v])
        # Value attribute indicates value of edge variable in solution.
        P.edges[u, v]['value'] = 1

    P.graph['dist'] = P.size(weight='weight')
    P.graph['type'] = '5-TT'

    return P


def double_tour_alg(H):
    """
    **2-CHR**: Double Tour Algorithm: *4*-factor approx.

    Produces two sub-tours:

        1) *T_o* of nodes *{* **s** *} U {o_1,...,o_n}*

        2) *T_d* of nodes *{* **t** *} U {d_1,...,d_n}*

    and links these by *e = {* **s**, **t** *}*.
    """

    s = H.graph['s']
    t = H.graph['t']
    requests = H.graph['requests']

    # Construct origin: G_o and destination: G_d sub-graphs.
    G_o = nx.Graph(H)
    G_d = nx.Graph(H)
    for node, req_info in requests.items():
        if node == s or req_info[1] == 'o':
            G_d.remove_node(node)
        else:
            G_o.remove_node(node)

    # Perform Christofides' algorithm on graphs with s
    # and t as respective source nodes of each tour.
    # Remove final s -> [s, o_1, ..., o_k].
    tour_o = helper.christofides_approx(G_o, s)[0:-1]
    # Remove initial t -> [d_1, ..., d_k, t].
    tour_d = helper.christofides_approx(G_d, t)[1:]

    # Combine tours and append s. -> [s, o_i, d_i, t, s]
    total_tour = tour_o + tour_d + [s]

    # Compose final tour P from total_tour.
    P = nx.Graph()
    P.graph.update(H.graph)
    nx.add_path(P, total_tour)

    # Include all the attributes from H in P.
    for u, v in P.edges():
        P.nodes[u].update(H.nodes[u])
        P.nodes[v].update(H.nodes[v])
        P.edges[u, v].update(H.edges[u, v])
        # Value attribute indicates value of edge variable in solution.
        P.edges[u, v]['value'] = 1

    P.graph['dist'] = P.size(weight='weight')
    P.graph['type'] = '2-CHR'

    return P


def branch_cut_integer_prog(H, warm_start=None, timeout=None,
                            ignore_timeout=False):
    """
    **IP**: Branch and Cut Integer Program for the PDP.

    Parameters
    ----------
    H: graph
        :ref:`Request graph<Request (PDP) Graph>` instance.

    warm_start: list
        Starts solver with a solution for faster computation.
        Each entry is an edge *(u, v)* in the warm-start solution.
        The model will scrap the start if found to be infeasible.

    timeout: int
        Maximum computation time (sec.) for model.

    ignore_timeout: bool
        Ignore timeout until first solution found.
        ``False`` (default) - timeout parameter immediately enforced, may find
        no solution.

    Returns
    -------
    P: graph
        :ref:`PDP tour<Tour Graph>` solution.

    Raises
    ------
    GurobiError
        If no solution is produced before model runtime times out.
    """

    s = H.graph['s']
    t = H.graph['t']

    # Model initialization.
    m = gp.Model()
    m.setParam('OutputFlag', 0)   # Silences output.

    # If timeout is provided and ignore_timeout is True
    # make model find at least 1 solution before timing out.
    if timeout is not None:
        if ignore_timeout:
            old_solution_limit = m.Params.SolutionLimit
            m.Params.SolutionLimit = 1  # Minimum 1 solution if timeout ignored.
        else:
            m.Params.TimeLimit = timeout    # timeout enforced from start.

    m.Params.LazyConstraints = 1    # Enables cutting planes method.

    # Addition of edge weights to costs.
    costs = [e[2]['weight'] for e in H.edges(data=True)]

    # Addition of edge decision variables.
    # Corresponding costs are coefficients in the objective function.
    mvars = m.addVars(H.edges(), obj=costs, vtype=gp.GRB.BINARY, name='x')
    # Key order is interchangeable since H is undirected.
    for i, j in mvars.keys():
        mvars[j, i] = mvars[i, j]
        # Set initial values to the warm_start values
        if warm_start is not None:
            if (i, j) in warm_start or (j, i) in warm_start:
                mvars[i, j].Start = 1.0
            else:
                mvars[i, j].Start = 0.0

    # Addition of constraints.
    # Each node must have 1 inbound and 1 outbound edge.
    m.addConstrs(mvars.sum(i, '*') == 2 for i in H.nodes())
    # (s,t) edge must exist.
    m.addConstr(mvars[s, t] == 1)

    # Storing values in created model fields.
    m._vars = mvars
    m._H = H

    # Update changes and call the branch and cut method in helper.py.
    m.update()
    m.optimize(helper.separation_oracle)

    # If timeout is provided and 1 solution has been found,
    # continue optimization if there is more time left.
    if timeout is not None and ignore_timeout:
        runtime = m.getAttr(gp.GRB.Attr.Runtime)
        if runtime < timeout:
            m.Params.TimeLimit = timeout - runtime
            m.Params.SolutionLimit = old_solution_limit
            # Continue optimization.
            m.update()
            m.optimize(helper.separation_oracle)

    # Retrieve model solution and variable and values.
    # Raise GurobiError and return None if unsuccessful.
    try:
        opt_bc = m.objval
        x_bc = m.getAttr('x', mvars).items()

        # Construct tour from model solution through method in helper.py.
        P = helper.construct_tour(H, x_bc)
        P.graph['dist'] = opt_bc
        P.graph['type'] = 'IP'

        return P

    except gp.GurobiError:
        print('Model timed out before any solution was found.')

        return

    except AttributeError:
        print('Model was found to be infeasible.')

        return


def linear_prog(H):
    """
    **LP**: Linear Program Relaxation for the PDP.

    Used as a lower bound for comparison of algorithmic performance.
    """

    s = H.graph['s']
    t = H.graph['t']

    # Model initialization.
    m = gp.Model()
    m.setParam('OutputFlag', 0)  # Silences output.
    m.Params.LazyConstraints = 1  # Enables cutting planes method.
    # Addition of edge weights to costs.
    costs = [e[2]['weight'] for e in H.edges(data=True)]

    # Addition of edge decision variables.
    # Corresponding costs are coefficients in the objective function.
    mvars = m.addVars(H.edges(), obj=costs, vtype=gp.GRB.CONTINUOUS,
                      lb=0.0, ub=1.0, name='x')
    # Key order is interchangeable since H is undirected.
    for i, j in mvars.keys():
        mvars[j, i] = mvars[i, j]

    # cbLazy in Gurobi only works for a mixed integer linear program:
    # add dummy binary variable to LP model.
    v = m.addVar(obj=0, vtype=gp.GRB.BINARY, lb=0.0, ub=0.0)

    # Addition of constraints.
    # Each node must have 1 inbound and 1 outbound edge.
    m.addConstrs(mvars.sum(i, '*') == 2 for i in H.nodes())
    # (s,t) edge must exist.
    m.addConstr(mvars[s, t] == 1)

    # Storing values in created model fields.
    m._vars = mvars
    m._H = H

    # Update changes and call the branch and cut method in helper.py.
    m.update()
    m.optimize(helper.separation_oracle)

    # Remove dummy binary variable form model.
    m.remove(v)
    # Update model without dummy variable.
    m.update()
    m.optimize()

    # Retrieve model solution and variable and values.
    opt_lp = m.objval
    x_lp = m.getAttr('x', mvars).items()

    # Construct tour from model solution through method in helper.py.
    # Parameter integer=False means P is the graphical support of the solution.
    P = helper.construct_tour(H, x_lp, integer=False)
    P.graph['dist'] = opt_lp
    P.graph['type'] = 'LP'

    return P

