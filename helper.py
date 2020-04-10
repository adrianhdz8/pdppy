"""
Helper
======

``pdppy.helper``

Helper file contains all necessary helper functions for package.
The functions are separated by module in which they are used.
"""

import networkx as nx
import osmnx as ox
import numpy as np
import gurobipy as gp
from itertools import combinations
from datetime import datetime

__author__ = 'Adrian Hernandez <ah695@cornell.edu>'

__all__ = ['choose_s_t_requests', 'random_geo_graph_generator',
           'city_graph_generator',
           'christofides_approx', 'preorder_st_traversal', 'valid_minor',
           'contract_tree', 'separation_oracle', 'construct_tour',
           'node_labels', 'node_marker', 'graph_routes', 'plot_routes']


# Instances Module

def choose_s_t_requests(k, nodes, seed):
    """
    Choose **s**, **t**, and request pairs pseudo-randomly from given nodes.

    Parameters
    ----------
    k: int
        *2 **k** *+ 2* is the number of nodes in the request graph.

    nodes: list
        Nodes to select from.

    seed: int
        Seed for random choice of nodes.

    Returns
    -------
    s: node
        Source node.

    t: node
        Target node.

    request_pairs: list
        Request pairs from **G**.
        Each request_pair entry is (origin, destination).
        No node can be repeated in multiple request pairs.
    """

    # Pick random s and t from nodeset.
    np.random.seed(seed=seed)
    s_t = np.random.choice(nodes, 2, replace=False)
    s = s_t[0]
    t = s_t[1]
    nodes.remove(s)
    nodes.remove(t)
    # Pick random request pairs from remaining nodeset.
    requests = np.random.choice(nodes, (k, 2), replace=False)
    request_pairs = [(r[0], r[1]) for r in requests]

    return s, t, request_pairs


def random_geo_graph_generator(k, seed=None):
    """
    Generate a random geometric graph with *2 **k** *+ 2* nodes.

    Graph is generated in the unit square with random node positions. The
    graph is undirected, made complete (metric closure), and edges are
    weighted by their Euclidian distance.

    Parameters
    ----------
    k:  int
        *2 **k** *+ 2* is the number of nodes in graph.

    seed: int
        Seed for random graph generator.
        System time taken as default.

    Returns
    -------
    G: graph
        Randomly generated graph on unit square of 2*k+2 nodes.
    """

    # Generate seed if none provided.
    if seed is None:
        time = datetime.now()
        seed = 1000000 * time.second + time.microsecond

    # Generate random graph with seed and store seed in attribute.
    G = nx.random_geometric_graph((2 * k + 2), radius=5, dim=2, seed=seed)
    G.graph['seed'] = seed
    G.graph['name'] = 'random_geo'

    for u in G.nodes():
        G.nodes[u]['x'] = G.nodes[u]['pos'][0]
        G.nodes[u]['y'] = G.nodes[u]['pos'][1]

    # Complete metric closure of G, using Euclidean distance as weight.
    for u, v in combinations(list(G.nodes()), 2):
        wgt = np.sqrt((G.nodes[u]['pos'][0] - G.nodes[v]['pos'][0]) ** 2 +
                      (G.nodes[u]['pos'][1] - G.nodes[v]['pos'][1]) ** 2)
        G.edges[u, v]['weight'] = wgt

    return G


def city_graph_generator(city):
    """
    Generate an OSMnx graph of **city**.

    Graph modified to be undirected, weighted, and meeting the triangle
    inequality.

    Parameters
    ----------
    city: str
        The city to query road network from. *'City, Country'*.

    Returns
    -------
    G: graph
        OSMnx generated graph of **city**.
        Graph is weighted, undirected, and satisfies the triangle inequality.

    Raises
    ------
    AssertionError
        If **city** is not a string.
    """

    # Assert city input is string.
    assert type(city) is str, "City name not valid, must be string."

    # Query and generate OSMNX graph from city.
    O = ox.graph_from_place(city, network_type='drive', simplify=True,
                            retain_all=False, timeout=60)

    # Keep only the largest strongly connected component.
    keep = max(nx.kosaraju_strongly_connected_components(O), key=len)
    G_sub = O.subgraph(keep).copy()
    G_dir = nx.DiGraph(G_sub)
    G = nx.Graph()
    G.graph.update(G_dir.graph)

    # Iterate over edges to keep maximum weight edges
    # in the case where there are parallel edges.
    for u, v in G_dir.edges():
        G.add_nodes_from([(u, G_dir.nodes[u]), (v, G_dir.nodes[v])])
        if (v, u) in G_dir.edges():
            wgt = max(G_dir.edges[u, v]['length'], G_dir.edges[v, u]['length'])
        else:
            wgt = G_dir.edges[u, v]['length']
        G.add_edges_from([(u, v, G_dir.edges[u, v])])
        G.edges[u, v]['length'] = wgt
        G.edges[u, v]['weight'] = wgt

    return G


# Algorithms Module


def christofides_approx(G, source):
    """
    Christofides' Algorithm: *1.5*-factor TSP approximation
        Author: Juan Carlos Martinez Mori `jm2638`

    Valid for metric traveling salesman problem.

    Parameters
    ----------
    G: graph
        The graph to apply algorithm to.

    source: node
        Node to be taken as root of MST.

    Returns
    -------
    tour: list
        Tour found by algorithm. *[* **source** *, nodes,* **source** *]*.
    """

    # Find MST of G.
    mst = nx.minimum_spanning_tree(G)

    # Compose sub-graph of nodes with odd parity.
    T = nx.subgraph(G, [p for p, deg in nx.degree(mst) if deg % 2 == 1])

    # Find minimum cost t-join.
    for p, q, attr in T.edges(data=True):
        attr['weight'] *= -1
    T_join = nx.max_weight_matching(T, maxcardinality=True, weight='weight')

    # Add minimum cost t-join to the MST.
    mst.add_edges_from(T_join)

    # Traverse and shortcut the MST and minimum cost t-join.
    tour = list(nx.dfs_preorder_nodes(mst, source))
    tour.append(source)

    return tour


def preorder_st_traversal(H, T):
    """
    Conducts a pre-ordered traversal of tree **T** with an **s-t** orientation.

    **s** is the source node and **t** is the target node for the PDP.
    Traversal of an **s-t** oriented tree involves traversing the tree in
    pre-order, but leaving the branch with **t** to be visited last at
    each step possible.

    Parameters
    ----------
    H: graph
        The :ref:`request graph<Request (PDP) Graph>`.

    T: graph
        A spanning tree of **H**.

    Returns
    -------
    P: graph
        :ref:`PDP tour<Tour Graph>` solution.
    """

    s = T.graph['s']
    t = T.graph['t']
    requests = T.graph['requests']

    P = nx.Graph()
    P.graph.update(T.graph)

    # Find successors and predecessors of nodes in T for traversing.
    V = set(T.nodes())
    T_successors = nx.dfs_successors(T, s)
    successors = set(T_successors.keys())
    # Update for nodes with empty successors.
    [T_successors.update({j: []}) for j in V.difference(successors)]
    T_predecessors = nx.dfs_predecessors(T, s)
    # Branch with t must be visited last. Store the nodes on this branch
    # for identification of traversal priority.
    last_branch = {t}
    b = T_predecessors[t]
    while not b == s:
        last_branch.add(b)
        b = T_predecessors[b]

    # boundary_nodes will keep track of the next nodes to be visited.
    boundary_nodes = []
    # If a node is visited and added to P, its request number is added.
    visited_requests = {requests[s][0]}
    # Initialize with source, s.
    P.add_nodes_from([(s, T.nodes[s])])
    u = s
    # While P does not contain all nodes, traverse T and add every
    # feasible node in order.
    while P.number_of_nodes() < T.number_of_nodes():
        v = s
        neighbors = T_successors[v].copy()
        # For each node w branching from v, determine if w is to be
        # visited last, if so, relocate to front of boundary_nodes
        # so that it is popped (visited) last.
        for w in neighbors:
            if w in last_branch:
                boundary_nodes = [w] + boundary_nodes
            else:
                boundary_nodes = boundary_nodes + [w]
        # Consider nodes while there are nodes in boundary_nodes.
        while len(boundary_nodes) > 0:
            v = boundary_nodes.pop()
            # If v has not been added to P, and:
            #   a) v is an origin type,
            #   b) v is a destination type with corresponding origin in P, or
            #   c) v is t and all other nodes are in P
            # Then add v to P and make an edge from last added u node.
            if v not in P and (requests[v][1] == 'o'
                               or (requests[v][1] == 'd'
                                   and requests[v][0] in visited_requests)
                               or (v == t and P.number_of_nodes()
                                   == T.number_of_nodes() - 1)):
                visited_requests.add(requests[v][0])
                P.add_nodes_from([(v, T.nodes[v])])
                P.add_edges_from([(u, v, H.edges[u, v])])
                u = v
            neighbors = T_successors[v].copy()
            # Update the boundary nodes.
            for w in neighbors:
                if w in last_branch:
                    boundary_nodes = [w] + boundary_nodes
                else:
                    boundary_nodes = boundary_nodes + [w]

    # Add in final edge from t to s to complete cycle.
    P.add_edges_from([(t, s, H.edges[t, s])])
    P.graph['dist'] = P.size(weight='weight')

    return P


def valid_minor(T, node_set):
    """
    Helper function for :ref:`two\_tree\_traversal<Algorithms>`.

    Determines if **node_set** produces an **s-t** oriented sub-tree of **T**.

    Parameters
    ----------
    T: graph
        The current spanning tree of the PDP instance.

    node_set: set
        Set of nodes to be in sub-tree of **T**.

    Returns
    -------
    bool
        ``False`` if sub-tree satisfies:

            a. all *3* edges of subset contain **s**, or

            b. all *3* edges of subset contain **t**, or

            c. subset contains both edges **(s, d)** and **(t, o)**

           ``True`` otherwise
    """

    s = T.graph['s']
    t = T.graph['t']
    requests = T.graph['requests']

    # For all nodes in node_set, classify by request type.
    for u in node_set:
        if requests[u][1] == 'o':  # origin node
            o = u
        elif requests[u][1] == 'd':  # destination node
            d = u

    # Call contract tree method to keep only relevant requests.
    sub_T = contract_tree(T, node_set)

    # Any of these 3 edge combinations leads to an unfeasible sub tree for POST.
    if ((s, o) in sub_T.edges and (s, d) in sub_T.edges
            and (s, t) in sub_T.edges):
        return False
    elif ((t, o) in sub_T.edges and (t, d) in sub_T.edges
            and (t, s) in sub_T.edges):
        return False
    elif (s, d) in sub_T.edges and (t, o) in sub_T.edges:
        return False
    else:
        return True


def contract_tree(T, node_set):
    """
    Helper function for :ref:`valid_minor<Helper>`.

    Contract the tree **T** to contain only the nodes in **node_set**:
        1. contract all edges containing no nodes in **node_set**,
        2. contract all edges containing *s* and a node not in **node_set**,
        3. contract all edges containing *o* and a node not in **node_set**,
        4. contract all edges containing *d* and a node not in **node_set**, and
        5. contract all edges containing *t* and a node not in **node_set**.

    Parameters
    ----------
    T: graph
        The current spanning tree of the PDP instance.

    node_set: set
        Set of nodes to be in sub-tree of **T**.

    Returns
    -------
    T_con: graph
        Contraction of **T** which contains only nodes in **node_set**.
    """

    s = T.graph['s']
    t = T.graph['t']
    requests = T.graph['requests']

    # For all nodes in node_set, classify by type.
    for u in node_set:
        if requests[u][1] == 'o':  # origin node
            o = u
        elif requests[u][1] == 'd':  # destination node
            d = u

    T_con = nx.Graph(T)

    # s-wise contraction.
    while not set(nx.neighbors(T_con, s)).issubset(node_set):
        for edge in list(T_con.edges(s)):
            u, v = edge  # u is always input, s in this case
            if v not in node_set:
                T_con = nx.contracted_nodes(T_con, u, v, self_loops=False)

    # o-wise contraction.
    while not set(nx.neighbors(T_con, o)).issubset(node_set):
        for edge in list(T_con.edges(o)):
            u, v = edge
            if v not in node_set:
                T_con = nx.contracted_nodes(T_con, u, v, self_loops=False)

    # d-wise contraction.
    while not set(nx.neighbors(T_con, d)).issubset(node_set):
        for edge in list(T_con.edges(d)):
            u, v = edge
            if v not in node_set:
                T_con = nx.contracted_nodes(T_con, u, v, self_loops=False)

    # t-wise contraction.
    while not set(nx.neighbors(T_con, t)).issubset(node_set):
        for edge in list(T_con.edges(t)):
            u, v = edge
            if v not in node_set:
                T_con = nx.contracted_nodes(T_con, u, v, self_loops=False)

    return T_con


def separation_oracle(m, where):
    """
    Helper function for :ref:`branch\_cut\_integer\_prog and
    linear\_prog<Algorithms>`.

    Performs the branching and cutting planes component for optimization model.

    Parameters
    ----------
    m: Gurobi model
        The Gurobi model being optimized.

    where:
        Location of Gurobi in optimization process.
    """

    if where == gp.GRB.Callback.MIPSOL:

        H = m._H  # retrieve graph from model field
        s = H.graph['s']
        t = H.graph['t']
        requests = H.graph['requests']

        orig_requests = dict()
        dest_requests = dict()
        for node, req_info in requests.items():
            if req_info[1] == 'o':
                orig_requests.update({req_info[0]: node})
            elif req_info[1] == 'd':
                dest_requests.update(({req_info[0]: node}))

        # set current solution as edge capacity
        values = m.cbGetSolution(m._vars)
        Q = H.copy()
        for (i, j) in m._vars.keys():
            Q[i][j]['capacity'] = values[i, j]

        # subtour elimination: check flow between u and every other node
        u = s
        for v in set(Q.nodes()).difference({u}):

            cut_value, (S, _) = nx.minimum_cut(Q, u, v, capacity='capacity')

            # add constraint to model if violated
            if cut_value < 2:
                cut = list(nx.edge_boundary(Q, S))
                m.cbLazy(gp.quicksum(m._vars[i, j] for (i, j) in cut) >= 2)

        # precedence respect: check flow for each request
        for req in orig_requests.keys():
            # add in proxy nodes to make cuts
            Q.add_edge(s, 's_prime')
            Q.add_edge(dest_requests[req], 's_prime')
            Q.add_edge(t, 't_prime')
            Q.add_edge(orig_requests[req], 't_prime')

            # obtain cut value and node partition
            cut_value, (S, _) = nx.minimum_cut(Q, 's_prime', 't_prime',
                                               capacity='capacity')
            Q.remove_nodes_from(['s_prime', 't_prime'])
            S -= {'s_prime', 't_prime'}

            # add constraint to model if violated
            if cut_value < 4:
                cut = list(nx.edge_boundary(Q, S))
                m.cbLazy(gp.quicksum(m._vars[i, j] for (i, j) in cut) >= 4)


def construct_tour(H, x_opt, integer=True):
    """
    Construct tour graph from decision variables from optimization.

    Tour graph **P** will have attribute 'value' with decision variable values.
    Linear programming solutions may have non-integer values in 'value'.
    All other solutions will only have 'value' of *1*.

    Parameters
    ----------
    H: graph
        :ref:`Request graph<Request (PDP) Graph>` of instance.

    x_opt: list
        Optimization model decision variables and values.
        Represents the edges in the solution set.

    integer: bool
        ``True`` if integer programming model.
        ``False`` if linear programming model.

    Returns
    -------
    P: graph
        :ref:`PDP tour<Tour Graph>` solution.
    """

    P = nx.Graph()
    P.graph.update(H.graph)

    for x in x_opt:
        val = x[1]
        if integer:
            val = round(val)
        if val is not 0:
            u, v = x[0]
            P.add_nodes_from([(u, H.nodes[u])])
            P.add_nodes_from([(v, H.nodes[v])])
            P.add_edges_from([(u, v, H.edges[u, v])])
            if 'capacity' in P.edges[u, v]:
                del P.edges[u, v]['capacity']
            P.edges[u, v].update({'value': val})

    return P


# Plot Module


def node_labels(dict_requests):
    """
    Create dict of nodes with node labels for plotting request number and type.

    Parameters
    ----------
    dict_requests: dict
        Dictionary of requests with request information.

    Returns
    -------
    labels: dict
        Node with label strings for plotting.
    """

    labels = dict()
    for key, val in dict_requests.items():
        if val[1] == 's' or val[1] == 't':
            labels[key] = '$' + val[1] + '$'
        else:
            # fix to be val[1] + '_' + str(val[0])
            # see if used in plot module
            labels[key] = '$' + val[1] + '_' + str(val[0]) + '$'

    return labels


def node_marker(H, node):
    """
    Return **node_marker** for plotting of nodes.

    Parameters
    ----------
    H: graph
        :ref:`Request graph<Request (PDP) Graph>` of instance.

    node: node
        Node to get marker from.

    Returns
    -------
    str
        Marker for plotting:
        source = 'c' = circle
        target = 's' = square
        origin = ^ = triangle-up
        destination = v = triangle-down
    """

    s = H.graph['s']
    t = H.graph['t']
    requests = H.graph['requests']

    if node == s:
        return 'o'      # circle
    elif node == t:
        return 's'      # square
    elif requests[node][1] == 'o':
        return '^'      # triangle-up
    elif requests[node][1] == 'd':
        return 'v'      # triangle-down


def graph_routes(G, P):
    """
    Find routes and position of nodes for each edge in **P**.

    The routes of **P** are returned as the shortest path route on **G**
    between each edge.

    Parameters
    ----------
    G: graph
        :ref:`Input graph<Input Graph>` from OSMnx query.

    P: graph
        :ref:`PDP tour<Tour Graph>` solution.

    Returns
    -------
    Groutes: list
        Routes from **G** for each edge in **P**.

    od_lats: list
        Latitude of origin-destination pairs in **P**.

    od_lons: list
        Longitude of origin-destination pairs in **P**.

    markers: list
        Marker shapes of nodes based on request type.

    labels: list
        Marker labels of nodes from request number and type information.
    """

    s = P.graph['s']
    t = P.graph['t']
    requests = P.graph['requests']

    Groutes = []
    od_lats = []
    od_lons = []
    markers = []
    labels = []
    nodes = set()

    # For each edge pair in P: find the route, position, and label information.
    for u, v in list(P.edges):
        u_original = P.nodes[u]['original']
        v_original = P.nodes[v]['original']
        # Do not plot edge between s and t for clarity on plot.
        # if not (u, v) == (s, t) and not (v, u) == (s, t):
        Groutes.append(nx.shortest_path(G, u_original, v_original,
                                        weight='weight'))

        # For u and v, store position (x, y), and marker and label info.
        if u not in nodes:
            nodes.add(u)
            od_lats.append(G.nodes[u_original]['y'])
            od_lons.append(G.nodes[u_original]['x'])
            markers.append(node_marker(P, u))
            if requests[u][0] == 0:
                label = requests[u][1]
            else:
                label = requests[u][1] + '_' + str(requests[u][0])
            labels.append(label)

        if v not in nodes:
            nodes.add(v)
            od_lats.append(G.nodes[v_original]['y'])
            od_lons.append(G.nodes[v_original]['x'])
            markers.append(node_marker(P, v))
            if requests[v][0] == 0:
                label = requests[v][1]
            else:
                label = requests[v][1] + '_' + str(requests[v][0])
            labels.append(label)

    return Groutes, od_lats, od_lons, markers, labels


def plot_routes(G, Groutes):
    """
    Produce overlain plot of **G** and **Groutes** with specified style below.

    Parameters
    ----------
    G: graph
        :ref:`Input graph<Input Graph>` from OSMnx query.

    Groutes: list
        Routes from **G** for each edge in **P**.

    Returns
    -------
    fig: figure object
        Figure handle for plotting

    ax: axes object
        Axes handle for plotting
    """

    fig, ax = ox.plot_graph_routes(G, Groutes, route_color='r',
                                   route_alpha=0.8, route_linewidth=2,
                                   node_size=0, edge_color='k',
                                   edge_linewidth=0.5, edge_alpha=0.5,
                                   orig_dest_node_size=0,
                                   orig_dest_node_alpha=0.0,
                                   fig_height=9, fig_width=6)

    return fig, ax
