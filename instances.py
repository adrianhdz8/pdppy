"""
Instances
=========

``pdppy.instances``

Graphical PDP instance management for a functional application of algorithms.

Prior to use of the ``algorithms.py`` methods, any graphical instance, given or
generated, must satisfy the :ref:`request graph<Request (PDP) Graph>`
criteria: metrically closed and contains the required PDP information stored
as graph attributes.
"""

import networkx as nx
from itertools import combinations
from pdppy import helper
from datetime import datetime

__author__ = 'Adrian Hernandez <ah695@cornell.edu>'

# List of all function names to be imported for when import * is called.
__all__ = ['request_graph',
           'random_geo_graph',
           'city_graph']

# TODO: Include code examples in docstring.
# TODO: Shorten parameter descriptions here. Use reference to other page of
#  explanations if necessary.


def request_graph(G, s, t, request_pairs):
    """
    Generates PDP version of input graph **G**.

    Output graph **H** contains additional information on the source, target,
    and request nodes and is metrically closed by use of the triangle
    inequality.

    Parameters
    ----------
    G: NetworkX graph
        The original NetworkX :ref:`input graph<Input Graph>`.
        Graph is weighted, undirected, and satisfies the triangle inequality.

    s: node
        Source node chosen from **G**.

    t: node
        Target node chosen from **G**.

    request_pairs: list
        Request pairs chosen from **G**.
        Each **request_pair** entry is *(origin, destination)*.
        No node can be repeated in multiple entries.

    Returns
    -------
    H: NetworkX graph
        :ref:`Request graph<Request (PDP) Graph>` version of **G**.

    Raises
    ------
    AssertionError
        If more nodes are given than exist in **G**.
        If **s**, **t**, or request nodes are not contained in **G**.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_weighted_edges_from([('ORD', 'LAX', 800), ('ORD', 'ATL', 500),
    >>> ...                        ('ATL', 'MIA', 400), ('MIA', 'DFW', 350),
    >>> ...                        ('LAX', 'SFO', 200), ('DFW', 'LAX', 500)])
    >>> H = request_graph(G, 'MIA', 'ORD', [('SFO', 'LAX'), ('DFW', 'ATL')])
    """

    V = set(G.nodes())
    k = len(request_pairs)

    assert len(V) >= (2*k + 2), "Insufficient nodes in graph."
    assert s in V, str(s) + " is not a node in G."
    assert t in V, str(t) + " is not a node in G."

    s_new = 0
    t_new = 2*k+1
    # Store requests by: key = new node name (using counter from 0 to 2*k+1),
    # value = (request number, request type).
    dict_requests = {s_new: (0, 's'), t_new: (0, 't')}
    # Store original node names by: key = new node name (counter),
    # value = original node name in G.
    dict_original = {s_new: s, t_new: t}
    i = 1   # Request counter.
    n = 1   # New node counter.

    # Iterate over all given request_pairs.
    for orig, dest in request_pairs:
        assert orig in V, str(orig) + ' is not a node in G.'
        assert dest in V, str(dest) + ' is not a node in G.'
        # Add in new node name along with request information.
        dict_requests[n] = (i, 'o')
        dict_requests[n+1] = (i, 'd')
        # Add in new node name with original node name from G.
        dict_original[n] = orig
        dict_original[n+1] = dest
        i += 1
        n += 2

    # H is a sub-graph of G containing s, t, and request_pairs;
    # this action preserves all attributes of G.
    H = nx.Graph()
    H.graph.update(G.graph)
    for node in dict_original.keys():
        original = dict_original[node]
        H.add_node(node)
        H.nodes[node].update(G.nodes[original])
        H.nodes[node]['original'] = original

    # Store new node information in attributes of H.
    H.graph['s'] = s_new
    H.graph['t'] = t_new
    H.graph['requests'] = dict_requests

    # Complete metric closure of H.
    for u, v in combinations(list(H.nodes()), 2):
        u_original = dict_original[u]
        v_original = dict_original[v]
        H.add_edge(u, v,
                   weight=nx.shortest_path_length(G, u_original, v_original,
                                                  weight='weight'))

    return H


def random_geo_graph(k, seed=None):
    """
    Generates a random geometric graph in the unit square.

    PDP graph contains **k** request pairs and an **s** and **t** node selected
    randomly through the **seed** parameter.

    Parameters
    ----------
    k: int
        *2* **k** *+ 2* is the number of nodes in graph.

    seed: int
        Seed for random graph generator.
        System time taken as default.
        **H** will always be the same, provided same **k** and **seed**
        parameters.

    Returns
    -------
    G: graph
        Randomly generated :ref:`input graph<Input Graph>` on unit square of
        *2* **k** *+ 2* nodes.

    H: graph
        :ref:`Request graph<Request (PDP) Graph>` version of G.

    Examples
    --------
    Can supply a seed.

    >>> G, H = random_geo_graph(3, 10001)

    Can have computer generate a seed; stored in H.graph['seed'].

    >>> G, H = random_geo_graph(4)
    """

    # Call random_geo_graph_generator for generation of G.
    G = helper.random_geo_graph_generator(k, seed=seed)
    V = list(G.nodes())
    seed = G.graph['seed']

    s, t, request_pairs = helper.choose_s_t_requests(k, V, seed=seed)

    H = request_graph(G, s, t, request_pairs)

    return G, H


def city_graph(city, G=None, k=None, seed=None, request_nodes=None):
    """
    Generate an OSMnx graph of **city**.

    PDP graph contains **k** request pairs and an **s** and **t** node selected
    randomly from original OSMnx graph through the **seed** parameter.

    Parameters
    ----------
    city: str
        The city to query road network from. *'City, Country'*.

    k: int
        *2* **k** *+ 2* is the number of nodes in graph.
        If ``None``, default is *3* request pairs.

    seed: int
        Seed for random graph generator.
        System time taken as default.
        **H** will always be the same, provided same **G** or city graph,
        **k**, and **seed** parameters.

    request_nodes: list
        **k** *+ 1* pairs of origin-destination elements to include from
        **G** in **H**.
        **(s, t)** pair is *0th* term and **(o_i, d_i)** pairs are remaining
        origin-destination pairs.

    Returns
    -------
    G: graph
        OSMnx generated :ref:`input graph<Input Graph>` of **city**.

    H: graph
        :ref:`Request graph<Request (PDP) Graph>` version of G.

    Examples
    --------
    Query a city's graph and make request graph a subgraph on *3* random
    request pairs set by the input seed.

    >>> G, H = city_graph('Manhattan, USA', k=3, seed=23451)

    Can also supply OSMnx graph **G** as a parameter in the case certain
    modifications were made by the user.

    >>> import osmnx as ox
    >>> G = ox.graph_from_place('Cornell University, USA', network_type='drive')
    >>> G, H = city_graph('Cornell University, USA', G=G, k=5)
    """

    # If no k value supplied, use 3 request nodes.
    if k is None:
        k = 3
    # If no seed value supplied, generate from system time.
    if seed is None:
        time = datetime.now()
        seed = time.hour * 10000 + time.minute * 100 + time.second

    # Call city_graph_generator for generation of G.
    if G is None:
        G = helper.city_graph_generator(city)
    if request_nodes is None:
        V = list(G.nodes())
        s, t, request_pairs = helper.choose_s_t_requests(k, V, seed=seed)
    else:
        V = request_nodes
        s, t = V[0]
        request_pairs = V[1:]

    G.graph['seed'] = seed

    H = request_graph(G, s, t, request_pairs)
    H.graph['seed'] = seed

    return G, H



