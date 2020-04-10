"""
Plot
====

``pdppy.plot``

Plotting and visualization method for algorithm and solver solutions.
"""

import networkx as nx
import matplotlib
from pdppy import helper

__all__ = ['plot_tour']


def plot_tour(P, G=None):
    """
    Plot tour **P** alone or overlay on **G** if provided.

    **G** must be an OSMnx graph.

    Parameters
    ----------
    P: graph
        :ref:`PDP tour<Tour Graph>` solution.

    G: graph
        :ref:`Input<Input Graph>` or :ref:`request<Request (PDP) Graph>` graph
        to overlay solution over.
    """

    # If G is provided and if it is an OSMNX graph, overlay tour P on it.
    if G is not None and 'crs' in G.graph.keys():
        G = nx.MultiDiGraph(G)
        num_reqs = (len(P.graph['requests']) - 2) // 2

        # Find shortest path routes on G from tour P.
        Groutes, od_lats, od_lons, markers, labels = helper.graph_routes(G, P)

        # Overlay these routes on plot of G and style with helper method.
        fig, ax = helper.plot_routes(G, Groutes)

        # Create a colormap for the markers of the nodes and plot on locations.
        norm = matplotlib.colors.Normalize(vmin=1, vmax=num_reqs)
        cmap = matplotlib.cm.get_cmap('Spectral')

        for i in range(len(od_lons)):
            if len(labels[i]) > 1:
                req = int(labels[i][-1])
                color = cmap(norm(req))
            else:
                color = 'turquoise'

            label = '$' + str(labels[i]) + '$'

            ax.scatter(od_lons[i], od_lats[i], s=150, c=color,
                       marker=markers[i], zorder=3)
            ax.scatter(od_lons[i], od_lats[i], s=100, c='k',
                       marker=label, zorder=4)
    # If G not provided or nor OSMNX graph, only plot tour P.
    else:
        labels = helper.node_labels(P.graph['requests'])
        nx.draw_networkx(P, with_labels=True, labels=labels)
