#!/usr/bin/env python3

"""
Module that implements the Reid/Daid/Hurley algorithm for community finding.
"""

import networkx
from ..utilities.log_utils import create_null_logger
from collections import defaultdict
from itertools import chain, combinations


def find_communities(graph: networkx.Graph, logger=None) -> set:
    """

    :param graph: a Graph instance from networkx
    :type graph: networkx.Graph

    :param logger: optional logger. A default null one will be created if none is provided.
    :type logger: (None|logging.Logger)

    This function is a wrapper around the networkX methods to find communities inside a graph.
    The method takes as input a precomputed graph and returns a set of the available communities.
    """
    if logger is None:
        logger = create_null_logger()

    logger.debug("Creating the communities for %s", logger.name)
    communities = [frozenset(comm) for comm in networkx.connected_components(graph)]

    logger.debug("Communities for %s:\n\t\t%s", logger.name, "\n\t\t".join(
        [str(_) for _ in communities]))
    # result = [frozenset(x) for x in communities.values()]
    for element in set.difference(set(graph.nodes()), set(chain(*communities[:]))):
        communities.append(frozenset([element]))
    return set(communities)


def define_graph(objects: dict, inters, **kwargs) -> networkx.Graph:
    """
    :param objects: a dictionary of objects to be grouped into a graph
    :type objects: dict

    :param inters: the intersecting function to be used to define the graph
    :type inters: callable

    :param kwargs: optional arguments to be passed to the inters function
    :type kwargs: dict

    This function will compute the graph which will later be used by find_communities.
    The method takes as mandatory inputs the following:
        - "objects" a dictionary of objects that form the graph
        - "inters" a function/method that determines whether two objects are connected or not.

    It will then return a graph.
    The method accepts also kwargs that can be passed to the inters function.
    WARNING: the kwargs option is really stupid and does not check
    for correctness of the arguments!
    """

    graph = networkx.Graph()

    # As we are using intern for transcripts, this should prevent
    # memory usage to increase too much
    graph.add_nodes_from(objects.keys())

    for obj, other_obj in combinations(objects.keys(), 2):
        if obj == other_obj:
            continue
        elif inters(objects[obj], objects[other_obj], **kwargs):
            # Connections are not directional
            graph.add_edge(*tuple(sorted([obj, other_obj])))

    return graph


def find_cliques(graph: networkx.Graph, logger=None) -> (networkx.Graph, list):
    """

    :param graph: graph to which it is necessary to call the cliques for.

    :param logger: optional logger for the function

    Wrapper for the BronKerbosch algorithm, which returns the maximal cliques in the graph.
    It is the new interface for the BronKerbosch function, which is not called directly
    from outside this class any longer.
    The "inters" keyword provides the function used to determine
    whether two vertices are connected or not in the graph.
    """

    if logger is None:
        logger = create_null_logger()

    logger.debug("Creating cliques for %s", logger.name)
    cliques = [frozenset(x) for x in networkx.find_cliques_recursive(graph)]
    logger.debug("Created %d cliques for %s", len(cliques), logger.name)

    return cliques
