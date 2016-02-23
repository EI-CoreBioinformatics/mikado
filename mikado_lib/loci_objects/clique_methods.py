#!/usr/bin/env python3

"""
Module that implements the Reid/Daid/Hurley algorithm for community finding.
"""

import networkx
from ..utilities.log_utils import create_null_logger
from collections import defaultdict

__all__ = ["reid_daid_hurley"]


def reid_daid_hurley(graph, k, cliques=None, logger=None):

    """
    Implementation of the Reid-Daid-Hurley algorithm for clique percolation
    published in http://arxiv.org/pdf/1205.0038.pdf

    :param graph:
    :param k:
    :param cliques:
    :param logger: optional logger for the function
    :return:
    """

    if k < 2:
        raise networkx.NetworkXError("k=%d, k must be greater than 1." % k)
    if cliques is None:
        cliques = [frozenset(x) for x in networkx.find_cliques_recursive(graph)]

    if logger is None:
        logger = create_null_logger("null")

    nodes_to_clique_dict = defaultdict(set)
    # Create the dictionary that links each node to its clique
    logger.debug("Creating the node dictionary")
    cliques = [_ for _ in cliques if len(_) >= k]
    for clique in cliques:
        for node in clique:
            nodes_to_clique_dict[node].add(clique)

    if len(nodes_to_clique_dict) > 100 or len(cliques) > 500:
        logger.debug("Complex locus at %s, with %d nodes and %d cliques with length >= %d",
                     logger.name, len(nodes_to_clique_dict), len(cliques), k)

    current_component = 0

    logger.debug("Starting to explore the clique graph")
    cliques_to_components_dict = dict()
    counter = 0
    for clique in cliques:
        # visited = set()
        counter += 1
        logger.debug("Exploring clique %d out of %d", counter, len(cliques))
        if clique not in cliques_to_components_dict:
            current_component += 1
            cliques_to_components_dict[clique] = current_component
            frontier = set()
            frontier.add(clique)
            cycle = 0
            while len(frontier) > 0:
                current_clique = frontier.pop()
                # if current_clique in visited_cliques:
                #     continue
                cycle += 1
                logger.debug("Cycle %d for clique %d with %d nodes",
                             cycle,
                             counter,
                             len(current_clique))

                for neighbour in _get_unvisited_neighbours(current_clique, nodes_to_clique_dict):
                    if len(frozenset.intersection(current_clique, neighbour)) >= (k-1):
                        cliques_to_components_dict[neighbour] = current_component
                        frontier.add(neighbour)
                        for node in neighbour:
                            nodes_to_clique_dict[node].remove(neighbour)

                logger.debug("Found %d neighbours of clique %d in cycle %d",
                             len(frontier), counter, cycle)

    logger.debug("Finished exploring the clique graph")
    communities = dict()
    for clique in cliques_to_components_dict:
        if cliques_to_components_dict[clique] not in communities:
            communities[cliques_to_components_dict[clique]] = set()
        communities[cliques_to_components_dict[clique]].update(set(clique))

    logger.debug("Reporting the results")
    return [frozenset(x) for x in communities.values()]


def _get_unvisited_neighbours(current_clique, nodes_to_clique_dict):

    """

    :param current_clique:
    :param nodes_to_clique_dict:
    :return:
    """

    neighbours = set()
    for node in current_clique:
        neighbours.update(nodes_to_clique_dict[node])
    return frozenset(neighbours - {current_clique})
