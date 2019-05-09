import unittest

import networkx

from ..transcripts.clique_methods import find_cliques, find_communities
from ..transcripts.clique_methods import reid_daid_hurley


class TestCliques(unittest.TestCase):

    """"""

    def setUp(self):

        r"""
        We are going to create a graph composed of three disjoint subgraphs:


        0            4                8
          \        /   \            / | \
            2 -- 3  --  5          7 --- 9        12
          /       \   /             \ | /  \
        1           6                10 --  11


        :return:
        """

        self.graph = networkx.Graph()
        self.graph.add_nodes_from([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

        # First graph
        self.graph.add_edges_from([(0,2), (1,2), (2,3), (3,4),
                                   (3, 5), (3, 6), (4, 5), (5, 6)])
        # Second graph
        self.graph.add_edges_from([(7, 8), (7, 9), (7, 10),
                                   (8, 9), (8, 10), (9, 10),
                                   (9, 11), (10, 11)])
        # Third graph is a singleton

        self.maxDiff = None

        self.correct_cliques = {
            frozenset([0, 2]),
            frozenset([1, 2]),
            frozenset([2, 3]),
            frozenset([3, 4, 5]),
            frozenset([3, 5, 6]),
            frozenset([7, 8, 9, 10]),
            frozenset([9, 10, 11]),
            frozenset([12])}

        self.correct_communities = {
            frozenset([0, 1, 2, 3, 4, 5, 6]),
            frozenset([7, 8, 9, 10, 11]),
            frozenset([12])}

    def test_find_cliques(self):

        cliques = find_cliques(self.graph)
        self.assertEqual(set(cliques),
                         self.correct_cliques,
                         cliques)

    def test_reid(self):

        comms = reid_daid_hurley(self.graph, 2)
        self.assertEqual(comms,
                         self.correct_communities,
                         comms)

    def test_comms(self):
        self.maxDiff = None
        comms = find_communities(self.graph)
        self.assertEqual(comms,
                         self.correct_communities,
                         (comms, self.correct_communities))

    def test_error_reid(self):

        with self.assertRaises(networkx.NetworkXError):
            _ = reid_daid_hurley(self.graph, 1)

if __name__ == '__main__':
    unittest.main()
