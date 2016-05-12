from Mikado.loci.transcriptchecker import TranscriptChecker
import unittest


class TChekerTester(unittest.TestCase):

    def test_translation_table(self):

        self.assertEqual(TranscriptChecker.get_translation_table(),
                         {65: 84, 67: 71, 71: 67, 84: 65})

    def test_rev_complement(self):

        string = "AGTCGTGCAGNGTCGAAGTGCAACAGTGC"

        self.assertEqual(TranscriptChecker.rev_complement(string),
                         "GCACTGTTGCACTTCGACNCTGCACGACT")

    def test_init(self):

        """"""

        seq = "AGTCAGTGCCAGTCGATGTNNNNNNNNAGTGTCAGTAGTCAGTGCAGTNNNNNNNNAGTGTTGACA"




if __name__ == '__main__':
    unittest.main()