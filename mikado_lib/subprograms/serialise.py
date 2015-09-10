import argparse
import glob
import sqlalchemy
from mikado_lib import json_utils
from mikado_lib.serializers import orf, blast_utils, junction, dbutils
import os
from Bio import SeqIO

__author__ = 'Luca Venturini'


def to_seqio(string):
    """
    Convert a string to a SeqIO index.

    :param string
    :type string: str
    """
    assert os.path.exists(string) and os.path.isfile(string) and os.stat(string).st_size > 0
    return SeqIO.index(string, "fasta")


def serialise(args):

    if args.db is not None:
        args.json_conf["db"] = args.db
        args.json_conf["dbtype"] = "sqlite"

    if args.force is True:
        engine = dbutils.connect(args.json_conf)
        meta = sqlalchemy.MetaData(bind=engine)
        meta.reflect(engine)
        for tab in reversed(meta.sorted_tables):
            tab.drop()
        dbutils.dbBase.metadata.create_all(engine)

    if args.orfs is not None:
        for orf_file in args.orfs.split(","):
            serializer = orf.OrfSerializer(orf_file, fasta_index=args.transcript_fasta, maxobjects=args.max_objects,
                                           json_conf=args.json_conf)
            serializer.serialize()

    if args.junctions is not None:
        for junction_file in args.junctions.split(","):
            serializer = junction.JunctionSerializer(junction_file, args.db,
                                                     fai=args.genome_fai, json_conf=args.json_conf,
                                                     maxobjects=args.max_objects)
            serializer()

    if args.xml is not None:
        for xml in args.xml.split(","):
            if os.path.isdir(xml):
                candidates = [os.path.join(xml, x) for x in
                              filter(
                              lambda x: x.endswith(".xml") or x.endswith(".xml.gz") or x.endswith(".asn.gz"),
                              os.listdir(xml)
                              )
                              ]
            elif "*" in xml:
                candidates = glob.glob(xml)
            else:
                candidates = [xml]

            list(
                map(
                    lambda xml_candidate: blast_utils.XmlSerializer(
                        xml_candidate,
                        discard_definition=args.discard_definition,
                        max_target_seqs=args.max_target_seqs,
                        maxobjects=args.max_objects,
                        target_seqs=args.target_seqs,
                        query_seqs=args.transcript_fasta,
                        json_conf=args.json_conf
                    )(),
                    candidates
                )
            )


def serialise_parser():
    """
    Parser function for the serialisation step.
    :return: argparse.Namespace
    """

    parser = argparse.ArgumentParser("Serialisation utility of the Mikado suite.")
    orfs = parser.add_argument_group()
    orfs.add_argument("--orfs", type=str, default=None,
                      help="ORF BED file(s), separated by commas")
    orfs.add_argument("--transcript_fasta", default=None,
                      help="""Transcript FASTA file(s) used for ORF calling and BLAST queries, separated by commas.
                      If multiple files are given, they must be in the same order of the ORF files.
                      E.g. valid command lines are:

                      --transcript_fasta all_transcript1.fasta --orfs all_orfs.bed
                      --transcript_fasta transcript1.fasta,transcript2.fasta --orfs orfs1.bed,orf2.bed
                      --transcript_fasta all_transcript.fasta --orfs orfs1.bed,orf2.bed

                      These are invalid instead:

                      # Inverted order
                      --transcript_fasta transcript1.fasta,transcript2.fasta --orfs orfs2.bed,orf1.bed
                      #Two transcript files, one ORF file
                      --transcript_fasta transcript1.fasta,transcript2.fasta --orfs all_orfs.bed

                      """)

    blast = parser.add_argument_group()
    blast.add_argument("--max_target_seqs", type=int, default=float("Inf"),
                       help="Maximum number of target sequences.")
    blast.add_argument("--target_seqs", default=None, type=to_seqio, help="Target sequences")
    blast.add_argument("--discard-definition", action="store_true", default=False,
                       help="""Flag. If set, the sequences IDs instead of their definition
                       will be used for serialisation.""")
    blast.add_argument("--xml", type=str, help="XML file(s) to parse, separated by a comma.")

    junctions = parser.add_argument_group()
    junctions.add_argument("--genome_fai", default=None)
    junctions.add_argument("--junctions", type=str)
    generic = parser.add_argument_group()
    generic.add_argument("-mo", "--max-objects", dest="max_objects", type=int, default=10 ** 5,
                         help="Maximum number of objects to cache in memory.")
    generic.add_argument("-f", "--force", action="store_true", default=False,
                         help="Flag. If set, an existing databse will be dropped/deleted before serialisation.")
    generic.add_argument("--json-conf", default=None, dest="json_conf", type=json_utils.to_json,
                         required=True)
    generic.add_argument("db", type=str, default=None,
                         nargs='?', help="Optional output database. Default: derived from json_conf")
    parser.set_defaults(func=serialise)
    return parser
