import tempfile
import multiprocessing as mp
from ...utilities.namespace import Namespace
import os
import functools
from ..assignment.assigner import Assigner
from ..assignment.distributed import Assigners, FinalAssigner
from .transmission import get_best_result, transmit_transcript, _transmit_transcript
from ..reference_preparation.gene_dict import GeneDict
from ..resultstorer import ResultStorer
from ...transcripts import Transcript, Gene
from ...parsers.GFF import GFF3
from ...parsers.GTF import GTF
from ...parsers.bed12 import Bed12Parser
from ...parsers.bam_parser import BamParser
from .parse_bam_prediction import parse_prediction_bam
from .parse_gff3_prediction import parse_prediction_gff3
from .parse_gtf_prediction import parse_prediction_gtf
from .parse_bed12 import parse_prediction_bed12
from ...parsers import to_gff
import gzip
import csv
import itertools


def parse_prediction(args, index, queue_logger):

    """
    This function performs the real comparison between the reference and the prediction.
     It needs the following inputs:
    :param args: the Namespace with the necessary parameters
    :param genes: Dictionary with the reference genes, of the form
    dict[chrom][(start,end)] = [gene object]
    :param positions: Dictionary with the positions of the reference genes, of the form
    dict[chrom][IntervalTree]
    :param queue_logger: Logger
    :return:
    """

    # start the class which will manage the statistics
    if hasattr(args, "self") and args.self is True:
        args.prediction = to_gff(args.reference.name)
    __found_with_orf = set()
    queue = mp.JoinableQueue(1000)
    returnqueue = mp.JoinableQueue(1000)

    queue_logger.info("Starting to parse the prediction")
    if args.processes > 1:
        log_queue = args.log_queue
        dargs = dict()
        doself = False
        for key, item in args.__dict__.items():
            if key in ("log_queue", "queue_handler"):
                continue
            elif key == "self":
                doself = item
            else:
                dargs[key] = item
        nargs = Namespace(default=False, **dargs)
        nargs.self = doself
        procs = [Assigners(index, nargs, queue, returnqueue, log_queue, counter)
                 for counter in range(1, args.processes)]
        [proc.start() for proc in procs]
        final_proc = FinalAssigner(index, nargs, returnqueue, log_queue=log_queue, nprocs=len(procs))
        final_proc.start()
        transmitter = functools.partial(transmit_transcript, queue=queue)
        assigner_instance = None
    else:
        procs = []
        assigner_instance = Assigner(index, args, printout_tmap=True, )
        transmitter = functools.partial(get_best_result, assigner_instance=assigner_instance)

    transmit_wrapper = functools.partial(_transmit_transcript,
                                         transmitter=transmitter,
                                         queue_logger=queue_logger)

    constructor = functools.partial(Transcript,
                                    logger=queue_logger, trust_orf=True, accept_undefined_multi=True)

    if args.prediction.__annot_type__ == BamParser.__annot_type__:
        annotator = parse_prediction_bam
    elif args.prediction.__annot_type__ == Bed12Parser.__annot_type__:
        annotator = parse_prediction_bed12
    elif args.prediction.__annot_type__ == GFF3.__annot_type__:
        annotator = parse_prediction_gff3
    elif args.prediction.__annot_type__ == GTF.__annot_type__:
        annotator = parse_prediction_gtf
    else:
        raise ValueError("Unsupported input file format")

    done, lastdone = annotator(args, queue_logger, transmit_wrapper, constructor)

    queue_logger.info("Finished parsing, %s transcripts in total", done)
    if assigner_instance is None:
        queue.join()
        queue.put("EXIT")
        [proc.join() for proc in procs]
        returnqueue.join()
        returnqueue.put("EXIT")
        final_proc.join()
    else:
        assert not procs, procs
        assert isinstance(assigner_instance, Assigner)
        assigner_instance.finish()


def parse_self(args, gdict: GeneDict, queue_logger):

    """
    This function is called when we desire to compare a reference against itself.
    :param args:
    :param gdict: gene MIDX database
    :param queue_logger:
    :return:
    """

    queue_logger.info("Performing a self-comparison.")

    if args.gzip is False:
        tmap_out = open("{0}.tmap".format(args.out), 'wt')
    else:
        tmap_out = gzip.open("{0}.tmap.gz".format(args.out), 'wt')
    tmap_rower = csv.DictWriter(tmap_out, ResultStorer.__slots__, delimiter="\t")
    tmap_rower.writeheader()

    for gid, gene in gdict.items():
        assert isinstance(gene, Gene)
        if len(gene.transcripts) == 1:
            continue

        combinations = itertools.combinations(gene.transcripts.keys(), 2)
        for combination in combinations:
            result, _ = Assigner.compare(gene.transcripts[combination[0]],
                                         gene.transcripts[combination[1]],
                                         fuzzy_match=args.fuzzymatch)
            tmap_rower.writerow(result.as_dict())

    queue_logger.info("Finished.")
