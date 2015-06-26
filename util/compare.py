import sys,argparse,os
sys.path.append(os.path.dirname(os.path.dirname( os.path.abspath(__file__)  )))
import collections
#import random
import multiprocessing
import logging
from logging import handlers as log_handlers


from shanghai_lib.loci_objects.transcript import transcript
from shanghai_lib.compare.assigner import assigner
from shanghai_lib.compare.reference_gene import gene
from shanghai_lib.parsers.GTF import GTF
from shanghai_lib.parsers.GFF import GFF3
from shanghai_lib.compare.stat_storer import stat_storer 


'''This is still an embryo. Ideally, this program would perform the following functions:

1- Define precision/recall for the annotation
2- Use a "flexible" mode to determine accuracy
3- Detect gene *fusions* as well as gene *splits*  

'''

def main():
    
    def to_gtf(string):
        '''Function to recognize the input file type and create the parser.'''
        
        if string.endswith(".gtf"):
            return GTF(string)
        elif string.endswith('.gff') or string.endswith('.gff3'):
            return GFF3(string)
        else:
            raise ValueError('Unrecognized file format.')
        
    
    parser=argparse.ArgumentParser('Tool to define the spec/sens of predictions vs. references.')
    input_files=parser.add_argument_group('Prediction and annotation files.')
    input_files.add_argument('-r', '--reference', type=to_gtf, help='Reference annotation file.', required=True)
    input_files.add_argument('-p', '--prediction', type=to_gtf, help='Prediction annotation file.', required=True)
    parser.add_argument('--distance', type=int, default=2000, 
                        help='Maximum distance for a transcript to be considered a polymerase run-on. Default: %(default)s')
    parser.add_argument('-pc', '--protein-coding', dest="protein_coding", action="store_true", default=False,
                        help="Flag. If set, only transcripts with a CDS (both in reference and prediction) will be considered.")
#     parser.add_argument("-t", "--threads", default=1, type=int)
    parser.add_argument("-o","--out", default="shangai_compare", type = str,
                        help = "Prefix for the output files. Default: %(default)s" )
    parser.add_argument("-l","--log", default=None, type = str)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)

    
    args=parser.parse_args()

    #Flags for the parsing
    if type(args.reference) is GFF3:
        ref_gff = True
    else:
        ref_gff = False
    
    if os.path.dirname(args.out)!='' and os.path.dirname(args.out)!=os.path.dirname(os.path.abspath(".")):
        dirname = os.path.dirname(args.out)
        if os.path.exists(dirname):
            assert os.path.isdir(dirname)
        else:
            os.makedirs(dirname)
        
    
    context = multiprocessing.get_context() #@UndefinedVariable
    manager = context.Manager()

    logger = logging.getLogger("main")
    formatter = logging.Formatter("{asctime} - {name} - {levelname} - {message}", style="{")
    args.log_queue = manager.Queue()
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate=False
    log_queue_listener.start()

    if args.log is None:
        handler = logging.StreamHandler()
    else:
        if os.path.exists(args.log):
            os.remove(args.log)
        
        handler = logging.FileHandler(args.log, mode='a')
    handler.setFormatter(formatter)
     
    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    
    
    queue_logger = logging.getLogger("main_queue")
    queue_logger.setLevel(logging.INFO)
    main_queue_handler = log_handlers.QueueHandler(args.log_queue)
    queue_logger.addHandler(main_queue_handler)
    
    queue_logger.propagate=False
    queue_logger.info("Start")
    args.commandline = " ".join(sys.argv)
    queue_logger.info("Command line: {0}".format(args.commandline))
    
#     refmap_queue = manager.Queue(-1)
    
#     pp=threading.Thread(target=printer, args=(args,), name="printing_thread", daemon=True)
#     pp.start()
#     args.refmap_queue = refmap_queue

    queue_logger.info("Starting parsing the reference")

    genes = dict()
    positions = collections.defaultdict(dict)

    transcript2gene = dict()
    logger.handlers[0].flush()
    
    for row in args.reference:
        #Assume we are going to use GTF for the moment
        if row.header is True:
            continue
#         logger.debug(str(row))
        if row.is_transcript is True:
            queue_logger.debug("Transcript\n{0}".format(str(row)))
            tr = transcript(row)
            transcript2gene[row.id]=row.gene
            if row.gene not in genes:
                genes[row.gene] = gene(tr, gid=row.gene)
            genes[row.gene].add(tr)
            assert tr.id in genes[row.gene].transcripts
        elif row.is_exon is True:
#             logger.debug(str(row))
#             assert type(row.transcript) is list
#             logger.debug("Exon found: {0}, {1}".format(row.transcript, row.parent))
            if ref_gff is True:
                for tr in row.transcript:
#                     logger.debug(tr)
                    gid = transcript2gene[tr]
                    genes[gid][tr].addExon(row)
            else:
#                 logger.debug(row.transcript)
                try:
                    genes[row.gene][row.transcript].addExon(row)
                except KeyError as exc:
                    assert row.gene in genes
                    queue_logger.exception(exc)
                    queue_logger.exception( "Keys for {0}: {1}".format(row.gene,  genes[row.gene].transcripts.keys() ))
                    raise
        else:
            continue
#     logger.info("Finished parsing the reference")
    genes_to_remove = set()
    for gid in genes:
        genes[gid].finalize()
        if args.protein_coding is True:
            to_remove=[]
            for tid in genes[gid].transcripts:
                if genes[gid].transcripts[tid].combined_cds_length==0:
                    to_remove.append(tid)
                    logger.debug("No CDS for {0}".format(tid))
            if len(to_remove)==len(genes[gid].transcripts):
                genes_to_remove.add(gid)
                logger.debug("Noncoding gene: {0}".format(tid))
                continue
            elif len(to_remove)>0:
                for tid in to_remove: genes[gid].remove(tid) 
        key = (genes[gid].start,genes[gid].end)
        if key not in positions[genes[gid].chrom]: 
            positions[genes[gid].chrom][key]=[]
        positions[genes[gid].chrom][key].append(genes[gid])
    
    for gid in genes_to_remove: del genes[gid]

    

    #Needed for refmap

    queue_logger.info("Finished preparation")


    
    stat_storer_instance = stat_storer(genes, args) #start the class which will manage the statistics
    assigner_instance = assigner(genes, positions, args, stat_storer_instance)


    currentTranscript = None
    for row in args.prediction:
        if row.header is True:
            continue
#         queue_logger.debug("Row:\n{0:>20}".format(str(row)))
        if row.is_transcript is True:
            queue_logger.debug("Transcript row:\n{0}".format(str(row)))
            if currentTranscript is not None:
                try:
                    result=assigner_instance.get_best(currentTranscript)
                    stat_storer_instance.store(currentTranscript, result, genes)
                except Exception as err:
                    queue_logger.exception(err)
                    log_queue_listener.enqueue_sentinel()
                    handler.close()
                    log_queue_listener.stop()
                    args.queue_handler.close()
                    raise
                
            currentTranscript=transcript(row)
        elif row.is_exon is True:
            try:
                currentTranscript.addExon(row)
            except Exception as err:
                queue_logger.exception(err)
                #In case of error, signal the threads to exit
                args.queue.put("EXIT")
                args.queue.all_tasks_done = True
                args.refmap_queue.put("EXIT")
                break
        else:
            continue

    if currentTranscript is not None:
        try:
            result=assigner_instance.get_best(currentTranscript)
            stat_storer_instance.store(currentTranscript, result, genes)
        except Exception as err:
            queue_logger.exception(err)
            log_queue_listener.enqueue_sentinel()
            handler.close()
            log_queue_listener.stop()
            args.queue_handler.close()
            raise
 
    assigner_instance.finish()
#     stat_storer_instance.print_stats(args)
    
    
    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    handler.close()
    log_queue_listener.stop()
    args.queue_handler.close()
    return


if __name__=='__main__': main()