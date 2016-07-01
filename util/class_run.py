#!/usr/bin/env python3

import subprocess,sys,argparse,os
import shutil
import logging

def main():

    logging.basicConfig(format=" %(asctime)s - %(levelname)s - %(message)s",
                        level=logging.INFO)

    parser=argparse.ArgumentParser("Quick utility to rewrite the wrapper for CLASS.")
    parser.add_argument("--clean", default=False, action="store_true",
                        help="Flag. If set, remove tepmorary files.")
    parser.add_argument("--force", default=False, action="store_true",
                        help="Flag. If set, it forces recalculation of all intermediate files.")
    parser.add_argument("-c","--class_options", type=str, default='',
                        help="Additional options to be passed to CLASS. Default: no additional options.")
    parser.add_argument("-p", "--processors", type=int, default=1,
                        help="Number of processors to use with class.")
    parser.add_argument("--class_help", action="store_true",
                        default=False,
                        help="If called, the wrapper will ask class to display its help and exit.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)                        
    parser.add_argument("bam", type=argparse.FileType('rb'), default=None, nargs="?", help="Input BAM file.")
    parser.add_argument("out", nargs="?", type=argparse.FileType('wt'), default=sys.stdout,
                        help="Optional output file.")
    args=parser.parse_args()

    if args.class_help:
        print("Calling CLASS help..", file=sys.stderr)
        subprocess.call("class", shell=True)
        sys.exit(0)

    if args.bam is None:
        parser.error("The input BAM is required as an argument.")

    args.class_options+=" -p {0} ".format(args.processors) #Add the processors as argument
    logging.info("CLI: {0}".format(
        " ".join(sys.argv)))
    
    args.bam.close() #Quick and dirty check that the file exists.
    args.bam=os.path.abspath(args.bam.name) #Absolute position
    prefix=os.path.splitext(args.bam)[0] #Prefix without the .bam extension (comprehensive of dirs)

    #    if shutil.which("samtools") is None:
    # logging.debug("Loading the SAMTOOLS utility")
    # subprocess.call('source samtools-1.1', shell=True) #Load samtools if it is not present already
    #    else: pass
        
    # #    if shutil.which("junc") is None or shutil.which("class") is None:
    # logging.debug("Loading CLASS")
    # subprocess.call("source class-2.12", shell=True)
    # #   else: pass

    depth_file="{0}.depth".format(prefix)
    if not os.path.exists(depth_file) or args.force:
        if os.path.exists(depth_file):
            logging.warning("Deleting old depth file, because of --force option")
        with open(depth_file, 'wt') as depth_buffer:
            logging.info("Calculating depth with samtools")
            subprocess.call('samtools depth {0}'.format(args.bam), stdout=depth_buffer, shell=True)
    else:
        logging.warning("Depth file already present. Skipping this phase.")
    
    splice_file="{0}.splice".format(prefix)
    if not os.path.exists(splice_file) or args.force:
        if os.path.exists(splice_file):
            logging.warning("Deleting old splice file, because of --force option")
        with open(splice_file, 'wt') as splice_buffer:
            logging.info("Recovering junctions with the internal utility junc")
            subprocess.call("junc {0}".format(args.bam), stdout=splice_buffer, shell=True)
    else:
        logging.warning("Splice file already present. Skipping this phase.")


    logging.info("Launching the main process")
    class_cli="class {0} {1}".format(prefix, args.class_options)
    if args.verbose:
        class_cli+=" --verbose "
    logging.info("CLASS CLI:\n\t{0}".format(class_cli))
    class_sub=subprocess.Popen(class_cli, shell=True,
                    stdout=args.out, stderr=subprocess.PIPE)

    for line in class_sub.stderr:
        line=line.decode().rstrip()
        logging.info("CLASS message:\n{0}".format(line))

    if args.clean:
        logging.info("Removing temporary files..")
        os.remove(depth_file)
        logging.debug("Removed the DEPTH file")
        os.remove(splice_file)
        logging.debug("Removed the SPLICE file")

    logging.info("CLASS finished!")
    return

if __name__=='__main__': main()
