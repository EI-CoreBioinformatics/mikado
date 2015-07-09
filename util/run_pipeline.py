#!/usr/bin/env python3

#import sys
#from logging import Logger
import argparse
import sys,os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import shanghai_lib.loci_objects

def main():
        
    parser=argparse.ArgumentParser("Launcher of the Shanghai pipeline.")
    parser.add_argument("-p", "--procs", type=int, default=None, help="Number of processors to use. Default: look in the configuration file (1 if undefined)")
    parser.add_argument("--json_conf", type=argparse.FileType("r"), required=True, help="JSON configuration file for scoring transcripts.")
    parser.add_argument("--subloci_out", type=str, default=None)
    parser.add_argument("--monoloci_out", type=str, default=None)
    parser.add_argument("--loci_out", type=str, default=None,
                        help="This output file is mandatory. If it is not specified in the configuration file, it must be provided here.")
    parser.add_argument("--no_cds", action="store_true", default=None,
                        help="Flag. If set, not CDS information will be printed out in the GFF output files.")
    parser.add_argument('--source', type=str, default=None,
                        help='Source field to use for the output files.')
    parser.add_argument('--purge', action='store_true', default=None,
                        help='Flag. If set, the pipeline will suppress any loci whose transcripts do not pass the requirements set in the JSON file.'
                        )
    log_options=parser.add_argument_group("Log options")
    log_options.add_argument("-l", "--log", default=None, help="File to write the log to. Default: decided by the configuration file.")
    verbosity=log_options.add_mutually_exclusive_group()
    verbosity.add_argument("--verbose", action="store_true", default=False, help="Flag. If set, the debug mode will be activated.")
    verbosity.add_argument("--noverbose", action="store_true", default=False, help="Flag. If set, the debug mode will be activated.")
    log_options.add_argument("-lv", "--log-level", dest="log_level", choices=["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"], default=None,
                             help="Loggin level. Default: retrieved by the configuration file.")
    parser.add_argument("gff", type=argparse.FileType("r"), nargs="?", default=None)
    
    args=parser.parse_args()

    args.json_conf.close()
    args.json_conf = shanghai_lib.json_utils.to_json(args.json_conf.name)
    
    if args.procs is not None:
        args.json_conf["run_options"]["threads"] = args.procs
        
    if args.log == "stderr":
        args.json_conf["log_settings"]['log']=None
    elif args.log is not None:
        args.json_conf["log_settings"]['log']=args.log
        
    if args.log_level is not None:
        args.json_conf["log_settings"]['log_level']=args.log_level
    elif args.verbose is True:
        args.json_conf["log_settings"]['log_level']="DEBUG"
    elif args.noverbose is True:
        args.json_conf["log_settings"]['log_level']="ERROR"
        
    if args.monoloci_out is not None:
        args.json_conf["monoloci_out"] = args.monoloci_out
    if args.subloci_out is not None:
        args.json_conf["subloci_out"] = args.subloci_out
    if args.loci_out is not None:
        args.json_conf["loci_out"] = args.loci_out
    if args.no_cds is not None:
        args.json_conf["run_options"]["exclude_cds"]=True
    if args.source is not None:
        args.json_conf["source"]=args.source
    if args.purge is not None:
        args.json_conf["run_options"]["purge"]=True
        
    if args.gff is not None:
        args.gff.close()
        args.gff=args.gff.name
        args.json_conf["input"]=args.gff
    
    creator = shanghai_lib.loci_objects.Creator.Creator(args.json_conf, commandline = " ".join(sys.argv))
    creator() #Run
       
if __name__=="__main__": main()
