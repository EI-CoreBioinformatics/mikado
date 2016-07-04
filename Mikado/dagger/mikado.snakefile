import sys
import os,sys                                                              
import glob
import itertools
import subprocess
from os import listdir
from os.path import isfile, join
import gzip
from snakemake import logger as snake_logger

CFG=workflow.overwrite_configfile


REF = config["reference"]["fasta"]
# TODO: this is hack that should be solved more neatly
if "out_dir" in config:
    OUT_DIR = config["out_dir"]
else:
    OUT_DIR = config["prepare"]["files"]["output_dir"]
if "threads" in config:
    THREADS = int(config["threads"])
else:
    THREADS = 1


# Directory shortcuts
OUT_DIR_FULL = os.path.abspath(OUT_DIR)
MIKADO_DIR = OUT_DIR + "/5-mikado"
MIKADO_DIR_FULL = os.path.abspath(MIKADO_DIR)
BLAST_DIR = MIKADO_DIR + "/blast"
BLAST_DIR_FULL = os.path.abspath(BLAST_DIR)
TDC_DIR = MIKADO_DIR + "/transdecoder"
TDC_DIR_FULL = os.path.abspath(TDC_DIR)

CWD = os.getcwd()


BLAST_DB_LUT = {}
BLASTX_PROTEIN_DB_PATHS = config["blastx"]["prot_db"]  # config["serialise"]["files"]["blast_targets"]
BLAST_DB_OUT_LIST = list()
BLAST_DB_IN_LIST = list()
for bdbp in BLASTX_PROTEIN_DB_PATHS:
	bdb, ext = os.path.splitext(os.path.basename(bdbp))
	BLAST_DB_LUT[bdb] = bdbp
	BLAST_DB_OUT_LIST.append(BLAST_DIR+"/mikado-" + bdb + "-proteins.xml.gz")
	BLAST_DB_IN_LIST.append(bdb)

BLASTX_PROTEIN_DB_LIST = ','.join(BLAST_DB_OUT_LIST)
BLASTX_TARGET_LIST = ','.join(BLASTX_PROTEIN_DB_PATHS)
BLASTX_MAX_TARGET_SEQS = config["blastx"]["max_target_seqs"]
BLASTX_EVALUE = config["blastx"]["evalue"]

#########################
# Rules

rule all:
	input: 
		mikado=MIKADO_DIR+"/mikado.loci.gff3"

rule clean:
	shell: "rm -rf {OUT_DIR}"

rule mikado_prepare:
	input: 
		ref=REF,
		cfg=CFG
	output:
		gtf=MIKADO_DIR+"/mikado_prepared.gtf",
		fa=MIKADO_DIR+"/mikado_prepared.fasta"
	params:
		load=config["load"]["mikado"]
	threads: THREADS
	message: "Preparing transcripts using mikado"
	shell: "{params.load} && mikado prepare --start-method=spawn --procs={threads} --fasta={input.ref} --json-conf={input.cfg} -od {MIKADO_DIR} 2>&1"

rule make_blast:
	input: fa="{bdb}.fasta"
	output: BLAST_DIR+"/index/{bdb}-proteins.pog"
	params:
		db=BLAST_DIR+"/index/{bdb}-proteins",
		load=config["load"]["blast"]
	log: BLAST_DIR+"/blast.index.log"
	message: "Making BLAST protein database for: {input.fa}"
	shell: "{params.load} && makeblastdb -in {input.fa} -out {params.db} -dbtype prot -parse_seqids > {log} 2>&1"

rule blastx:
	input: 
		db=rules.make_blast.output,
		tr=rules.mikado_prepare.output.fa
	output: BLAST_DIR+"/mikado-{bdb}-proteins.xml.gz"
	params: 
		db=BLAST_DIR + "/index/{bdb}-proteins",
		load=config["load"]["blast"]
	threads: THREADS
	message: "Running BLASTX for mikado transcripts against: {input.tr}"
	run:
	    command="{load} && blastx -num_threads {threads} -query {query} -outfmt 5 -db {db} -evalue {BLASTX_EVALUE} -max_target_seqs {BLASTX_MAX_TARGET_SEQS}".format(load=params.load, threads=threads, query=input.tr, db=params.db, BLASTX_MAX_TARGET_SEQS=BLASTX_MAX_TARGET_SEQS, BLASTX_EVALUE=BLASTX_EVALUE)
	    outfile=gzip.open("{}".format(output), "w")
	    error = open(BLAST_DIR+"/blast-{bdb}.err".format(bdb=bdb), "w")
	    snake_logger.info(command)
	    outfile.writelines(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=error).stdout)
	    outfile.close()
	    error.close()

rule blast_all:
	input: expand(BLAST_DIR + "/mikado-{bdb}-proteins.xml.gz", bdb=BLAST_DB_IN_LIST)
	output: BLAST_DIR + "/blastx.all.done"
	shell: "touch {output}"


rule transdecoder:
	input: rules.mikado_prepare.output.fa
	output: TDC_DIR+"/transcripts.fasta.transdecoder.bed"
	params: outdir=TDC_DIR_FULL,
		tr_in=MIKADO_DIR_FULL+"/mikado_prepared.fasta",
		lolog=TDC_DIR_FULL+"/transdecoder.longorf.log",
		plog=TDC_DIR_FULL+"/transdecoder.predict.log",
		tr="transcripts.fasta",
		load=config["load"]["transdecoder"]
		# ss="-S" if MIKADO_STRAND else ""
	threads: THREADS
	message: "Running transdecoder on Mikado prepared transcripts"
	shell: "{params.load} && cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && TransDecoder.LongOrfs -t {params.tr} > {params.lolog} 2>&1 && TransDecoder.Predict -t {params.tr} --cpu {threads} > {params.plog} 2>&1"


rule genome_index:
	input: os.path.abspath(REF)
	output: MIKADO_DIR+"/"+os.path.basename(REF)+".fai"
	params: load=config["load"]["samtools"],
		fa=MIKADO_DIR+"/"+os.path.basename(REF)
	threads: 1
	message: "Using samtools to index genome"
	shell: "ln -sf {input} {params.fa} && touch -h {params.fa} && {params.load} && samtools faidx {params.fa}"

rule mikado_serialise:
	input: 
		cfg=CFG,
		blast=rules.blast_all.output,
		orfs=rules.transdecoder.output,
		fai=rules.genome_index.output,
		transcripts=rules.mikado_prepare.output.fa
	output: db=MIKADO_DIR+"/mikado.db"
	log: MIKADO_DIR+"/mikado_serialise.err"
	params:
		blast="--xml=" + BLASTX_PROTEIN_DB_LIST if len(BLASTX_PROTEIN_DB_LIST) > 0 else "",
		load=config["load"]["mikado"],
		blast_target="--blast_targets=" + BLASTX_TARGET_LIST if len(BLASTX_TARGET_LIST) > 0 else ""
	threads: 1
	message: "Running Mikado serialise to move numerous data sources into a single database"
	shell: "{params.load} && mikado serialise {params.blast} {params.blast_target} --start-method=spawn --transcripts={input.transcripts} --genome_fai={input.fai} --json-conf={input.cfg} --force --orfs {input.orfs} -od {MIKADO_DIR} --max-objects=200000 > {log} 2>&1"

rule mikado_pick:
	input:
		cfg=CFG,
		gtf=rules.mikado_prepare.output.gtf,
		db=rules.mikado_serialise.output
	output:
		loci=MIKADO_DIR+"/mikado.loci.gff3"
	log: MIKADO_DIR+"/mikado_pick.err"
	params:
		load=config["load"]["mikado"]
	threads: THREADS
	message: "Running mikado picking stage"
	shell: "{params.load} && mikado pick --procs={threads} --start-method=spawn --json-conf={input.cfg} -od {MIKADO_DIR} -lv INFO {input.gtf} -db {input.db} > {log} 2>&1"

rule complete:
  input: rules.mikado_pick.output.loci
