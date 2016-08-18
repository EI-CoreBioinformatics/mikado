import sys
import os,sys                                                              
import glob
import itertools
import subprocess
from math import log, ceil
from os import listdir
from os.path import isfile, join
from snakemake import logger as snake_logger
from shutil import which

CFG=workflow.overwrite_configfile


REF = config["reference"]["genome"]
# TODO: this is hack that should be solved more neatly
if "out_dir" in config:
    OUT_DIR = config["out_dir"]
else:
    OUT_DIR = config["prepare"]["files"]["output_dir"]
if "threads" in config:
    THREADS = int(config["threads"])
else:
    THREADS = 1

if "mikado" in config:
    MIKADO_MODES=config["mikado"]["modes"]
else:
    MIKADO_MODES = ["split", "nosplit", "permissive", "lenient", "stringent"]

# Directory shortcuts
OUT_DIR_FULL = os.path.abspath(OUT_DIR)
MIKADO_DIR = OUT_DIR + "/5-mikado"
MIKADO_DIR_FULL = os.path.abspath(MIKADO_DIR)
BLAST_DIR = MIKADO_DIR + "/blast"
BLAST_DIR_FULL = os.path.abspath(BLAST_DIR)
TDC_DIR = MIKADO_DIR + "/transdecoder"
TDC_DIR_FULL = os.path.abspath(TDC_DIR)

CWD = os.getcwd()

BLASTX_TARGET = config["blastx"]["prot_db"]
BLASTX_MAX_TARGET_SEQS = config["blastx"]["max_target_seqs"]
BLASTX_EVALUE = config["blastx"]["evalue"]
BLASTX_CHUNKS = int(config["blastx"]["chunks"])
ASM_COLLECT = which("asm_collect.py")

CHUNK_ARRAY = []
if len(BLASTX_TARGET) > 0:
	for a in range(1,BLASTX_CHUNKS+1):
		val=str(a).zfill(3)
		CHUNK_ARRAY.append(val)


def loadPre(command):
	cc = command.strip()
	if not cc:
		return ""
	else:
		return "set +u && {} &&".format(cc)


#########################
# Rules

rule all:
	input:
		MIKADO_DIR+"/pick/comparison.stats"
	output: touch(os.path.join(MIKADO_DIR, "all.done"))

rule clean:
	shell: "rm -rf {OUT_DIR}"

rule mikado_prepare:
	input: 
		ref=REF
	output:
		gtf=MIKADO_DIR+"/mikado_prepared.gtf",
		fa=MIKADO_DIR+"/mikado_prepared.fasta"
	params:
		load=loadPre(config["load"]["mikado"]),
		cfg=CFG
	threads: THREADS
	message: "Preparing transcripts using mikado"
	shell: "{params.load} mikado prepare --start-method=spawn --procs={threads} --fasta={input.ref} --json-conf={params.cfg} -od {MIKADO_DIR} 2>&1"

rule create_blast_database:
    input: fa=BLASTX_TARGET
    output: BLAST_DIR+"/index/blastdb-proteins.fa"
    message: "Creating the BLASTX database"
    params:
         fastas=" ".join(BLASTX_TARGET)
    shell: """sanitize_blast_db.py --out {output} {params.fastas}"""

rule make_blast:
	input: fa=BLAST_DIR+"/index/blastdb-proteins.fa"
	output: BLAST_DIR+"/index/blastdb-proteins.pog"
	params:
		db=BLAST_DIR+"/index/blastdb-proteins",
		load=loadPre(config["load"]["blast"])
	log: BLAST_DIR+"/blast.index.log"
	message: "Making BLAST protein database for: {input.fa}"
	shell: "{params.load} makeblastdb -in {input.fa} -out {params.db} -dbtype prot -parse_seqids > {log} 2>&1"

rule split_fa:
	input: tr=rules.mikado_prepare.output.fa
	output: BLAST_DIR+"/fastas/split.done"
	params: 
		outdir=BLAST_DIR+"/fastas/chunk",
		chunks=config["blastx"]["chunks"],
		load=loadPre(config["load"]["mikado"])
	threads: 1
	message: "Splitting fasta: {input.tr}"
	shell: "{params.load} split_fasta.py -m {params.chunks} {input.tr} {params.outdir} && touch {output}"

rule blastx:
	input: 
		db=rules.make_blast.output,
		split=rules.split_fa.output
	output: BLAST_DIR+"/xmls/chunk-{chunk_id}-proteins.xml.gz"
	params: 
		tr=BLAST_DIR+"/fastas/chunk_{chunk_id}.fasta",
		db=BLAST_DIR + "/index/blastdb-proteins",
		load=loadPre(config["load"]["blast"]),
		uncompressed=BLAST_DIR+"/xmls/chunk-{chunk_id}-proteins.xml",
	log: BLAST_DIR + "/logs/chunk-{chunk_id}.blastx.log"
	threads: THREADS
	message: "Running BLASTX for mikado transcripts against: {params.tr}"
	shell: "{params.load} if [ -s {params.tr} ]; then blastx -num_threads {threads} -outfmt 5 -query {params.tr} -db {params.db} -evalue {BLASTX_EVALUE} -max_target_seqs {BLASTX_MAX_TARGET_SEQS} > {params.uncompressed} 2> {log}; else touch {params.uncompressed}; fi && gzip {params.uncompressed}"


rule blast_all:
	input: expand(BLAST_DIR + "/xmls/chunk-{chunk_id}-proteins.xml.gz", chunk_id=CHUNK_ARRAY)
	output: BLAST_DIR + "/blastx.all.done"
	shell: "touch {output}"


rule transdecoder_lo:
	input: rules.mikado_prepare.output.fa
	output: TDC_DIR+"/transcripts.fasta.transdecoder_dir/longest_orfs.gff3"
	params: outdir=TDC_DIR_FULL,
		tr="transcripts.fasta",
		tr_in=MIKADO_DIR_FULL+"/mikado_prepared.fasta",
		load=loadPre(config["load"]["transdecoder"]),
		minprot=config["transdecoder"]["min_protein_len"]
	log: TDC_DIR_FULL+"/transdecoder.longorf.log",
		# ss="-S" if MIKADO_STRAND else ""
	threads: 1
	message: "Running transdecoder longorf on Mikado prepared transcripts: {input}"
	shell: "{params.load} cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && TransDecoder.LongOrfs -m {params.minprot} -t {params.tr} > {log} 2>&1"

rule transdecoder_pred:
	input: 
		mikado=rules.mikado_prepare.output.fa,
		trans=rules.transdecoder_lo.output		
	output: TDC_DIR+"/transcripts.fasta.transdecoder.bed"
	params: outdir=TDC_DIR_FULL,
		tr_in=MIKADO_DIR_FULL+"/mikado_prepared.fasta",
		lolog=TDC_DIR_FULL+"/transdecoder.longorf.log",
		plog=TDC_DIR_FULL+"/transdecoder.predict.log",
		tr="transcripts.fasta",
		load=loadPre(config["load"]["transdecoder"])
		# ss="-S" if MIKADO_STRAND else ""
	log: TDC_DIR_FULL+"/transdecoder.predict.log"
	threads: THREADS
	message: "Running transdecoder predict on Mikado prepared transcripts: {input}"
	shell: "{params.load} cd {params.outdir} && TransDecoder.Predict -t {params.tr} --cpu {threads} > {log} 2>&1"

rule genome_index:
	input: os.path.abspath(REF)
	output: MIKADO_DIR+"/"+os.path.basename(REF)+".fai"
	params: load=loadPre(config["load"]["samtools"]),
		fa=MIKADO_DIR+"/"+os.path.basename(REF)
	threads: 1
	message: "Using samtools to index genome"
	shell: "ln -sf {input} {params.fa} && touch -h {params.fa} && {params.load} samtools faidx {params.fa}"

rule mikado_serialise:
	input: 
		blast=rules.blast_all.output,
		orfs=rules.transdecoder_pred.output,
		fai=rules.genome_index.output,
		transcripts=rules.mikado_prepare.output.fa
	output: db=MIKADO_DIR+"/mikado.db"
	log: MIKADO_DIR+"/mikado_serialise.err"
	params:
	    cfg=CFG,
		blast="--xml=" + BLAST_DIR+"/xmls" if len(BLASTX_TARGET) > 0 else "",
		load=loadPre(config["load"]["mikado"]),
		blast_target="--blast_targets=" + BLAST_DIR+"/index/blastdb-proteins.fa" if len(BLASTX_TARGET) > 0 else ""
	threads: THREADS
	message: "Running Mikado serialise to move numerous data sources into a single database"
	shell: "{params.load} mikado serialise {params.blast} {params.blast_target} --start-method=spawn --transcripts={input.transcripts} --genome_fai={input.fai} --json-conf={params.cfg} --force --orfs {input.orfs} -od {MIKADO_DIR} --procs={threads} > {log} 2>&1"

rule mikado_pick:
	input:
		gtf=rules.mikado_prepare.output.gtf,
		db=rules.mikado_serialise.output
	output:
		loci=os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.gff3")
	log: os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.pick.err")
	params:
		cfg=CFG,
		load=loadPre(config["load"]["mikado"]),
		outdir=os.path.join(MIKADO_DIR, "pick", "{mode}")
	threads: THREADS
	message: "Running mikado picking stage"
	shell: "{params.load} mikado pick --source Mikado_{wildcards.mode} --mode={wildcards.mode} --procs={threads} --start-method=spawn --json-conf={params.cfg} -od {params.outdir} --loci_out mikado-{wildcards.mode}.loci.gff3 -lv INFO {input.gtf} -db {input.db} > {log} 2>&1"

rule mikado_stats:
    input:
        rules.mikado_pick.output.loci
    output:
        stats=os.path.join(MIKADO_DIR,
                           "pick", "{mode}",
                           "mikado-{mode}.loci.stats")
    params: load=loadPre(config["load"]["mikado"])
    shell: "{params.load} mikado util stats {input} > {output}"

rule mikado_collect_stats:
	input:
	    mikado=expand(MIKADO_DIR+"/pick/{mode}/mikado-{mode}.loci.stats", mode=MIKADO_MODES)
	output: MIKADO_DIR+"/pick/comparison.stats"
	params:
	    load=loadPre(config["load"]["mikado"])
	threads: 1
	message: "Collecting mikado statistics"
	shell: "{params.load} {ASM_COLLECT} {input.mikado} > {output}"

rule complete:
  input: rules.mikado_collect_stats.output
