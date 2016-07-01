import sys
import os,sys                                                              
import glob
import itertools
import subprocess
import yaml
from os import listdir
from os.path import isfile, join
 

CFG=workflow.overwrite_configfile

# Get shortcuts from configuration file
R1 = config["r1"]
R2 = config["r2"]
SAMPLES = config["samples"]
STRANDEDNESS = config["strandedness"]

REF = config["ref"]
NAME = config["name"]
OUT_DIR = config["out_dir"]
MIN_INTRON = config["min_intron"]
MAX_INTRON = config["max_intron"]
THREADS = int(config["threads"])
TGG_MAX_MEM = config["tgg_max_mem"]
MIKADO_MODE = config["mikado_mode"]

# List of alignment and assembly methods to test
ALIGNMENT_METHODS = config["align_methods"]
ASSEMBLY_METHODS = config["asm_methods"]


# Directory shortcuts
OUT_DIR_FULL = os.path.abspath(OUT_DIR)
READS_DIR = OUT_DIR + "/1-reads"
READS_DIR_FULL = os.path.abspath(READS_DIR)
ALIGN_DIR = OUT_DIR + "/2-alignments"
ALIGN_DIR_FULL = os.path.abspath(ALIGN_DIR)
ASM_DIR = OUT_DIR + "/3-assemblies"
ASM_DIR_FULL = os.path.abspath(ASM_DIR)
PORTCULLIS_DIR = OUT_DIR + "/4-portcullis"
PORTCULLIS_DIR_FULL = os.path.abspath(PORTCULLIS_DIR)

CWD = os.getcwd()
ALIGN_COLLECT = os.path.dirname(os.path.abspath(workflow.snakefile)) + "/align_collect.py"
ASM_COLLECT = os.path.dirname(os.path.abspath(workflow.snakefile)) + "/asm_collect.py"
CLASS = os.path.dirname(os.path.abspath(workflow.snakefile)) + "/class_run.py"

TRINITY_STRAND = "--SS_lib_type=RF" if STRANDEDNESS == "fr-firststrand" else "--SS_lib_type=FR" if STRANDEDNESS == "fr-secondstrand" else ""
PORTCULLIS_STRAND = "firststrand" if STRANDEDNESS == "fr-firststrand" else "secondstrand" if STRANDEDNESS == "fr-secondstrand" else "unstranded"


SAMPLE_MAP = {}
INPUT_1_MAP = {}
INPUT_2_MAP = {}
EXT_MAP = {}
for i in range(len(SAMPLES)):
	SAMPLE_MAP[SAMPLES[i]] = STRANDEDNESS[i]
	name, ext = os.path.splitext(R1[i])
	r1 = READS_DIR+"/"+SAMPLES[i]+".R1.fq"
	r2 = READS_DIR+"/"+SAMPLES[i]+".R2.fq"
	if ext == ".gz" or ext == ".bz2":
		r1 += ext
		r2 += ext
	EXT_MAP[SAMPLES[i]] = ext
	INPUT_1_MAP[SAMPLES[i]] = r1
	INPUT_2_MAP[SAMPLES[i]] = r2


SAMPLE_STR = ",".join(SAMPLES)

RUN_PORTCULLIS = config["portcullis"]["do"]
CANONICAL_JUNCS = config["portcullis"]["canonical_juncs"]

PORTCULLIS_ALIGNMENT_METHODS = ALIGNMENT_METHODS if RUN_PORTCULLIS else []

PORTCULLIS_IN = {}
PORTCULLIS_OUT = []
for a in PORTCULLIS_ALIGNMENT_METHODS :
	s=""
	for samp in SAMPLES :
		s += ALIGN_DIR+"/output/" + a + "-" + samp + ".sorted.bam "
	PORTCULLIS_IN[a] = s.strip()
	PORTCULLIS_OUT.append(PORTCULLIS_DIR+"/output/"+a+".pass.junctions.bed")

PORTCULLIS_OUT_STR_VAL = ",".join(PORTCULLIS_OUT) if RUN_PORTCULLIS else ""
PORTCULLIS_OUT_STR = "--junctions " + PORTCULLIS_OUT_STR_VAL if RUN_PORTCULLIS else ""


aln_abrv = {"tophat":"tph", "star":"sta", "gsnap":"gsp", "hisat":"hst"}
asm_abrv = {"cufflinks":"cuf", "stringtie":"stn", "class":"cls", "trinity":"trn"}
	
GTF_ASSEMBLY_METHODS = []
GFF_ASSEMBLY_METHODS = []
for asm in ASSEMBLY_METHODS:
	if asm == "cufflinks" or asm == "stringtie" or asm == "class":
		GTF_ASSEMBLY_METHODS.append(asm)
	elif asm == "trinity":
		GFF_ASSEMBLY_METHODS.append(asm)

TRANSCRIPT_ARRAY=[]
LABEL_ARRAY=[]
SS_ARRAY=[]
for samp in SAMPLES:
	for aln in ALIGNMENT_METHODS:
		for gtf in GTF_ASSEMBLY_METHODS:
			TRANSCRIPT_ARRAY.append(ASM_DIR + "/output/" + gtf + "-" + aln + "-" + samp + ".gtf")
			LABEL_ARRAY.append(asm_abrv[gtf] + "-" + aln_abrv[aln] + "-" + samp)
			if not SAMPLE_MAP[samp] == "fr-unstranded":
				SS_ARRAY.append(TRANSCRIPT_ARRAY[-1])
		for gff in GFF_ASSEMBLY_METHODS:
			TRANSCRIPT_ARRAY.append(ASM_DIR + "/output/" + gff + "-" + aln + "-" + samp + ".gff")
			LABEL_ARRAY.append(asm_abrv[gff] + "-" + aln_abrv[aln] + "-" + samp)
			if not SAMPLE_MAP[samp] == "fr-unstranded":
				SS_ARRAY.append(TRANSCRIPT_ARRAY[-1])

TRANSCRIPTS_STR = ",".join(TRANSCRIPT_ARRAY)
LABEL_STR = ",".join(LABEL_ARRAY)
SS_STR = ",".join(SS_ARRAY)


#########################
Rules

localrules: mikado_cfg

rule all:
	input:
		mikado=OUT_DIR + "/mikado.yaml",
		align=ALIGN_DIR+"/alignment.stats",
		asm=ASM_DIR+"/assembly.stats"

rule clean:
	shell: "rm -rf {OUT_DIR}"


rule align_tophat_index:
	input: ref=REF
	output: ALIGN_DIR + "/tophat/index/" + NAME + ".4.bt2"
	params: 
		idxdir=ALIGN_DIR + "/tophat/index/" + NAME,
		load=config["load"]["tophat"]
	log: ALIGN_DIR + "/tophat.index.log"
	threads: 1
	message: "Indexing genome with tophat"
	shell: "{params.load} && bowtie2-build {input.ref} {params.idxdir} > {log} 2>&1"


rule align_tophat:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		r2=lambda wildcards: INPUT_2_MAP[wildcards.sample],
		index=rules.align_tophat_index.output
	output:
		link=ALIGN_DIR+"/output/tophat-{sample}.bam"
	params: 
		outdir=ALIGN_DIR+"/tophat/{sample}",
		bam=ALIGN_DIR+"/tophat/{sample}/accepted_hits.bam",
		indexdir=ALIGN_DIR+"/tophat/index/"+NAME,
		load=config["load"]["tophat"],
		extra=config["extra"]["tophat"],
		link_src="../tophat/{sample}/accepted_hits.bam"
	log: ALIGN_DIR + "/tophat-{sample}.log"
	threads: THREADS
	message: "Aligning RNAseq data with tophat: {input.r1}; {input.r2}"
	shell: "{params.load} && tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={STRANDEDNESS} {params.extra} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule align_gsnap_index:
	input: REF
	output: ALIGN_DIR + "/gsnap/index/" + NAME + "/" + NAME + ".sachildguide1024"
	params: load=config["load"]["gmap"]
	log: ALIGN_DIR +"/gsnap.index.log"
	threads: 1
	message: "Indexing genome with gsnap"
	shell: "{params.load} && gmap_build --dir={ALIGN_DIR}/gsnap/index --db={NAME} {input} > {log} 2>&1"


rule align_gsnap:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		r2=lambda wildcards: INPUT_2_MAP[wildcards.sample],
		index=rules.align_gsnap_index.output
	output:
		bam=ALIGN_DIR+"/gsnap/{sample}/gsnap.bam",
		link=ALIGN_DIR+"/output/gsnap-{sample}.bam"
	params: 
		load=config["load"]["gmap"],
		load_sam=config["load"]["samtools"],
		extra=config["extra"]["gsnap"],
		link_src="../gsnap/{sample}/gsnap.bam"
	log: ALIGN_DIR+"/gsnap-{sample}.log"
	threads: THREADS
	message: "Aligning RNAseq with gsnap: {input.r1}; {input.r2}"
	shell: "{params.load} && {params.load_sam} && gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {params.extra} --novelsplicing=1 --localsplicedist={MAX_INTRON} --nthreads={threads} --format=sam --npaths=20 {input.r1} {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"



rule align_star_index:
	input: os.path.abspath(REF)
	output: ALIGN_DIR +"/star/index/SAindex"
	params: 
		indexdir=ALIGN_DIR_FULL+"/star/index",
		load=config["load"]["star"],
		extra=config["extra"]["star_index"]
	log: ALIGN_DIR_FULL+"/star.index.log"
	threads: THREADS
	message: "Indexing genome with star"
	shell: "{params.load} && cd {ALIGN_DIR_FULL}/star && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input} {params.extra} > {log} 2>&1 && cd {CWD}"


rule align_star:
	input:
		r1=lambda wildcards: os.path.abspath(INPUT_1_MAP[wildcards.sample]),
		r2=lambda wildcards: os.path.abspath(INPUT_2_MAP[wildcards.sample]),
		index=rules.align_star_index.output
	output:
		bam=ALIGN_DIR+"/star/{sample}/Aligned.out.bam",
		link=ALIGN_DIR+"/output/star-{sample}.bam"
	params:
		outdir=ALIGN_DIR_FULL+"/star/{sample}",
		indexdir=ALIGN_DIR_FULL+"/star/index",
		load=config["load"]["star"],
		extra=config["extra"]["star"],
		link_src="../star/{sample}/Aligned.out.bam"
	log: ALIGN_DIR_FULL+"/star-{sample}.log"
	threads: int(THREADS)
	message: "Aligning input with star: {input.r1}; {input.r2}"
	run: 
		rfc = ""
		if EXT_MAP[wildcards.sample] == ".gz":
			rfc = "--readFilesCommand=zcat"
		elif EXT_MAP[wildcards.sample] == ".gz":
			rfc = "--readFilesCommand=bzcat"
		shell("{params.load} && cd {params.outdir}; STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} {rfc} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} --alignMatesGapMax 20000 --outFileNamePrefix {params.outdir}/ {params.extra} > {log} 2>&1 && cd {CWD} && ln -sf {params.link_src} {output.link} && touch -h {output.link}")



rule align_hisat_index:
	input: REF
	output: ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
	params: load=config["load"]["hisat"]
	log: ALIGN_DIR+"/hisat.index.log"
	threads: 1
	message: "Indexing genome with hisat"
	shell: "{params.load} && hisat2-build {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"



rule align_hisat:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		r2=lambda wildcards: INPUT_2_MAP[wildcards.sample],
		index=rules.align_hisat_index.output
	output:
		bam=ALIGN_DIR+"/hisat/{sample}/hisat.bam",
		link=ALIGN_DIR+"/output/hisat-{sample}.bam"
	params:
		indexdir=ALIGN_DIR+"/hisat/index/"+NAME,
		load=config["load"]["hisat"],
		load_samtools=config["load"]["samtools"],
		extra=config["extra"]["hisat"],
		link_src="../hisat/{sample}/hisat.bam"
	log: ALIGN_DIR+"/hisat-{sample}.log"
	threads: THREADS
	message: "Aligning input with hisat: {input.r1}; {input.r2}"
    	run:
		strand = "--rna-strandness=RF" if SAMPLE_MAP[wildcards.sample] == "fr-firststrand" else "--rna-strandness=FR" if SAMPLE_MAP[wildcards.sample] == "fr-secondstrand" else ""
		shell("{params.load} && {params.load_samtools} && hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {strand} {params.extra} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}")


rule bam_sort:
	input: ALIGN_DIR+"/output/{align_method}-{sample}.bam"
	output: ALIGN_DIR+"/output/{align_method}-{sample}.sorted.bam"
	params:
		load=config["load"]["samtools"],
	threads: int(THREADS)
	message: "Using samtools to sort {input}"
	shell: "{params.load} && samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.align_method}_{wildcards.sample} -@ {threads} {input}"



rule bam_index:
	input: rules.bam_sort.output
	output: ALIGN_DIR+"/output/{align_method}-{sample}.sorted.bam.bai"
	params: load=config["load"]["samtools"]
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{params.load} && samtools index {input}"

rule bam_stats:
	input:
		bam=rules.bam_sort.output,
		idx=rules.bam_index.output
	output: ALIGN_DIR+"/output/{align_method}-{sample}.sorted.bam.stats"
	params: 
		load=config["load"]["samtools"],
		plot_out=ALIGN_DIR+"/output/plots/{align_method}-{sample}/{align_method}-{sample}"
	threads: 1
	message: "Using samtools to collected stats for: {input}"
	shell: "{params.load} && samtools stats {input.bam} > {output} && plot-bamstats -p {params.plot_out} {output}"


rule align_all:
	input: expand(ALIGN_DIR+"/output/{align_method}-{sample}.sorted.bam.stats", align_method=ALIGNMENT_METHODS, sample=SAMPLES)
	output: ALIGN_DIR+"/output/all.done"
	threads: 1
	shell: "touch {output}"

rule align_collect_stats:
	input: rules.align_all.output
	output: ALIGN_DIR+"/alignment.stats"
	params: stats=ALIGN_DIR+"/output/*.sorted.bam.stats"
	threads: 1
	message: "Collecting alignment stats"
	shell: "{ALIGN_COLLECT} {params.stats} > {output}"
		

rule asm_cufflinks:
	input: 
		bam=rules.bam_sort.output,
		ref=REF
	output: 
		gtf=ASM_DIR+"/output/cufflinks-{align_method}-{sample}.gtf"
	params: 
		outdir=ASM_DIR+"/cufflinks-{align_method}-{sample}",
		gtf=ASM_DIR+"/cufflinks-{align_method}-{sample}/transcripts.gtf",
		link_src="../cufflinks-{align_method}-{sample}/transcripts.gtf",
		load=config["load"]["cufflinks"],
		extra=config["extra"]["cufflinks"]
	log: ASM_DIR+"/cufflinks-{align_method}-{sample}.log"
	threads: int(THREADS)
	message: "Using cufflinks to assemble: {input.bam}"
	shell: "{params.load} && cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --no-update-check {params.extra} {input.bam} > {log} 2>&1 && ln -sf {params.link_src} {output.gtf} && touch -h {output.gtf}"



rule asm_trinitygg:
	input:
		bam=rules.bam_sort.output,
		ref={REF}
	output: ASM_DIR+"/trinity-{align_method}-{sample}/Trinity-GG.fasta"
	params: 
		outdir=ASM_DIR+"/trinity-{align_method}-{sample}",
		load=config["load"]["trinity"],
		extra=config["extra"]["trinity"]
	log: ASM_DIR+"/trinity-{align_method}-{sample}.log"
	threads: int(THREADS)
	message: "Using trinity in genome guided mode to assemble: {input.bam}"
	shell: "{params.load} && Trinity --seqType=fq {TRINITY_STRAND} --output={params.outdir} --genome_guided_bam={input.bam} {params.extra} --full_cleanup --genome_guided_max_intron={MAX_INTRON} --max_memory={TGG_MAX_MEM} --CPU={threads} > {log} 2>&1"


rule asm_map_trinitygg:
	input: 
		transcripts=rules.asm_trinitygg.output,
		index=ASM_DIR +"/gmap_index/"+NAME+"/"+NAME+".sachildguide1024"
	output: 
		gff=ASM_DIR+"/output/trinity-{align_method}-{sample}.gff"
	params: 
		load=config["load"]["gmap"],
		gff=ASM_DIR+"/trinity-{align_method}-{sample}/trinity-{align_method}-{sample}.gff",
		link_src="../trinity-{align_method}-{sample}/trinity-{align_method}-{sample}.gff"
	threads: int(THREADS)
	message: "Mapping trinity transcripts to the genome: {input.transcripts}"
	shell: "{params.load} && gmap --dir={ASM_DIR}/gmap_index --db={NAME} --min-intronlength={MIN_INTRON} --intronlength={MAX_INTRON} --format=3 {input.transcripts} > {params.gff} && ln -sf {params.link_src} {output.gff} && touch -h {output.gff}"


rule asm_gmap_index:
        input: REF
        output: ASM_DIR +"/gmap_index/"+NAME+"/"+NAME+".sachildguide1024"
	params:	load=config["load"]["gmap"]
	threads: 1
        log: ASM_DIR+"/trinity_gmap_index.log"
        message: "Indexing genome with gmap"
        shell: "{params.load} && gmap_build --dir={ASM_DIR}/gmap_index --db={NAME} {input} > {log} 2>&1"


rule asm_stringtie:
	input: 
		bam=rules.bam_sort.output
	output: 
		link=ASM_DIR+"/output/stringtie-{align_method}-{sample}.gtf",
		gtf=ASM_DIR+"/stringtie-{align_method}-{sample}/stringtie-{align_method}-{sample}.gtf"
	params:
		load=config["load"]["stringtie"],
		extra=config["extra"]["stringtie"],
		gtf=ASM_DIR+"/stringtie-{align_method}-{sample}/stringtie-{align_method}-{sample}.gtf",
		link_src="../stringtie-{align_method}-{sample}/stringtie-{align_method}-{sample}.gtf"
	log: ASM_DIR+"/stringtie-{align_method}-{sample}.log"
	threads: int(THREADS)
	message: "Using stringtie to assemble: {input.bam}"
	shell: "{params.load} && stringtie {input.bam} -l Stringtie_{wildcards.align_method}_{wildcards.sample} -f 0.05 -m 200 {params.extra} -o {params.gtf} -p {threads} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule asm_class:
        input: 
                bam=rules.bam_sort.output,
                ref=REF
        output: 
                link=ASM_DIR+"/output/class-{align_method}-{sample}.gtf",
                gtf=ASM_DIR+"/class-{align_method}-{sample}/class-{align_method}-{sample}.gtf"
        params: 
                outdir=ASM_DIR+"/class-{align_method}-{sample}",
                load=config["load"]["class"],
                extra=config["extra"]["class"],
		link_src="../class-{align_method}-{sample}/class-{align_method}-{sample}.gtf"
        log: ASM_DIR+"/class-{align_method}-{sample}.log"
        threads: int(THREADS)
        message: "Using class to assemble: {input.bam}"
        shell: "{params.load} && {CLASS} --clean --force -c \"{params.extra}\" -p {threads} {input.bam} > {output.gtf} 2> {log} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule asm_gff_stats:
	input: ASM_DIR + "/output/{gff_method}-{align_method}-{sample}.gff"
	output: ASM_DIR + "/output/{gff_method}-{align_method}-{sample}.gff.stats"
	params: load=config["load"]["mikado"]
	threads: 1
	message: "Computing assembly stats for: {input}"
	shell: "{params.load} && mikado util stats {input} > {output}"

rule asm_gtf_stats:
	input: ASM_DIR + "/output/{gtf_method}-{align_method}-{sample}.gtf"
	output: ASM_DIR + "/output/{gtf_method}-{align_method}-{sample}.gtf.stats"
	params: load=config["load"]["mikado"]
	threads: 1
	message: "Computing assembly stats for: {input}"
	shell: "{params.load} && mikado util stats {input} > {output}"

rule asm_all:
	input: 	rules.align_all.output,
		expand(ASM_DIR + "/output/{gtf_method}-{align_method}-{sample}.gtf.stats", gtf_method=GTF_ASSEMBLY_METHODS, align_method=ALIGNMENT_METHODS, sample=SAMPLES),
		expand(ASM_DIR + "/output/{gff_method}-{align_method}-{sample}.gff.stats", gff_method=GFF_ASSEMBLY_METHODS, align_method=ALIGNMENT_METHODS, sample=SAMPLES)	
	output:	ASM_DIR+"/output/all.done"
	threads: 1
	shell: "touch {output}"

rule asm_collect_stats:
	input: rules.asm_all.output
	output: ASM_DIR+"/assembly.stats"
	params: asms=ASM_DIR+"/output/*.stats"
	threads: 1
	message: "Collecting assembly statistics"
	shell: "{ASM_COLLECT} {params.asms} > {output}"

rule portcullis_prep:
	input:	ref=REF,
		aln_done=rules.align_all.output
	output: PORTCULLIS_DIR+"/portcullis_{aln_method}/1-prep/portcullis.sorted.alignments.bam.bai"
	params: 
		outdir=PORTCULLIS_DIR+"/portcullis_{aln_method}/1-prep",
		load=config["load"]["portcullis"],
		files=lambda wildcards: PORTCULLIS_IN[wildcards.aln_method]
	log: PORTCULLIS_DIR+"/portcullis_{aln_method}-prep.log"
	threads: int(THREADS)
	message: "Using portcullis to prepare: {input}"
	shell: "{params.load} && portcullis prep -o {params.outdir} --strandedness={PORTCULLIS_STRAND} -t {threads} {input.ref} {params.files} > {log} 2>&1"


rule portcullis_junc:
	input: 
		bai=rules.portcullis_prep.output
	output: PORTCULLIS_DIR+"/portcullis_{aln_method}/2-junc/{aln_method}.junctions.tab"
	params: 
		prepdir=PORTCULLIS_DIR+"/portcullis_{aln_method}/1-prep",
		outdir=PORTCULLIS_DIR+"/portcullis_{aln_method}/2-junc",
		load=config["load"]["portcullis"]
	log: PORTCULLIS_DIR+"/portcullis_{aln_method}-junc.log"
	threads: int(THREADS)
	message: "Using portcullis to analyse potential junctions: {input}"
	shell: "{params.load} && portcullis junc -o {params.outdir}/{wildcards.aln_method} --strandedness={PORTCULLIS_STRAND} -t {threads} {params.prepdir} > {log} 2>&1"

rule portcullis_filter:
	input: rules.portcullis_junc.output
	output:
		link=PORTCULLIS_DIR+"/output/portcullis_{aln_method}.pass.junctions.bed"
	params: 
		outdir=PORTCULLIS_DIR+"/portcullis_{aln_method}/3-filt",
		prepdir=PORTCULLIS_DIR+"/portcullis_{aln_method}/1-prep/",
		load=config["load"]["portcullis"],
		bed=PORTCULLIS_DIR+"/portcullis_{aln_method}/3-filt/{aln_method}.pass.junctions.bed",
		link_src="../portcullis_{aln_method}/3-filt/{aln_method}.pass.junctions.bed"
	log: PORTCULLIS_DIR+"/portcullis_{aln_method}-filter.log"
	threads: 1
	message: "Using portcullis to filter invalid junctions: {input}"
	shell: "{params.load} && portcullis filter -o {params.outdir}/{wildcards.aln_method} --canonical={CANONICAL_JUNCS} --max_length={MAX_INTRON} {params.prepdir} {input} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"



rule portcullis_merge:
	input: expand(PORTCULLIS_DIR + "/output/portcullis_{aln_method}.pass.junctions.bed", aln_method=PORTCULLIS_ALIGNMENT_METHODS)
	output: bed=PORTCULLIS_DIR + "/output/portcullis.merged.bed"
	params: load=config["load"]["portcullis"]
	log: PORTCULLIS_DIR + "/output/portcullis.merged.log"
	threads: 1
	message: "Taking intersection of portcullis results"
	shell: "{params.load} && bed_merge.py --prefix=portcullis_merged --output={output.bed} {input} > {log} || touch {output.bed}"


#### Now run mikado
rule mikado_cfg:
	input:
		asm=rules.asm_all.output,
		portcullis=rules.portcullis_merge.output,
		cfg=CFG
	output:
		mikado=OUT_DIR + "/mikado.yaml"
	params: 
		load=config["load"]["mikado"],
		mikado=OUT_DIR + "/mikado.cfg",
		scoring=config["mikado_scoring"]
	log: OUT_DIR + "/mikado.yaml.log"
	threads: 1
	message: "Creating Mikado configuration file"
	shell: "{params.load} && mikado.py configure --full --mode={MIKADO_MODE} --gff={TRANSCRIPTS_STR} --labels={LABEL_STR} --strand-specific-assemblies={SS_STR} --junctions={input.portcullis} --reference={REF} > {params.mikado} 2> {log} && cat {input.cfg} {params.mikado} > {output} && rm {params.mikado}"

