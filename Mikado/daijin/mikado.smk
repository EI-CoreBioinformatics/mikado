import os
from shutil import which
import pkg_resources
from Bio.Data import CodonTable
from Mikado.serializers.blast_serializer.tabular_utils import blast_keys
import functools
import subprocess
import re


diamond_pat = re.compile("^diamond version (\S*)[$|\s]*")

@functools.lru_cache(maxsize=4, typed=True)
def diamond_to_correct(command):
    """This will always return False until https://github.com/bbuchfink/diamond/issues/334 is fixed."""
    if workflow.use_conda is True:
        return False
    cmd = "{} diamond --version && set -u".format(command)
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    version = None
    for line in output.split("\n"):
        m = diamond_pat.search(line)
        if m:
            version = m.groups()[0]
            break
    if version is None:
        return False
    else:
        major, minor, micro = [int(_) for _ in version.split(".")]
        if major > 0 or minor > 9 or micro > 100:
            return True
        else:
            return False


CodonTable.ambiguous_dna_by_id[0] = CodonTable.ambiguous_dna_by_id[1]


try:
    CFG=workflow.overwrite_configfiles[0]
except AttributeError:
    CFG=workflow.overwrite_configfile

envdir = pkg_resources.resource_filename("Mikado.daijin", "envs")


REF = config["reference"]["genome"]
# TODO: this is hack that should be solved more neatly
if "out_dir" in config:
    OUT_DIR = config["out_dir"]
else:
    OUT_DIR = config["prepare"]["files"]["output_dir"]
THREADS = int(config["threads"])

if "mikado" in config:
    MIKADO_MODES=config["mikado"]["modes"]
else:
    MIKADO_MODES = ["split", "nosplit", "permissive", "lenient", "stringent"]

# Directory shortcuts
OUT_DIR_FULL = os.path.abspath(OUT_DIR)
MIKADO_DIR = os.path.join(OUT_DIR, "5-mikado")
MIKADO_DIR_FULL = os.path.abspath(MIKADO_DIR)
BLAST_DIR = os.path.join(MIKADO_DIR, "blast")
BLAST_DIR_FULL = os.path.abspath(BLAST_DIR)
TDC_DIR = os.path.join(MIKADO_DIR, "transdecoder")
TDC_DIR_FULL = os.path.abspath(TDC_DIR)
PROD_DIR = os.path.join(MIKADO_DIR, "prodigal")
PROD_DIR_FULL = os.path.abspath(PROD_DIR)

CWD = os.getcwd()

BLASTX_TARGET = config["blastx"]["prot_db"]
BLASTX_MAX_TARGET_SEQS = config["blastx"].get("max_target_seqs", 10)
BLASTX_EVALUE = config["blastx"]["evalue"]
BLASTX_CHUNKS = max(int(config["blastx"]["chunks"]), THREADS)

ASM_COLLECT = which("asm_collect.py")

CHUNK_ARRAY = []
if len(BLASTX_TARGET) > 0:
    for a in range(1,BLASTX_CHUNKS+1):
        val=str(a).zfill(3)
        CHUNK_ARRAY.append(val)


def loadPre(config, program):
    if workflow.use_conda is True:
            return ""
    cc = config.get("load", dict()).get(program, "")
    if not cc:
        return ""
    else:
        return "set +u && {} &&".format(cc)


td_codes = {
    "Acetabularia": ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5'],
    "Candida": ['Alternative Yeast Nuclear'],
    "Ciliate": ["Condylostoma Nuclear"],
    "Dasycladacean": ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5'],
    "Euplotid": ['Euplotid Nuclear', 'SGC9'],
    "Hexamita": ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5'],
    "Mesodinium": ["Mesodinium Nuclear"],
    "Mitochondrial-Ascidian": ['Ascidian Mitochondrial'],
    "Mitochondrial-Chlorophycean": ['Chlorophycean Mitochondrial'],
    "Mitochondrial-Echinoderm": ['Echinoderm Mitochondrial', 'Flatworm Mitochondrial', 'SGC8'],
    "Mitochondrial-Flatworm": ['Echinoderm Mitochondrial', 'Flatworm Mitochondrial', 'SGC8'],
    "Mitochondrial-Invertebrates": ['Invertebrate Mitochondrial', 'SGC4'],
    "Mitochondrial-Protozoan": ['Mold Mitochondrial', 'Protozoan Mitochondrial',
                                'Coelenterate Mitochondrial', 'Mycoplasma', 'Spiroplasma', 'SGC3'],
    "Mitochondrial-Pterobranchia": ['Pterobranchia Mitochondrial'],
    "Mitochondrial-Scenedesmus_obliquus": ['Scenedesmus obliquus Mitochondrial'],
    "Mitochondrial-Thraustochytrium": ['Thraustochytrium Mitochondrial'],
    "Mitochondrial-Trematode": ['Trematode Mitochondrial'],
    "Mitochondrial-Vertebrates": ['Vertebrate Mitochondrial', "SCG1"],
    "Mitochondrial-Yeast": ['Yeast Mitochondrial', 'SGC2'],
    "Pachysolen_tannophilus": ["Pachysolen tannophilus Nuclear"],
    "Peritrich": ["Peritrich Nuclear"],
    "SR1_Gracilibacteria": ['Candidate Division SR1', 'Gracilibacteria'],
    "Tetrahymena": ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5'],
    "Universal": ['Standard', 'SGC0']
}


td_codes_inverted = dict()
for key, items in td_codes.items():
    for item in items:
        td_codes_inverted[item] = key


def get_codon_table(return_id=True):

    table = config["serialise"]["codon_table"]
    if table == 0:
        table = 1
    if isinstance(table, int):
        table = CodonTable.ambiguous_dna_by_id[table]
    elif isinstance(table, (str, bytes)):
        if isinstance(table, bytes):
            table = table.decode()
        table = CodonTable.ambiguous_dna_by_name[table]

    if return_id is True:
        return table.id
    else:
        found = None
        for name in table.names:
            found = td_codes_inverted.get(name, None)
            if found is not None:
                break
        if found is None:
            found = "Universal"
        return found


#########################
# Rules

rule all:
    input:
        os.path.join(MIKADO_DIR, "pick", "comparison.stats")
    output: touch(os.path.join(MIKADO_DIR, "all.done"))

rule clean:
    shell: "rm -rf {OUT_DIR}"

rule mikado_prepare:
    input: 
        ref=REF
    output:
        gtf=os.path.join(MIKADO_DIR, "mikado_prepared.gtf"),
        fa=os.path.join(MIKADO_DIR, "mikado_prepared.fasta")
    params:
        load=loadPre(config, "mikado"),
        cfg=CFG
    log: os.path.join(MIKADO_DIR, "mikado_prepare.log")
    threads: THREADS
    message: "Preparing transcripts using mikado"
    shell: "{params.load} mikado prepare -l {log} --start-method=spawn --fasta={input.ref} --json-conf={params.cfg} -od {MIKADO_DIR} 2>&1"

rule create_blast_database:
    input: fa=BLASTX_TARGET
    output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
    message: "Creating the BLASTX database"
    params:
        fastas=" ".join(BLASTX_TARGET)
    shell: """sanitize_blast_db.py --out {output} {params.fastas}"""

rule split_fa:
    input: tr=rules.mikado_prepare.output.fa
    output: os.path.join(BLAST_DIR, "fastas", "split.done")
    params: 
        outdir=os.path.join(BLAST_DIR, "fastas", "chunk"),
        chunks=BLASTX_CHUNKS,
        load=loadPre(config, "mikado")
    threads: 1
    message: "Splitting fasta: {input.tr}"
    shell: "{params.load} split_fasta.py -m {params.chunks} {input.tr} {params.outdir} && touch {output}"

if config["mikado"]["use_diamond"] is False:
    # Default, use BLAST
    rule make_blast:
        input:
            fa=os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
        output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.pog")
        params:
            db=os.path.join(BLAST_DIR, "index", "blastdb-proteins"),
            load=loadPre(config, "blast")
        log: os.path.join(BLAST_DIR, "blast.index.log")
        message: "Making BLAST protein database for: {input.fa}"
        conda: os.path.join(envdir, "blast.yaml")
        shell: "{params.load} makeblastdb -in {input.fa} -out {params.db} -dbtype prot -parse_seqids > {log} 2>&1"

    rule blastx:
        input:
            db=rules.make_blast.output,
            split=rules.split_fa.output
        output: os.path.join(BLAST_DIR, "tsv" , "chunk-{chunk_id}-proteins.tsv.gz")
        params:
            tr=os.path.join(BLAST_DIR, "fastas", "chunk_{chunk_id}.fasta"),
            db=os.path.join(BLAST_DIR, "index", "blastdb-proteins"),
            load=loadPre(config, "blast"),
            uncompressed=os.path.join(BLAST_DIR, "tsv", "chunk-{chunk_id}-proteins.tsv"),
            blast_keys=" ".join(blast_keys)
        log: os.path.join(BLAST_DIR, "logs", "chunk-{chunk_id}.blastx.log")
        threads: THREADS
        conda: os.path.join(envdir, "blast.yaml")
        message: "Running BLASTX for mikado transcripts against: {params.tr}"
        shell: "{params.load} if [ -s {params.tr} ]; then blastx -num_threads {threads} "\
"-outfmt \"6 {params.blast_keys}\" "\
" -query {params.tr} -db {params.db} -evalue {BLASTX_EVALUE} -max_target_seqs {BLASTX_MAX_TARGET_SEQS} "\
"> {params.uncompressed} 2> {log}; else touch {params.uncompressed}; fi && gzip {params.uncompressed}"

else:
    rule diamond_index:
        input:
            fa=os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
        output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.dmnd")
        params:
            db=os.path.join(BLAST_DIR, "index", "blastdb-proteins"),
            load=loadPre(config, "diamond")
        log: os.path.join(BLAST_DIR, "diamond.index.log")
        message: "Making DIAMOND protein database for: {input.fa}"
        threads: THREADS
        conda: os.path.join(envdir, "diamond.yaml")
        shell: "{params.load} diamond makedb --threads THREADS --in {input.fa} --db {params.db} 2> {log} > {log}"

    rule diamond:
        input:
            db=rules.diamond_index.output,
            split=rules.split_fa.output
        output: os.path.join(BLAST_DIR, "tsv", "chunk-{chunk_id}-proteins.tsv.gz")
        params:
            load=loadPre(config, "diamond"),
            tr=os.path.join(BLAST_DIR, "fastas", "chunk_{chunk_id}.fasta"),
            blast_keys=" ".join(blast_keys),
            matrix=config["serialise"]["substitution_matrix"],
        threads: THREADS
        log: os.path.join(BLAST_DIR, "logs", "chunk-{chunk_id}.blastx.log")
        conda: os.path.join(envdir, "diamond.yaml")
        shell: "{params.load} if [ -s {params.tr} ]; then diamond blastx --threads {threads} "\
"--outfmt 6 {params.blast_keys} "\
"--max-target-seqs {BLASTX_MAX_TARGET_SEQS} --matrix {params.matrix} "\
"--evalue {BLASTX_EVALUE} --db {input.db} --salltitles --query {params.tr} --sensitive "\
" --compress 1 --out {output} 2> {log} > {log}; else touch {output}; fi"

rule blast_all:
    input: expand(os.path.join(BLAST_DIR, "tsv", "chunk-{chunk_id}-proteins.tsv.gz"), chunk_id=CHUNK_ARRAY)
    output: os.path.join(BLAST_DIR, "blastx.all.done")
    shell: "touch {output}"

if config.get("mikado", dict()).get("use_prodigal", False) is False and config.get("orf_calling", dict()).get("execute", True) is True:
    rule transdecoder_lo:
        input: rules.mikado_prepare.output.fa
        output: os.path.join(TDC_DIR, "transcripts.fasta.transdecoder_dir", "longest_orfs.gff3")
        params:
            outdir=TDC_DIR_FULL,
            tr="transcripts.fasta",
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            load=loadPre(config, "transdecoder"),
            minprot=config["orf_calling"]["min_protein_len"],
            table=get_codon_table(return_id=False)
        log: os.path.join(TDC_DIR_FULL, "transdecoder.longorf.log")
        threads: 1
        message: "Running transdecoder longorf on Mikado prepared transcripts: {input}"
        conda: os.path.join(envdir, "transdecoder.yaml")
        shell: "{params.load} cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && \
TransDecoder.LongOrfs -G {params.table} -m {params.minprot} -t {params.tr} > {log} 2>&1"

    rule transdecoder_pred:
        input:
            mikado=rules.mikado_prepare.output.fa,
            trans=rules.transdecoder_lo.output
        output: os.path.join(TDC_DIR, "transcripts.fasta.transdecoder.bed")
        params:
            outdir=TDC_DIR_FULL,
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            lolog=os.path.join(TDC_DIR_FULL, "transdecoder.longorf.log"),
            plog=os.path.join(TDC_DIR_FULL, "transdecoder.predict.log"),
            tr="transcripts.fasta",
            load=loadPre(config, "transdecoder"),
            table=get_codon_table(return_id=False)
            # ss="-S" if MIKADO_STRAND else ""
        log: os.path.join(TDC_DIR_FULL, "transdecoder.predict.log")
        threads: 1
        message: "Running transdecoder predict on Mikado prepared transcripts: {input}"
        conda: os.path.join(envdir, "transdecoder.yaml")
        shell: "{params.load} cd {params.outdir} && TransDecoder.Predict -t {params.tr} -G {params.table} > {log} 2>&1"
    orf_out = rules.transdecoder_pred.output
elif config.get("mikado", dict()).get("use_prodigal", False) is True and config.get("orf_calling", dict()).get("execute", False) is True:
    rule prodigal:
        input: rules.mikado_prepare.output.fa
        output: os.path.join(PROD_DIR, "transcripts.fasta.prodigal.gff3")
        params:
            outdir=PROD_DIR_FULL,
            tr="transcripts.fasta",
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            tr_out="transcripts.fasta.prodigal.gff3",
            load=loadPre(config, "prodigal"),
            minprot=config["orf_calling"]["min_protein_len"],
            table=get_codon_table(return_id=True)
        log: os.path.join(PROD_DIR_FULL, "prodigal.log")
        threads: 1
        message: "Running PRODIGAL on Mikado prepared transcripts: {input}"
        conda: os.path.join(envdir, "prodigal.yaml")
        shell: "{params.load} mkdir -p {params.outdir} && cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && prodigal -f gff -g {params.table} -i {params.tr} -o {params.tr_out} > {log} 2>&1"
    orf_out = rules.prodigal.output

else:
    rule mock_orfs:
        input:
            mikado = rules.mikado_prepare.output.fa,
            trans = rules.transdecoder_lo.output
        output: os.path.join(TDC_DIR, "transcripts.fasta.transdecoder.bed")
        params:
            outdir = TDC_DIR_FULL,
            tr_in = os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            lolog = os.path.join(TDC_DIR_FULL, "transdecoder.longorf.log"),
            plog = os.path.join(TDC_DIR_FULL, "transdecoder.predict.log"),
            tr = "transcripts.fasta",
            load = loadPre(config, "transdecoder")
        log: os.path.join(TDC_DIR_FULL, "transdecoder.predict.log")
        threads: 1
        message: "Running transdecoder predict on Mikado prepared transcripts: {input}"
        shell: "mkdir -p $(dirname {output}) && touch {output}"


rule genome_index:
    input: os.path.abspath(REF)
    output: os.path.join(MIKADO_DIR, os.path.basename(REF)+".fai")
    params:
        load=loadPre(config, "samtools"),
        fa=os.path.join(MIKADO_DIR, os.path.basename(REF))
    threads: 1
    message: "Using samtools to index genome"
    conda: os.path.join(envdir, "samtools.yaml")
    shell: "ln -sf {input} {params.fa} && touch -h {params.fa} && {params.load} samtools faidx {params.fa}"


if config.get("mikado", dict()).get("use_prodigal") is True:
    orfs = "--orfs={}".format(rules.prodigal.output)
elif config.get("transdecoder", dict()).get("execute", True) is True:
    orfs = "--orfs={}".format(rules.transdecoder_pred.output)
else:
    orfs = ""

rule mikado_serialise:
    input: 
        blast=rules.blast_all.output,
        orfs=orf_out,
        fai=rules.genome_index.output,
        transcripts=rules.mikado_prepare.output.fa
    output: db=os.path.join(MIKADO_DIR, "mikado.db")
    log: os.path.join(MIKADO_DIR, "mikado_serialise.log")
    params:
        cfg=CFG,
        blast="--tsv={}".format(os.path.join(BLAST_DIR, "tsv")) if len(BLASTX_TARGET) > 0 else "",
        load=loadPre(config, "mikado"),
        blast_target="--blast-targets={}".format(os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")) if len(BLASTX_TARGET) > 0 else "",
        orfs=orfs,
        no_start_adj="-nsa" if config.get("mikado", dict()).get("use_prodigal", False) is False else ""
    threads: THREADS
    # conda: os.path.join(envdir, "mikado.yaml")
    shell: "{params.load} mikado serialise {params.blast} {params.blast_target} --start-method=spawn \
--transcripts={input.transcripts} --genome_fai={input.fai} --json-conf={params.cfg} {params.no_start_adj} \
--force {params.orfs} -od {MIKADO_DIR} --procs={threads} -l {log}"

rule mikado_pick:
    input:
        gtf=rules.mikado_prepare.output.gtf,
        db=rules.mikado_serialise.output
    output:
        loci=os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.gff3")
    log: os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.pick.log")
    params:
        cfg=CFG,
        load=loadPre(config, "mikado"),
        outdir=os.path.join(MIKADO_DIR, "pick", "{mode}")
    threads: THREADS
    message: "Running mikado picking stage"
    shell: "{params.load} mikado pick --source Mikado_{wildcards.mode} --mode={wildcards.mode} \
--procs={threads} --start-method=spawn --json-conf={params.cfg} -od {params.outdir} -l {log} \
--loci-out mikado-{wildcards.mode}.loci.gff3 -lv INFO -db {input.db} {input.gtf}"

rule mikado_stats:
    input:
        rules.mikado_pick.output.loci
    output:
        stats=os.path.join(MIKADO_DIR,
                           "pick", "{mode}",
                           "mikado-{mode}.loci.stats")
    params: load=loadPre(config, "mikado")
    shell: "{params.load} mikado util stats {input} > {output}"

rule mikado_collect_stats:
    input:
        mikado=expand(os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.stats"), mode=MIKADO_MODES)
    output: os.path.join(MIKADO_DIR, "pick", "comparison.stats")
    params:
        load=loadPre(config, "mikado")
    threads: 1
    message: "Collecting mikado statistics"
    shell: "{params.load} {ASM_COLLECT} {input.mikado} > {output}"

rule complete:
  input: rules.mikado_collect_stats.output
