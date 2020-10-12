import sqlite3
from ...utilities.file_type import filetype
from ...exceptions import CorruptIndex
import collections
import msgpack
from ...transcripts import Gene
import os
import rapidjson as json
import sys
import tempfile
import shutil
from .prepare_reference import prepare_reference

# Hack to give the rapidjson library this exception class
# This becomes necessary when we happen to have a corrupted index
if not hasattr(json, "decoder"):

    class Decoder:
        class JSONDecodeError(TypeError):
            pass
    json.decoder = Decoder


def load_index(args, queue_logger):

    """
    Function to load the genes and positions from the indexed GFF.
    :param args:
    :param queue_logger:
    :return: genes, positions
    :rtype: ((None|collections.defaultdict),(None|collections.defaultdict))
    """

    # genes, positions = None, None

    # New: now we are going to use SQLite for a faster experience
    if filetype("{0}.midx".format(args.reference.name)) == b"application/gzip":
        queue_logger.warning("Old index format detected. Starting to generate a new one.")
        raise CorruptIndex("Invalid index file")

    try:
        conn = sqlite3.connect("{0}.midx".format(args.reference.name))
        cursor = conn.cursor()
        tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
        if sorted(tables) != sorted([("positions",), ("genes",)]):
            raise CorruptIndex("Invalid database file")
        # Integrity check
        res = cursor.execute("PRAGMA integrity_check;").fetchone()
        if res[0] != "ok":
            raise CorruptIndex("Corrupt database, integrity value: {}".format(res[0]))

    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid database file")

    positions = dict()
    try:
        for counter, obj in enumerate(cursor.execute("SELECT * from positions")):
            chrom, start, end, gid = obj
            if chrom not in positions:
                positions[chrom] = collections.defaultdict(list)
            positions[chrom][(start, end)].append(gid)
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid index file. Rebuilding.")

    genes = dict()
    for gid, obj in cursor.execute("SELECT * from genes"):
        try:
            gene = Gene(None, logger=queue_logger)
            gene.load_dict(msgpack.loads(obj, raw=False),
                           exclude_utr=args.exclude_utr,
                           protein_coding=args.protein_coding,
                           trust_orf=True)
            if len(gene.transcripts) > 0:
                genes[gid] = gene
            else:
                queue_logger.warning("No transcripts for %s", gid)
        except (EOFError, json.decoder.JSONDecodeError) as exc:
            queue_logger.exception(exc)
            raise CorruptIndex("Invalid index file")
        except (TypeError, ValueError) as exc:
            queue_logger.exception(exc)
            raise CorruptIndex("Corrupted index file; deleting and rebuilding.")

    return genes, positions


def check_index(reference, queue_logger):

    if reference.endswith("midx"):
        reference = reference
    else:
        reference = "{}.midx".format(reference)

    if filetype(reference) == b"application/gzip":
        queue_logger.warning("Old index format detected. Starting to generate a new one.")
        raise CorruptIndex("Invalid index file")

    try:
        conn = sqlite3.connect(reference)
        cursor = conn.cursor()
        tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
        if sorted(tables) != sorted([("positions",), ("genes",)]):
            raise CorruptIndex("Invalid database file")
        gid, obj = cursor.execute("SELECT * from genes").fetchone()
        try:
            obj = msgpack.loads(obj, raw=False)
        except TypeError:
            try:
                obj = json.loads(obj)
            except (ValueError, TypeError, json.decoder.JSONDecodeError):
                raise CorruptIndex("Corrupt index")
            raise CorruptIndex("Old index, deleting and rebuilding")

        gene = Gene(None)
        try:
            gene.load_dict(obj)
        except:
            raise CorruptIndex("Invalid value for genes, indicating a corrupt index. Deleting and rebuilding.")

    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid database file")


def create_index(reference, queue_logger, index_name, ref_gff=False,
                 exclude_utr=False, protein_coding=False):

    """Method to create the simple indexed database for features."""

    temp_db = tempfile.mktemp(suffix=".db")
    queue_logger.info("Starting to create an index for %s", reference.name)
    if os.path.exists("{0}.midx".format(reference.name)):
        queue_logger.warning("Removing the old index")
        try:
            os.remove("{0}.midx".format(reference.name))
        except (OSError, PermissionError) as exc:
            queue_logger.critical(exc)
            queue_logger.critical(
                "I cannot delete the old index, due to permission errors. Please investigate and relaunch.")
            sys.exit(1)

    conn = sqlite3.connect(temp_db)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE positions (chrom text, start integer, end integer, gid text)")
    genes, positions = prepare_reference(reference, queue_logger, ref_gff=ref_gff,
                                         exclude_utr=exclude_utr, protein_coding=protein_coding)

    gid_vals = []
    for chrom in positions:
        for key in positions[chrom]:
            start, end = key
            for gid in positions[chrom][key]:
                gid_vals.append((chrom, start, end, gid))
    cursor.executemany("INSERT INTO positions VALUES (?, ?, ?, ?)",
                       gid_vals)
    cursor.execute("CREATE INDEX pos_idx ON positions (chrom, start, end)")
    cursor.execute("CREATE TABLE genes (gid text, json blob)")

    gobjs = []
    for gid, gobj in genes.items():
        gobjs.append((gid, msgpack.dumps(gobj.as_dict())))
    cursor.executemany("INSERT INTO genes VALUES (?, ?)", gobjs)
    cursor.execute("CREATE INDEX gid_idx on genes(gid)")
    cursor.close()
    conn.commit()
    conn.close()

    shutil.copy(temp_db, index_name)
    os.remove(temp_db)
    queue_logger.info("Finished to create an index for %s in %s", reference.name, index_name)
    return
