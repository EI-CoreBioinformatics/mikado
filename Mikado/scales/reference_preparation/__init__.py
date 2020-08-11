from ...parsers.GFF import GFF3
from .indexing import create_index, check_index
from ...exceptions import CorruptIndex
import os
import tempfile


def prepare_index(args, queue_logger):
    index_name = os.path.abspath("{0}.midx".format(args.reference.name))
    ref_gff = isinstance(args.reference, GFF3)
    if args.index is True:
        create_index(args.reference, queue_logger=queue_logger, index_name=index_name,
                     ref_gff=ref_gff,
                     exclude_utr=False, protein_coding=False)
        assert os.path.exists(index_name)
        return index_name

    if os.path.exists(index_name) and os.stat(args.reference.name).st_mtime >= os.stat(index_name).st_mtime:
        queue_logger.warning("Reference index obsolete, deleting and rebuilding.")
        try:
            os.remove("{0}.midx".format(args.reference.name))
        except (OSError, PermissionError) as exc:
            queue_logger.error(
                "I cannot delete the old index due to permission errors. I will create a temporary one instead.")
            __index = tempfile.NamedTemporaryFile(suffix=".midx")
            index_name = __index.name

        ref_gff = isinstance(args.reference, GFF3)
        create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                     protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
    elif os.path.exists(index_name):
        # queue_logger.info("Starting loading the indexed reference")
        queue_logger.info("Index found")
        try:
            check_index(args.reference.name, queue_logger)
            queue_logger.info("Index valid, proceeding.")
        except CorruptIndex as exc:
            queue_logger.warning(exc)
            queue_logger.warning("Reference index corrupt, deleting and rebuilding.")
            try:
                os.remove("{0}.midx".format(args.reference.name))
            except (OSError, PermissionError) as exc:
                queue_logger.error(
                    "I cannot delete the old index due to permission errors. I will create a temporary one instead.")
                __index = tempfile.NamedTemporaryFile(suffix=".midx")
                index_name = __index.name
            create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                         protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
    else:
        if args.no_save_index is True:
            __index = tempfile.NamedTemporaryFile(suffix=".midx", delete=False)
            index_name = __index.name
        create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                     protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
        assert os.path.exists(index_name)

    assert os.path.exists(index_name), (args.reference, index_name)
    return index_name
