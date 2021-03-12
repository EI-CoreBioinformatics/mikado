import argparse
import sys
from ..utilities.log_utils import create_default_logger
from ..configuration.daijin_configurator import create_daijin_config
from .mikado import mikado_pipeline
from .assemble import assemble_transcripts_pipeline
from . import create_parser, create_config_parser


def main(call_args=None):
    """
    Main call function.
    :param call_args: Arguments to use to launch the pipeline. If unspecified, the default behaviour
    (using CL arguments) will be adopted.
    :return:
    """

    if call_args is None:
        call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        """A Directed Acyclic pipeline for gene model reconstruction from RNA seq data.
        Basically, a pipeline for driving Mikado. It will first align RNAseq reads against
        a genome using multiple tools, then creates transcript assemblies using multiple tools,
        and find junctions in the alignments using Portcullis.
        This input is then passed into Mikado.

        WARNING: the "assemble" part of this pipeline will be soon DEPRECATED. 
        """)

    subparsers = parser.add_subparsers(
        title="Pipelines",
        help="""These are the pipelines that can be executed via daijin.""")

    subparsers.add_parser("configure",
                          help="Creates the configuration files for Daijin execution.")
    subparsers.choices["configure"] = create_config_parser()
    subparsers.choices["configure"].prog = "daijin configure"
    subparsers.choices["configure"].set_defaults(func=create_daijin_config)

    subparsers.add_parser("assemble",
                          description="Creates transcript assemblies from RNAseq data.",
                          help="""A pipeline that generates a variety of transcript assemblies
                          using various aligners and assemblers, as well a producing
                          a configuration file suitable for driving Mikado.
                          WARNING: this part of the Daijin pipeline will be DEPRECATED in future releases
                          as a new, more complete annotation pipeline is currently in development.""")
    subparsers.choices["assemble"] = create_parser()
    subparsers.choices["assemble"].add_argument(
        "config",
        help="Configuration file to use for running the transcript assembly pipeline.")
    subparsers.choices["assemble"].prog = "daijin assemble"
    subparsers.choices["assemble"].set_defaults(func=assemble_transcripts_pipeline)

    subparsers.add_parser("mikado",
                          description="Run full mikado pipeline",
                          help="""Using a supplied configuration file that describes
                          all input assemblies to use, it runs the Mikado pipeline,
                          including prepare, BLAST, transdecoder, serialise and pick.""")
    subparsers.choices["mikado"] = create_parser()
    subparsers.choices["mikado"].add_argument(
        "config",
        help="Configuration file to use for running the Mikado step of the pipeline.")
    subparsers.choices["mikado"].prog = "daijin mikado"
    subparsers.choices["mikado"].set_defaults(func=mikado_pipeline)

    try:
        args = parser.parse_args(call_args)
        if hasattr(args, "func"):
            args.func(args)
        else:
            parser.print_help()
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except BrokenPipeError:
        pass
    except Exception as exc:
        logger = create_default_logger("main")
        logger.error("Daijin crashed, cause:")
        logger.exception(exc)
        sys.exit(1)


if __name__ == '__main__':
    # pylint: disable=redefined-builtin
    # noinspection PyShadowingBuiltins
    __spec__ = "Daijin"
    # pylint: enable=redefined-builtin
    main()
