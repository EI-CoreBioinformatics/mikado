import argparse
import sqlite3

import pandas as pd
import pyfaidx

from Mikado.parsers import parser_factory

__doc__ = """Script to generate a test case for Mikado from a run"""

from Mikado.utilities import to_region


def generate_in_argument(id_list):
    return ','.join(f"'{name}'" for name in id_list)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-p", "--prepared", type=parser_factory, help="Mikado prepare output file")
    parser.add_argument("-db", "--db", help="Mikado populated database (after serialise)")
    parser.add_argument("-ref", "--reference", type=pyfaidx.Fasta, help="Fasta reference file")
    parser.add_argument("-r", "--region", type=str, help="Region of interest to extract")
    parser.add_argument("-w", "--window", type=int, default=0,
                        help="Window up/down stream of the gene start/end to extract")
    parser.add_argument("-o", "--output_prefix", type=str, help="Output files prefix", default="test")
    args = parser.parse_args()


    id_list = []
    in_region = False
    target_chrom, target_start, target_end = to_region(args.region)

    # Extract region from genome
    min_seq_pos = max(1, target_start - 1 - args.window)
    max_seq_pos = min(len(args.reference[target_chrom]), target_end + args.window)
    seq = args.reference.get_seq(target_chrom,
                                 min_seq_pos + 1,
                                 max_seq_pos).seq
    print(f">{args.region}\n{seq}", file=open(args.output_prefix+'.fasta', 'w'))

    # Extract transcripts from prepare
    with open(args.output_prefix+".gtf", 'w') as gtf_output:
        for row in args.prepared:
            if row.header is True:
                continue
            else:
                if target_chrom == row.chrom and row.start >= target_start and row.end <= target_end:
                    in_region = True
                    if row.is_transcript:
                        # Generate a list of transcripts and get their query.ID from db
                        id_list.append(row.id)
                    row.start -= min_seq_pos
                    row.end -= min_seq_pos
                    row.chrom = args.region
                    print(row, file=gtf_output)
                elif in_region and (row.start >= target_end or target_chrom != row.chrom):
                    break

    input_db = sqlite3.connect(database=args.db)
    output_db = sqlite3.connect(database=args.output_prefix+'.db')

    chroms = pd.read_sql_query(f"SELECT * FROM chrom", input_db)
    chroms.loc[(chroms.name == target_chrom), 'name'] = args.region
    chroms.to_sql("chrom", output_db, if_exists='replace')

    # Extract db junctions from region remembering to shift by the same amount the genome is 'shifted'
    junctions = pd.read_sql_query(f"SELECT junctions.* FROM junctions "
                                  f"JOIN chrom ON junctions.chrom_id == chrom.chrom_id "
                                  f"WHERE "
                                  f"junctions.start >= {target_start} "
                                  f"AND junctions.end <= {target_end} AND "
                                  f"chrom.name == \"{target_chrom}\"", input_db)

    # Transform junctions to location in the extracted region including the window
    junctions['start'] -= min_seq_pos
    junctions['end'] -= min_seq_pos
    junctions['junction_start'] -= min_seq_pos
    junctions['junction_end'] -= min_seq_pos

    # Store junctions
    junctions.to_sql('junctions', output_db, if_exists='replace')

    # Extract db hsp,hit from transcript list using query.ID list
    query = pd.read_sql_query(f"SELECT * FROM query WHERE query_name IN ({generate_in_argument(id_list)})", input_db)
    query.to_sql('query', output_db, if_exists='replace')

    hsp = pd.read_sql_query(f"SELECT hsp.* FROM hsp WHERE query_id IN ({generate_in_argument(query['query_id'])})",
                            input_db)
    hsp.to_sql('hsp', output_db, if_exists='replace')

    hit = pd.read_sql_query(f"SELECT hit.* FROM hit WHERE hit.query_id IN ({generate_in_argument(query['query_id'])})",
                            input_db)
    hit.to_sql('hit', output_db, if_exists='replace')

    # Extract orf from db
    orf = pd.read_sql_query(f"SELECT * from orf WHERE query_id IN ({generate_in_argument(query['query_id'])})",
                            input_db)
    orf.to_sql('orf', output_db, if_exists='replace')

    target = pd.read_sql_query(f"SELECT * from target WHERE target_id IN ({generate_in_argument(hit['target_id'])})",
                               input_db)
    target.to_sql('target', output_db, if_exists='replace')


if __name__ == "__main__":
    main()
