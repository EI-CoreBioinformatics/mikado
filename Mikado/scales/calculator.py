import csv
from ..transcripts.transcriptcomputer import TranscriptComputer
from ..loci import Gene
from ..parsers.GFF import GFF3
from ..utilities.log_utils import create_default_logger
from collections import defaultdict
import sklearn.utils.extmath
import numpy
from collections import Counter
# pylint: disable=E1101
numpy.seterr(all="ignore")  # Suppress warnings
numpy.warnings.filterwarnings("ignore")
# pylint: enable=E1101


def itemize(counter):
    """
    Private static method to convert a counter into a 2-dimensional numpy array.
    :param counter: the counter to transform
    :type counter: dict
    :return: Numpy array of the counter
    :rtype: array
    """

    if len(counter) == 0:
        return numpy.array([[], []])  # Empty 2dimensional array
    else:
        return numpy.array(list(zip(*counter.items())))


def weighted_percentile(a, percentile=numpy.array([75, 25]), weights=None):
    """
    O(nlgn) implementation for weighted_percentile.
    From: http://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
    Kudos to SO user Nayyary http://stackoverflow.com/users/2004093/nayyarv

    :param a: array
    :type a: (dict|list|numpy.array|set|tuple)

    :param percentile: the percentiles to calculate.
    :type percentile: (numpy.array|list|tuple)

    """

    percentile = numpy.array(percentile)/100.0

    if isinstance(a, (dict, Counter)):
        a, weigths = itemize(a)
    else:
        assert isinstance(a, (list, set, tuple, numpy.ndarray)), (a, type(a))
        if not isinstance(a, type(numpy.array)):
            a = numpy.array(a)
        if weights is None:
            weights = numpy.ones(len(a))

    a_indsort = numpy.argsort(a)
    a_sort = a[a_indsort]
    weights_sort = weights[a_indsort]
    ecdf = numpy.cumsum(weights_sort)

    percentile_index_positions = percentile * (weights.sum() - 1) + 1
    # need the 1 offset at the end due to ecdf not starting at 0
    locations = numpy.searchsorted(ecdf, percentile_index_positions)

    out_percentiles = numpy.zeros(len(percentile_index_positions))

    for i, empiricalLocation in enumerate(locations):
        # iterate across the requested percentiles
        if ecdf[empiricalLocation-1] == numpy.floor(percentile_index_positions[i]):
            # i.e. is the percentile in between 2 separate values
            uppWeight = percentile_index_positions[i] - ecdf[empiricalLocation-1]
            lowWeight = 1 - uppWeight

            out_percentiles[i] = a_sort[empiricalLocation-1] * lowWeight + \
                                 a_sort[empiricalLocation] * uppWeight
        else:
            # i.e. the percentile is entirely in one bin
            out_percentiles[i] = a_sort[empiricalLocation]

    return out_percentiles


class Calculator:

    """
    This class has the purpose of parsing a reference file,
    calculating the statistics, and printing them out.
    """

    def __init__(self, parsed_args):

        """Constructor function"""

        self.gff = parsed_args.gff

        self.atype = self.gff.__annot_type__

        self.__logger = create_default_logger("calculator")
        if parsed_args.verbose is True:
            self.__logger.setLevel("DEBUG")
        self.only_coding = parsed_args.only_coding
        self.tab_handle = parsed_args.tab_stats
        if self.tab_handle is not None:
            __fieldnames = ["TID", "GID",
                            "Coordinates",
                            "Strand",
                            "Exon number",
                            "cDNA length",
                            "Intronic length",
                            "CDS length",
                            "# coding exons",
                            "cDNA/CDS ratio",
                            "5'UTR",
                            "5'UTR exons",
                            "3'UTR",
                            "3'UTR exons"]
            self.tab_writer = csv.DictWriter(self.tab_handle, __fieldnames, delimiter="\t")
            self.tab_writer.writeheader()

        self.out = parsed_args.out
        self.genes = dict()
        self.coding_genes = []
        self.__distances = numpy.array([])
        self.__positions = defaultdict(list)
        self.__coding_positions = defaultdict(list)
        self.__coding_distances = numpy.array([])
        self.__fieldnames = ['Stat', 'Total', 'Average', 'Mode', 'Min',
                             '1%', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', '99%', 'Max']
        self.__rower = csv.DictWriter(self.out,
                                      self.__fieldnames,
                                      delimiter="\t")
        self.__arrays = dict()
        self.__stores = dict()
        self.__prepare_stores()

    def parse_input(self):
        """
        Method to parse the input GTF/GFF file.
        """
        transcript2gene = dict()

        derived_features = set()

        current_gene = None

        for record in self.gff:
            if self.atype == "bed12":
                self.__store_gene(current_gene)
                if not record.parent:
                    record.parent = "{}.gene".format(record.id)
                current_gene = Gene(record, gid=record.parent[0], only_coding=self.only_coding,
                                    logger=self.__logger, use_computer=True)
                transcript2gene[record.id] = record.parent[0]
                current_gene.transcripts[record.id] = TranscriptComputer(record, logger=self.__logger)
            elif record.is_gene is True:
                self.__store_gene(current_gene)
                current_gene = Gene(record,
                                    only_coding=self.only_coding,
                                    logger=self.__logger, use_computer=True)
            elif record.is_transcript is True or (self.atype == GFF3.__annot_type__ and record.feature == "match"):
                if record.parent is None and record.feature != "match":
                    raise TypeError("No parent found for:\n{0}".format(str(record)))
                elif record.feature == "match":
                    gid = record.id
                elif len(record.parent) == 0:
                    self.__logger.warning("No gene found for %s, creating a mock one.", record.id)
                    gid = record.id + ".gene"
                else:
                    gid = record.parent[0]
                    # assert current_gene is not None, record
                if current_gene is None or gid != current_gene.id:
                    # Create a gene record
                    self.__store_gene(current_gene)
                    # if record.feature != "match":
                    #     new_record.feature = "gene"
                    current_gene = Gene(
                        record,
                        gid=gid,
                        only_coding=self.only_coding,
                        logger=self.__logger,
                        use_computer=True)
                transcript2gene[record.id] = gid
                current_gene.transcripts[record.id] = TranscriptComputer(record,
                                                                         logger=self.__logger)

            elif record.is_derived is True and record.is_exon is False:
                derived_features.add(record.id)
            elif record.is_exon is True:
                if self.atype != GFF3.__annot_type__:
                    if "transcript_id" not in record.attributes:
                        # Probably truncated record!
                        record.header = True
                        continue
                    if current_gene is None or record.gene != current_gene.id:
                        self.__store_gene(current_gene)
                        # new_record = record.copy()
                        # new_record.feature = "gene"
                        # new_record.id = new_record.gene
                        current_gene = Gene(
                            record,
                            only_coding=self.only_coding,
                            logger=self.__logger,
                            use_computer=True)
                        record.id = record.transcript
                        transcript2gene[record.transcript] = record.gene
                        current_gene.transcripts[record.transcript] = TranscriptComputer(record,
                                                                                         logger=self.__logger)
                    elif record.transcript not in current_gene:
                        assert record.transcript not in transcript2gene, record.transcript
                        transcript2gene[record.transcript] = record.gene
                        current_gene.transcripts[record.transcript] = TranscriptComputer(record,
                                                                                         logger=self.__logger)
                    else:
                        current_gene.transcripts[record.transcript].add_exon(record)
                elif self.atype == GFF3.__annot_type__ and "cDNA_match" in record.feature:
                    # Here we assume that we only have "cDNA_match" lines, with no parents
                    record.parent = record.id
                    if current_gene is None or record.id != current_gene.id:
                        self.__store_gene(current_gene)
                        # new_record = record.copy()
                        current_gene = Gene(
                            record,
                            gid=record.id,
                            only_coding=self.only_coding,
                            logger=self.__logger,
                            use_computer=True)
                        transcript2gene[record.id] = record.id
                        current_gene.transcripts[record.id] = TranscriptComputer(record,
                                                                                 logger=self.__logger,
                                                                                 )
                    elif record.id not in current_gene:
                        raise ValueError(
                            "cDNA_match instances should not have more than one transcript per \"gene\"!")
                    else:
                        current_gene.transcripts[record.id].add_exon(record)

                elif self.atype == GFF3.__annot_type__:
                    current_gene.add_exon(record)

            elif record.header is True:
                continue
            else:
                continue

        self.__store_gene(current_gene)

    def __call__(self):

        self.parse_input()
        self.gff.close()
        distances = []
        for chromosome in self.__positions:
            __ordered = sorted(self.__positions[chromosome])
            for index, position in enumerate(__ordered[:-1]):
                distances.append(__ordered[index + 1][0] - position[1])

        self.__distances = numpy.array(distances)

        distances = []
        for chromosome in self.__positions:
            __ordered = sorted(self.__coding_positions[chromosome])
            for index, position in enumerate(__ordered[:-1]):
                distances.append(__ordered[index + 1][0] - position[1])

        self.__coding_distances = numpy.array(distances)
        self.writer()
        self.out.close()

    @staticmethod
    def get_stats(row: dict, array: numpy.array) -> dict:
        """
        Method to calculate the necessary statistic from a row of values.
        :param row: the output dictionary row.
        :type row: dict

        :param array: an array of values.

        :rtype : dict
        """

        # Decimal to second digit precision

        if array is None or len(array) == 0:
            return row

        array, weights = array

        if len(weights) > 0 and sum(weights) != 0:
            row["Average"] = "{0:,.2f}".format(round(
                sum(_[0] * _[1] for _ in zip(array, weights)) / sum(weights), 2))

            mode = sklearn.utils.extmath.weighted_mode(array, weights)
            # mode = scipy.stats.mode(array, axis=None)
            row["Mode"] = mode[0][0]
        else:
            row["Average"] = "NA"
            row["Mode"] = "NA"

        keys = ['Min', '1%', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', '99%', 'Max']
        if len(array) == 0:
            quantiles = ["NA"]*len(keys)
        else:
            quantiles = weighted_percentile(array,
                                            [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100],
                                            weights)
        for key, val in zip(keys, quantiles):
            try:
                row[key] = "{0:,.0f}".format(val)  # No decimal
            except KeyError:
                row[key] = val
            except ValueError:
                row[key] = val
            except TypeError:
                row[key] = val
            except Exception:
                raise
        return row

    def __prepare_stores(self):

        self.__stores["genes"] = set()
        self.__stores["coding_genes"] = set()
        self.__stores["monoexonic_genes"] = set()
        self.__stores["transcripts_per_gene"] = dict()
        self.__stores["coding_transcripts_per_gene"] = dict()
        self.__stores["exons"], self.__stores["exons_coding"] = dict(), dict()
        self.__stores["exon_num"], self.__stores["exon_num_coding"] = dict(), dict()
        self.__stores["introns"], self.__stores["introns_coding"] = dict(), dict()
        self.__stores["cds_introns"] = dict()
        self.__stores["cds_exons"] = dict()
        self.__stores["cds_exon_num"] = dict()
        self.__stores["cds_exon_num_coding"] = dict()
        self.__stores["cdna_lengths"] = dict()  # Done
        self.__stores["cdna_lengths_coding"] = dict()
        self.__stores["cds_lengths"] = dict()  # Done
        self.__stores["cds_lengths_coding"] = dict()  # Done
        self.__stores["cds_ratio"] = dict()
        self.__stores["monoexonic_lengths"] = dict()
        self.__stores["multiexonic_lengths"] = dict()
        self.__stores["monocds_lengths"] = dict()

        self.__stores["five_utr_lengths"] = dict()
        self.__stores["five_utr_nums"] = dict()
        self.__stores["three_utr_lengths"] = dict()
        self.__stores["three_utr_nums"] = dict()
        self.__stores["end_distance_from_junction"] = dict()

    def __store_gene(self, gene):

        """

        :param gene:
        :type gene: (None|Mikado.loci.Gene)
        :return:
        """

        if gene is None:
            return

        gene.finalize()
        if self.only_coding is True and gene.is_coding is False:
            return
        if len(gene.transcripts) == 0:
            self.__logger.warning("Gene %s only has invalid, or no, transcripts. Ignoring it.", gene.id)
            return

        self.__positions[gene.chrom].append((gene.start, gene.end, str(gene.strand)))
        self.__stores["genes"].add(gene.id)
        if gene.is_coding is True:
            self.__coding_positions[gene.chrom].append((gene.start, gene.end, gene.strand))
            self.__stores["coding_genes"].add(gene.id)

        self.__stores["transcripts_per_gene"][gene.num_transcripts] = self.__stores["transcripts_per_gene"].get(
            gene.num_transcripts, 0) + 1
        self.__stores["coding_transcripts_per_gene"][
            gene.num_coding_transcripts] = self.__stores["coding_transcripts_per_gene"].get(
            gene.num_coding_transcripts, 0) + 1

        if gene.monoexonic is True:
            self.__stores["monoexonic_genes"].add(gene.id)

        for tid in gene.transcripts:

            if self.tab_handle is not None:
                # __fieldnames = ["TID", "GID",
                #                 "Exon number",
                #                 "cDNA length",
                #                 "CDS length",
                #                 "# coding exons",
                #                 "cDNA/CDS ratio",
                #                 "5'UTR",
                #                 "5'UTR exons",
                #                 "3'UTR",
                #                 "3'UTR exons"]
                row = dict()
                row["TID"] = tid
                row["GID"] = gene.id
                row["Coordinates"] = "{}:{}-{}".format(gene.chrom,
                                                       gene.transcripts[tid].start,
                                                       gene.transcripts[tid].end)
                row["Strand"] = gene.strand
                row["Exon number"] = len(gene.transcripts[tid].exon_lengths)
                row["cDNA length"] = gene.transcripts[tid].cdna_length
                row["Intronic length"] = sum([_[1] - _[0] for _ in gene.transcripts[tid].introns])
                row["CDS length"] = gene.transcripts[tid].selected_cds_length
                row["# coding exons"] = len(gene.transcripts[tid].cds_exon_lengths)
                row["cDNA/CDS ratio"] = "{:.2f}".format(100 * gene.transcripts[tid].selected_cds_length / gene.transcripts[tid].cdna_length)
                row["5'UTR"] = gene.transcripts[tid].five_utr_length
                row["5'UTR exons"] = gene.transcripts[tid].five_utr_num
                row["3'UTR"] = gene.transcripts[tid].three_utr_length
                row["3'UTR exons"] = gene.transcripts[tid].three_utr_num
                self.tab_writer.writerow(row)

            for el in gene.transcripts[tid].exon_lengths:
                self.__stores["exons"][el] = self.__stores["exons"].get(el, 0) + 1
            exon_number = len(gene.transcripts[tid].exon_lengths)
            self.__stores["exon_num"][exon_number] = self.__stores["exon_num"].get(exon_number, 0) + 1
            if exon_number == 1:
                self.__stores["monoexonic_lengths"][
                    gene.transcripts[tid].cdna_length] = self.__stores["monoexonic_lengths"].get(
                    gene.transcripts[tid].cdna_length, 0) + 1
            else:
                self.__stores["multiexonic_lengths"][
                    gene.transcripts[tid].cdna_length] = self.__stores["multiexonic_lengths"].get(
                        gene.transcripts[tid].cdna_length, 0) + 1

            for il in gene.transcripts[tid].intron_lengths:
                self.__stores["introns"][il] = self.__stores["introns"].get(il, 0) + 1

            for cdsil in gene.transcripts[tid].cds_intron_lengths:
                self.__stores["cds_introns"][cdsil] = self.__stores["cds_introns"].get(cdsil, 0) + 1

            for cdsel in gene.transcripts[tid].cds_exon_lengths:
                self.__stores["cds_exons"][cdsel] = self.__stores["cds_exons"].get(cdsel, 0) + 1

            cds_num = len(gene.transcripts[tid].cds_exon_lengths)
            sel_length = gene.transcripts[tid].selected_cds_length
            if (cds_num == 1) or (exon_number == 1 and gene.transcripts[tid].selected_cds_length > 0):
                self.__stores["monocds_lengths"][sel_length] = self.__stores["monocds_lengths"].get(sel_length, 0) + 1

            self.__stores["cds_exon_num"][cds_num] = self.__stores["cds_exon_num"].get(cds_num, 0) + 1
            self.__stores["cdna_lengths"][
                gene.transcripts[tid].cdna_length] = self.__stores["cdna_lengths"].get(
                gene.transcripts[tid].cdna_length, 0) + 1
            self.__stores["cds_lengths"][
                gene.transcripts[tid].selected_cds_length] = self.__stores["cds_lengths"].get(
                gene.transcripts[tid].selected_cds_length, 0) + 1

            if gene.transcripts[tid].selected_cds_length > 0:
                __cds_length = gene.transcripts[tid].selected_cds_length
                __cdna_length = gene.transcripts[tid].cdna_length
                assert __cds_length > 0
                self.__stores["cds_lengths_coding"][
                    gene.transcripts[tid].selected_cds_length] = self.__stores["cds_lengths_coding"].get(
                    gene.transcripts[tid].selected_cds_length, 0) + 1

                for key, attribute in (("five_utr_lengths", "five_utr_length"),
                                       ("three_utr_lengths", "three_utr_length"),
                                       ("five_utr_nums", "five_utr_num"),
                                       ("three_utr_nums", "three_utr_num"),
                                       ("end_distance_from_junction", "end_distance_from_junction")):
                    val = getattr(gene.transcripts[tid], attribute)
                    self.__stores[key][val] = self.__stores[key].get(val, 0) + 1

                cds_ratio = 100 * __cds_length / __cdna_length
                self.__stores["cds_ratio"][cds_ratio] = self.__stores["cds_ratio"].get(cds_ratio, 0) + 1

            if self.only_coding is False:
                if gene.transcripts[tid].selected_cds_length > 0:
                    cdl = gene.transcripts[tid].cdna_length
                    self.__stores["cdna_lengths_coding"][cdl] = self.__stores["cdna_lengths_coding"].get(cdl, 0) + 1
                    for el in gene.transcripts[tid].exon_lengths:
                        self.__stores["exons_coding"][el] = self.__stores["exons_coding"].get(el, 0) + 1
                    num_exons = len(gene.transcripts[tid].exon_lengths)

                    self.__stores["exon_num_coding"][num_exons] = self.__stores["exon_num_coding"].get(num_exons,
                                                                                                       0) + 1
                    cds_num_exons = len(gene.transcripts[tid].cds_exon_lengths)
                    self.__stores["cds_exon_num_coding"][cds_num_exons] = self.__stores[
                                                                              "cds_exon_num_coding"].get(cds_num_exons,
                                                                                                         0) + 1
                    for il in gene.transcripts[tid].intron_lengths:
                        self.__stores["introns_coding"][il] = self.__stores["introns_coding"].get(il, 0) + 1
        return

    def __finalize_arrays(self):

        self.__arrays["Transcripts per gene"] = itemize(self.__stores["transcripts_per_gene"])

        self.__arrays["Coding transcripts per gene"] = itemize(self.__stores["coding_transcripts_per_gene"])
        self.__arrays["Intergenic distances"] = itemize(Counter(self.__distances))
        self.__arrays["Intergenic distances (coding)"] = itemize(Counter(self.__coding_distances))

        self.__arrays['CDNA lengths'] = itemize(self.__stores["cdna_lengths"])
        self.__arrays["CDNA lengths (mRNAs)"] = itemize(self.__stores["cdna_lengths_coding"])
        self.__arrays['CDS lengths'] = itemize(self.__stores["cds_lengths"])
        if self.only_coding is False:
            self.__arrays["CDS lengths (mRNAs)"] = itemize(self.__stores["cds_lengths_coding"])
            self.__arrays['Exons per transcript (mRNAs)'] = itemize(self.__stores["exon_num_coding"])
            self.__arrays['Exon lengths (mRNAs)'] = itemize(self.__stores["exons_coding"])
            self.__arrays["CDS exons per transcript (mRNAs)"] = itemize(
                self.__stores["cds_exon_num_coding"])

        self.__arrays['Monoexonic transcripts'] = itemize(self.__stores["monoexonic_lengths"])
        self.__arrays['MonoCDS transcripts'] = itemize(self.__stores["monocds_lengths"])
        self.__arrays['Exons per transcript'] = itemize(self.__stores["exon_num"])
        self.__arrays['Exon lengths'] = itemize(self.__stores["exons"])
        self.__arrays["Intron lengths"] = itemize(self.__stores["introns"])
        self.__arrays["Intron lengths (mRNAs)"] = itemize(self.__stores["introns_coding"])
        self.__arrays["CDS exons per transcript"] = itemize(self.__stores["cds_exon_num"])
        self.__arrays["CDS exon lengths"] = itemize(self.__stores["cds_exons"])
        self.__arrays["CDS Intron lengths"] = itemize(self.__stores["cds_introns"])
        self.__arrays["5'UTR exon number"] = itemize(self.__stores["five_utr_nums"])
        self.__arrays["3'UTR exon number"] = itemize(self.__stores["three_utr_nums"])
        self.__arrays["5'UTR length"] = itemize(self.__stores["five_utr_lengths"])
        self.__arrays["3'UTR length"] = itemize(self.__stores["three_utr_lengths"])
        self.__arrays["Stop distance from junction"] = itemize(
            self.__stores["end_distance_from_junction"])
        self.__arrays["CDS/cDNA ratio"] = itemize(self.__stores["cds_ratio"])
    # pylint: enable=too-many-locals,too-many-statements

    def writer(self):
        """Method which creates the final output"""

        self.__finalize_arrays()
        self.__rower.writeheader()
        self.__write_statrow('Number of genes',
                             len(self.__stores["genes"]))
        self.__write_statrow("Number of genes (coding)",
                             len(self.__stores["coding_genes"]))
        self.__write_statrow("Number of monoexonic genes",
                             len(self.__stores["monoexonic_genes"])
                             )
        self.__write_statrow('Transcripts per gene',
                             total=numpy.dot)
        self.__write_statrow("Coding transcripts per gene", total=numpy.dot)

        self.__write_statrow('CDNA lengths', total=numpy.dot)
        self.__write_statrow("CDNA lengths (mRNAs)", total=numpy.dot)
        self.__write_statrow('CDS lengths', total=numpy.dot)
        if self.only_coding is False:
            self.__write_statrow("CDS lengths (mRNAs)", total=False)

        self.__write_statrow("CDS/cDNA ratio", total=False)

        self.__write_statrow('Monoexonic transcripts', total=sum)
        self.__write_statrow('MonoCDS transcripts', total=sum)
        self.__write_statrow('Exons per transcript', total=numpy.dot)

        if self.only_coding is False:
            self.__write_statrow('Exons per transcript (mRNAs)',
                                 total="Exon lengths (mRNAs)")
        self.__write_statrow('Exon lengths', total=False)
        if self.only_coding is False:
            self.__write_statrow('Exon lengths (mRNAs)', total=False)
        self.__write_statrow("Intron lengths", total=False)
        if self.only_coding is False:
            self.__write_statrow("Intron lengths (mRNAs)", total=False)

        self.__write_statrow("CDS exons per transcript",
                             total="CDS exon lengths")

        if self.only_coding is False:
            self.__write_statrow("CDS exons per transcript (mRNAs)",
                                 total="CDS exon lengths")

        self.__write_statrow("CDS exon lengths", total=numpy.dot)
        self.__write_statrow("CDS Intron lengths", total=numpy.dot)
        self.__write_statrow("5'UTR exon number", total=sum)
        self.__write_statrow("3'UTR exon number", total=sum)
        self.__write_statrow("5'UTR length", total=numpy.dot)
        self.__write_statrow("3'UTR length", total=numpy.dot)
        self.__write_statrow("Stop distance from junction",
                             total=False)
        self.__write_statrow("Intergenic distances", total=False)
        self.__write_statrow("Intergenic distances (coding)", total=False)

    def __write_statrow(self, stat, total=True):
        """
        Static method to write out a statistic to the
        output file.
        :param stat: the name of the row
        :type stat: str
        :param total: value to display in the "Total" column
        :type total: str | int | sum | bool
        """
        row = dict()
        for key in self.__fieldnames:
            row[key] = "NA"
        row["Stat"] = stat
        if total is False:
            total = "NA"
        elif total is True:
            if len(self.__arrays[stat]) == 2 and isinstance(
                        self.__arrays[stat], numpy.ndarray):
                total = len(self.__arrays[stat][0])
            else:
                total = len(self.__arrays[stat])
        else:
            if total is sum:
                try:
                    _, weights = self.__arrays[stat]
                    total = sum(weights)
                except ValueError:
                    total = "NA"
            elif total is numpy.dot:
                total = numpy.dot(self.__arrays[stat][0],
                                  self.__arrays[stat][1])

            elif not isinstance(total, int):
                assert total in self.__arrays
                if len(self.__arrays[total]) == 2 and isinstance(
                        self.__arrays[total], numpy.ndarray):

                    total = len(self.__arrays[total][0])
                else:
                    total = len(self.__arrays[total])
            else:
                pass  # Just keep the total that was passed from the external

        row["Total"] = total
        if stat in self.__arrays:
            current_array = self.__arrays[stat]
            # assert isinstance(current_array, Counter), type(current_array)
        else:
            current_array = None
        row = self.get_stats(row, current_array)
        self.__rower.writerow(row)
# pylint: enable=too-many-instance-attributes
