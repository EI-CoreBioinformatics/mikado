import luigi
import sciluigi as sl
import subprocess
import logging
import Mikado
import os


log = logging.getLogger("sciluigi-interface")
log.setLevel("DEBUG")

swissprot = "uniprot_sprot_plants.fasta"


class csWorkFlow(sl.WorkflowTask):

    # task = luigi.Parameter()

    def workflow(self):

        # uncompress = self.new_task("UncompressBlast", UncompressBlast)
        # uncompress.in_data = inp_data.output

        daijin = self.new_task("daijin", Daijin)
        # compare = self.new_task("compareloci", CompareLoci)

        # return compare
        return daijin


class BlastInput(sl.ExternalTask):

    def out(self):

        return sl.TargetInfo(self, "./" + swissprot + ".gz")


class UncompressBlast(sl.Task):

    # def requires(self):
    #     return swissprot + ".gz"

    # swissprot = None

    in_data = BlastInput.out

    def out(self):
        return sl.TargetInfo(self, "uniprot_sprot_plants.fasta")

    def run(self):
        cmd = "gzip -dc {} > {}".format(self.in_data().path, self.out_blast().path)
        print(cmd)
        log.info(cmd)
        subprocess.call(cmd, shell=True)


class Chrom(sl.ExternalTask):

    def out(self):

        return sl.TargetInfo(self, "chr5.fas.gz")


class UncompressChrom(sl.Task):

    data = Chrom.out

    def out(self):

        return sl.TargetInfo(self, "chr5.fas")

    def run(self):
        command = "gzip -dc {} > {}".format(self.data().path, self.out().path)
        subprocess.call(command)


class ConfigureDaijin(sl.Task):

    swissprot = UncompressBlast.out

    def out_put(self):

        return sl.TargetInfo(self, "configuration.yaml")

    def run(self):
        command = "mikado configure --list list.txt --reference chr5.fas --mode permissive --daijin \
                --scoring plants.yaml --junctions junctions.bed -bt {swiss} -bc 1 {configname}".format(
            configname=self.out_put().path, swiss=self.swissprot().path)
        subprocess.call(command, shell=True)


class RawInputFiles(sl.ExternalTask):

    def out_put(self):

        files = []
        for fil in ["class.gtf", "cufflinks.gtf", "stringtie.gtf", "trinity.gff3", "mikado.bed"]:
            files.append(sl.TargetInfo(self, fil))

        # files.append(UncompressChrom.out)

        return files


class Daijin(sl.Task):

    configuration = ConfigureDaijin.out_put
    raw_files = RawInputFiles.out_put
    chrom = UncompressChrom.out

    def out(self):

        # results = {"loci": sl.TargetInfo(self,
        #                                  os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3")),
        #            "subloci": sl.TargetInfo(self,
        #                                  os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado.subloci.gff3")),
        #            "prep": sl.TargetInfo(self,
        #                                  os.path.join("Daijin", "5-mikado", "mikado_prepared.gtf"))}
        return sl.TargetInfo(self, os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"))

    def validate(self):

        config = Mikado.configuration.configurator.to_json(self.configuration().path)

    def run(self):

        command = "daijin mikado --nolock -nd configuration.yaml"
        subprocess.call(command, shell=True)


# class ReferenceGff3(sl.ExternalTask):
#
#     def out(self):
#
#         return sl.TargetInfo(self, "reference.gff3")
#
#
# class IndexReference(sl.Task):
#
#     gff3 = ReferenceGff3.out
#
#     def out(self):
#
#         return sl.TargetInfo(self, "reference.gff3.midx")
#
#     def run(self):
#         log = "index.log"
#         command = "mikado compare -r {reference} --index --log {log}".format(reference=self.gff3().path,
#                                                                              log=log)
#         subprocess.call(command, shell=True)
#
#
# class CompareLoci(sl.Task):
#
#     results = Daijin.out
#     index = IndexReference.out
#     reference = ReferenceGff3.out
#
#     def out(self):
#
#         # results = [sl.TargetInfo(self, "compare_loci.tmap"),
#         #            sl.TargetInfo(self, "compare_loci.refmap"),
#         #            sl.TargetInfo(self, "compare_loci.stats")
#         #            ]
#
#         return sl.TargetInfo(self, "compare_loci.stats")
#
#     def run(self):
#
#         key = token = "loci"
#         command = "mikado compare -r {reference} -p {prediction} -o compare_{token} -l {log}".format(
#             reference=self.reference().path,
#             prediction=self.results()[key].path,
#             token=token,
#             log="compare_{}.log".format(token))
#         subprocess.call(command, shell=True)

# class Compare(sl.Task):
#
#     results = Daijin.out
#     index = IndexReference.out
#     reference = ReferenceGff3.out
#
#     def out(self):
#
#         results = []
#         for key in self.results():
#
#             if key == "prep":
#                 token = "input"
#             else:
#                 token = key
#             for suffix in ["refmap", "tmap", "stats"]:
#                 output = "compare_{}.{}".format(token, suffix)
#                 print(output)
#                 results.append(sl.TargetInfo(self, output))
#         return results
#
#     def run(self):
#
#         for key in self.results():
#             print(key)
#
#             if key == "prep":
#                 token = "input"
#             else:
#                 token = key
#
#             command = "mikado compare -r {reference} -p {prediction} -o compare_{token} -l {log}".format(
#                 reference = self.reference().path,
#                 prediction=self.results()[key].path,
#                 token=token,
#                 log="compare_{}.log".format(token))
#             subprocess.call(command, shell=True)
#


if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=csWorkFlow)
