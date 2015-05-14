import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.json_utils import check_json,to_json
from loci_objects.abstractlocus import abstractlocus
from collections import OrderedDict
import argparse
import csv

def calculate_score(rows, json_dict):
    
    transcripts=OrderedDict((row["tid"],row) for row in rows)

    new_rows=OrderedDict()
    for tid in transcripts:
        new_rows[tid]=dict()
        new_rows[tid]["tid"]=transcripts[tid]["tid"]
        new_rows[tid]["parent"]=transcripts[tid]["parent"]
        new_rows[tid]["original_score"]=transcripts[tid]["score"]
        new_rows[tid]["not_passing"]=False
        new_rows[tid]["recalculated_score"]=0

    
    not_passing=set()
    if "requirements" in json_dict:
        json_dict["requirements"]["compiled"]=compile(json_dict["requirements"]["expression"], "<json>", "eval")
        for tid in transcripts:
            evaluated=dict()
            for key in json_dict["requirements"]["parameters"]:
                name=json_dict["requirements"]["parameters"][key]["name"]
                evaluated[key]=abstractlocus.evaluate(transcripts[tid][name],json_dict["requirements"]["parameters"][key] )
            if eval(json_dict["requirements"]["compiled"]) is False:
                    new_rows[tid]["not_passing"]=True
                    not_passing.add(tid)
#                 if len(not_passing)==len(transcripts): #all transcripts in the locus fail to pass the filter
#                     continue
#                 else:
    
    for param in json_dict["parameters"]:
        rescaling = json_dict["parameters"][param]["rescaling"]
        metrics = [float(transcripts[tid][param]) for tid in transcripts ]
        assert len(metrics)>0, param
        if rescaling=="target":
            target = json_dict["parameters"][param]["value"]
            denominator = max( abs( x-target ) for x in metrics)
        else:
            denominator=(max(metrics)-min(metrics))
        if denominator==0:
            denominator=1
            
        if "requirements" in json_dict and param in json_dict["requirements"]:
            add_original=True
        else:
            add_original=False
            
        for tid in transcripts.keys():
            score=0
            tid_metric = float(transcripts[tid][param])
            if "filter" in json_dict["parameters"][param]:
                    result=abstractlocus.evaluate( tid_metric, json_dict["parameters"][param]["filter"])
                    print("Evaluating param {0} for {1}: {2}".format(param, tid, result), file=sys.stderr)
                    if result is False:
                        score=0
            
            if rescaling == "max":
                ##scoreAM = (rAM - min(rM))/(max(rM)-min(rM)) 
                score = abs( ( tid_metric - min(metrics) ) / denominator )
            elif rescaling=="min":
                score = abs( 1- ( tid_metric - min(metrics) ) / denominator )
            elif rescaling == "target":
                score = 1 - (abs( tid_metric  - target )/denominator )
            score*=json_dict["parameters"][param]["multiplier"]
            
            new_rows[tid]["recalculated_score"]+=score
            new_rows[tid][param]=round(score,2)
            if add_original is True:
                new_rows[tid]["{0}_original".format(param)]=transcripts[tid][param]
    
    for tid in new_rows:
#         if abs(new_rows[tid]["recalculated_score"]-float(new_rows[tid]["original_score"]))>0.5:
#             raise AssertionError(tid, new_rows[tid]["recalculated_score"], float(new_rows[tid]["original_score"]))
        new_rows[tid]["recalculated_score"]=round(new_rows[tid]["recalculated_score"],2)
    
    return new_rows.values()

def main():
    
    parser=argparse.ArgumentParser("Script to convert the metrics file to scores, given a JSON configuration.")
    parser.add_argument("--json_conf", type=argparse.FileType('r'), required=True,
                        help="JSON configuration file. Required.")
    parser.add_argument('metrics', type=argparse.FileType('r'),
                        help='Input metrics file.')
    parser.add_argument('out', nargs='?', default=sys.stdout, type=argparse.FileType('w'),
                        help='Output file. Default STDOUT.')
    args=parser.parse_args()
    
    args.json_conf = to_json(args.json_conf.name)
    check_json(args.json_conf)
    reader=csv.DictReader(args.metrics, delimiter="\t")
    
    fieldnames = ["tid", "parent", "original_score", "recalculated_score", "not_passing"]
    others = set.union(
                        set(args.json_conf["parameters"].keys()),
                        set(args.json_conf["requirements"]["parameters"][key]["name"] for key in args.json_conf["requirements"]["parameters"]) if "requirements" in args.json_conf else set()
                        )


    for key in others:
        assert key in reader.fieldnames, key
    
    requirements_and_score = set.intersection(
                                            set(args.json_conf["parameters"].keys()),
                                            set(args.json_conf["requirements"].keys()) if "requirements" in args.json_conf else set()
                                            )
    final_others=others.copy()
    
    if len(requirements_and_score)>0:
        for key in requirements_and_score:
            final_others.add("{0}_original".format(key))
    
    fieldnames+=sorted(final_others)
    writer=csv.DictWriter(args.out, fieldnames, delimiter="\t")
    writer.writeheader()
    
    if "modules" in args.json_conf:
        import importlib
        for mod in args.json_conf["modules"]:
            globals()[mod]=importlib.import_module(mod)
    
    current_parent=None
    current_rows=[]
    
    for row in reader:
        if row["parent"]!=current_parent:
            if current_parent is not None:
                for new_row in calculate_score(current_rows, args.json_conf):
                    writer.writerow(new_row)
            current_parent=row["parent"]
            current_rows=[]
        current_rows.append(row)
        
    if current_parent is not None:
        for row in calculate_score(current_rows, args.json_conf):
            writer.writerow(row)

    
if __name__=="__main__": main()