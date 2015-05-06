import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.json_utils import check_json,to_json
import argparse
import csv

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
    
    fieldnames = ["tid", "parent", "original_score", "recalculated_score", "not_passing", "offending_keys"]
    others = set.union(
                        set(args.json_conf["parameters"].keys()),
                        set(args.json_conf["requirements"].keys()) if "requirements" in args.json_conf else set()
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
    
    for row in reader:
        new_row=dict()
#         score=int(round(float(row["score"]),0))
        new_row["tid"]=row["tid"]
        new_row["parent"]=row["parent"]
        new_row["original_score"]=row["score"]
        new_row["not_passing"]=False
        new_row["offending_keys"]=[]
        
        new_score = 0
        for key in others:
            if key in ("tid","parent"):
                new_row[key]=row[key]
            elif key not in ("score", "not_passing", "offending_keys"):
                x=row[key]
                if x=="True": x=True
                elif x=="False": x=False
                else: x=float(x)
                if key in args.json_conf["requirements"]:
                    if eval(args.json_conf["requirements"][key]["expression"]) is False:
                        new_row["not_passing"]=True
                        new_row["offending_keys"].append(key)
                if key in args.json_conf["parameters"]:
                    if "expression" in args.json_conf["parameters"][key]:
                        x=round(eval(args.json_conf["parameters"][key]["expression"]),2)
                    x=round(x*args.json_conf["parameters"][key]["multiplier"],2)
                    if key in requirements_and_score:
                        new_row["{0}_original".format(key)]=row[key]
                    new_row[key]=x
                    new_score += new_row[key]
                else:
                    new_row[key]=row[key]
        if new_row["not_passing"] is True:
            new_score=0
        new_score=max(0,int(round(new_score)))
        #assert max(0, new_score) == score, (new_score, score,row["tid"])
        if len(new_row["offending_keys"])>0:
            new_row["offending_keys"]=",".join(new_row["offending_keys"])
        else:
            new_row["offending_keys"]="NA"
        new_row["recalculated_score"] = new_score
        writer.writerow(new_row)
    
if __name__=="__main__": main()