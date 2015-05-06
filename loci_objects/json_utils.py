import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.superlocus import superlocus
import json

def check_json(json_conf):
    '''Quick function to check that the JSON dictionary is well formed.'''
    
    parameters_not_found=[]
    parameters_found=set()
    double_parameters=[]
    for parameter in json_conf["parameters"]:
        if parameter not in superlocus.available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.add(parameter)
        parameters_found.add(parameter)
    
    import importlib    
    mods_not_found = [] 
    for mod in json_conf["modules"]:
        try:
            importlib.import_module(mod)
        except ImportError:
            mods_not_found.append(mod)

    if len(parameters_not_found)>0 or len(double_parameters)>0 or len(mods_not_found)>0:
        err_message=''
        if len(parameters_not_found)>0:
            err_message="The following parameters, present in the JSON file, are not available!\n\t{0}\n".format("\n\t".join(parameters_not_found))
        if len(double_parameters)>0:
            err_message+="The following parameters have been specified more than once, please correct:\n\t{0}".format("\n\t".join(list(double_parameters)))
        if len(mods_not_found)>0:
            err_message+="The following requested modules are unavailable:\n\t{0}\n".format("\n\t".join(mods_not_found))
        print(err_message, file=sys.stderr)
        sys.exit(1)
        
def to_json(string):
    
    '''Function to serialize the JSON for configuration and check its consistency.'''
    
    with open(string) as json_file:
        json_dict = json.load(json_file)
        
    import importlib
    if "modules" in json_dict:
        not_found=[]
        for mod in json_dict["modules"]:
            try:
                globals()[mod]=importlib.import_module(mod)
            except ImportError:
                not_found.append(mod)
        if len(not_found)>0:
            raise ImportError("The following modules have not been found:\n{0}".format("\n".join(
                                                                                                 not_found
                                                                                                 )))
        
        
    for param in json_dict["parameters"]:
        
        if "rescaling" not in json_dict["parameters"][param]:
            raise ValueError("No rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(param))
        elif json_dict["parameters"][param]["rescaling"] not in ("max","min", "target"):
            raise ValueError("Invalid rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(param))
        elif json_dict["parameters"][param]["rescaling"]=="target":
            if "value" not in json_dict["parameters"][param]:
                raise ValueError("Target rescaling requested for {0}, but no target value specified. Please specify it with the \"value\" keyword.".format(param))
            json_dict["parameters"][param]["value"]=float(json_dict["parameters"][param]["value"])
        
        if "multiplier" not in json_dict["parameters"][param]:
            json_dict["parameters"][param]["multiplier"]=1
        else:
            json_dict["parameters"][param]["multiplier"]=float(json_dict["parameters"][param]["multiplier"])
            
            
            
    if "requirements" in json_dict:
        for key in json_dict["requirements"]:
            conf = json_dict["requirements"][key]
            assert "value" in conf and "type" in conf, (conf)
            if conf["type"]=="gt":
                oper=">"
            if conf["type"]=="ge":
                oper=">="
            elif conf["type"]=="lt":
                oper="<"
            elif conf["type"]=="le":
                oper="<="
            elif conf["type"]=="eq":
                oper="=="
            elif conf["type"]=="ne":
                oper="!="
            else:
                raise TypeError("Cannot recognize this type: {0}".format(conf["type"]))
            evaluator = compile("x {oper} {value}".format(oper=oper, value=conf["value"]  ),
                                "<json>",
                                "eval",
                                optimize=2)
            json_dict["requirements"][key]["expression"]=evaluator

    return json_dict
