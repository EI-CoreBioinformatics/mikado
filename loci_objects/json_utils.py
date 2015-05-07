import sys,os.path,re
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.superlocus import superlocus
import json

class UnrecognizedOperator(ValueError):
    pass

class UnrecognizedRescaler(ValueError):
    pass

class InvalidJson(KeyError):
    pass

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
    
    mods_not_found = []
    if "modules" in json_conf:
        import importlib    
        for mod in json_conf["modules"]:
            try:
                importlib.import_module(mod)
            except ImportError:
                mods_not_found.append(mod)

    if "requirements" in json_conf:
        if "parameters" not in json_conf["requirements"]:
            raise InvalidJson("The requirements field must have a \"parameters\" subfield!")
        for key in json_conf["requirements"]["parameters"]:
            key_name=key.split(".")[0]
            if key_name not in superlocus.available_metrics:
                parameters_not_found.append(key_name) 

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
            raise UnrecognizedRescaler("No rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(param))
        elif json_dict["parameters"][param]["rescaling"] not in ("max","min", "target"):
            raise UnrecognizedRescaler("Invalid rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(param))
        elif json_dict["parameters"][param]["rescaling"]=="target":
            if "value" not in json_dict["parameters"][param]:
                raise UnrecognizedRescaler("Target rescaling requested for {0}, but no target value specified. Please specify it with the \"value\" keyword.".format(param))
            json_dict["parameters"][param]["value"]=float(json_dict["parameters"][param]["value"])
        
        if "multiplier" not in json_dict["parameters"][param]:
            json_dict["parameters"][param]["multiplier"]=1
        else:
            json_dict["parameters"][param]["multiplier"]=float(json_dict["parameters"][param]["multiplier"])
            
            
            
    if "requirements" in json_dict:
        if "parameters" not in json_dict["requirements"]:
            raise InvalidJson("The requirements field must have a \"parameters\" subfield!")
        for key in json_dict["requirements"]["parameters"]:
            key_name=key.split(".")[0]

            if "operator" not in json_dict["requirements"]["parameters"][key]:
                raise InvalidJson("No operator provided for requirement {0}".format(key))
            elif "value" not in json_dict["requirements"]["parameters"][key]:
                raise InvalidJson("No value provided for requirement {0}".format(key))
            elif json_dict["requirements"]["parameters"][key]["operator"] not in ("gt","ge","eq","lt","le", "ne","in", "not in"):
                raise UnrecognizedOperator("Unrecognized operator: {0}".format(json_dict["parameters"][param]["operator"]))
            json_dict["requirements"]["parameters"][key]["name"]=key_name
            
        if "expression" not in json_dict["requirements"]:
            json_dict["requirements"]["expression"]=" and ".join(list(json_dict["requirements"]["parameters"].keys()))
            keys=json_dict["requirements"]["parameters"].keys()
            newexpr=json_dict["requirements"]["expression"][:]
        else:
            if type(json_dict["requirements"]["expression"]) is list:
                json_dict["requirements"]["expression"]=" ".join(json_dict["requirements"]["expression"])
            newexpr=json_dict["requirements"]["expression"][:]
            keys = list(filter(lambda x: x not in ("and","or", "not", "xor"), re.findall("([^ ()]+)", json_dict["requirements"]["expression"])))
            diff_params=set.difference(set(keys), set(json_dict["requirements"]["parameters"].keys()))
            if len(diff_params)>0:
                raise InvalidJson("Expression and required parameters mismatch:\n\t{0}".format("\n\t".join(list(diff_params))))

        for key in keys:
            newexpr=re.sub(key, "evaluated[\"{0}\"]".format(key), newexpr)
        json_dict["requirements"]["expression"]=newexpr
        print(newexpr)
        newexpr=compile(newexpr, "<json>", "eval")
        json_dict["requirements"]["compiled"]=newexpr

    return json_dict
