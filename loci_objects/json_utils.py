import sys,os.path,re
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.transcript import transcript
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
    invalid_filter=set()
    available_metrics = transcript.get_available_metrics()
    for parameter in json_conf["parameters"]:
        if parameter not in available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.add(parameter)
        if "filter" in json_conf["parameters"][parameter]:
            conf=json_conf["parameters"][parameter]["filter"]
            if "operator" not in conf or "value" not in conf:
                invalid_filter.add(parameter)
            elif conf["operator"] not in ("gt","ge","eq","lt","le", "ne","in", "not in"):
                invalid_filter.add(parameter)
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
            if key_name not in available_metrics:
                parameters_not_found.append(key_name) 

    if "chimera_split" in json_conf:
        if "blast_check" not in json_conf["chimera_split"] or type(json_conf["chimera_split"]["blast_check"]) is not bool:
            raise InvalidJson("A boolean value must be specified for chimera_split/blast_check!")
        if json_conf["chimera_split"]["blast_check"] is True:
            if "blast" not in  json_conf["chimera_split"]:
                raise InvalidJson("A blast program must be specified to execute the blast check.")
            if "blast_prefix" not in json_conf["chimera_split"]:
                import shutil
                if shutil.which(json_conf["chimera_split"]["blast"]) is None: #@UndefinedVariable
                    raise InvalidJson("I have not been able to find the requested blast program in PATH: {0}".format(json_conf["chimera_split"]["blast"]))
            else:
                if not os.path.exists(os.path.join(json_conf["chimera_split"]["blast_prefix"],json_conf["chimera_split"]["blast"])):
                    raise InvalidJson("I have not been able to find the requested blast program: {0}".format(
                                                                                                             os.path.join(json_conf["chimera_split"]["blast_prefix"],json_conf["chimera_split"]["blast"])
                                                                                                            ))
            if "evalue" not in json_conf["chimera_split"] or type(json_conf["chimera_split"]["evalue"]) is not float:
                raise InvalidJson("I need a maximum e-value for the blast")
            if "database" not in json_conf["chimera_split"] or not os.path.exists(json_conf["chimera_split"]["database"]):
                raise InvalidJson("I need a valid BLAST database!")
 
    if len(parameters_not_found)>0 or len(double_parameters)>0 or len(mods_not_found)>0 or len(invalid_filter)>0:
        err_message=''
        if len(parameters_not_found)>0:
            err_message="The following parameters, present in the JSON file, are not available!\n\t{0}\n".format("\n\t".join(parameters_not_found))
        if len(double_parameters)>0:
            err_message+="The following parameters have been specified more than once, please correct:\n\t{0}".format("\n\t".join(list(double_parameters)))
        if len(mods_not_found)>0:
            err_message+="The following requested modules are unavailable:\n\t{0}\n".format("\n\t".join(mods_not_found))
        if len(invalid_filter)>0:
            err_message+="The following parameters have an invalid filter, please correct:\n\t{0}".format("\n\t".join(list(invalid_filter)))
        print(err_message, file=sys.stderr)
        sys.exit(1)

       
def to_json(string):
    
    '''Function to serialize the JSON for configuration and check its consistency.'''
    
    with open(string) as json_file:
        json_dict = json.load(json_file)
    check_json(json_dict)
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
#         newexpr=compile(newexpr, "<json>", "eval")
#         json_dict["requirements"]["compiled"]=newexpr

    return json_dict
