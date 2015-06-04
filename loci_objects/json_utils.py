import sys,os.path,re
import shutil
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.transcript import transcript
# from loci_objects.exceptions import *
import loci_objects.exceptions
import json
import subprocess

def check_chimera_split(json_conf):

    '''Function to check the "chimera_split" section of the json'''

    assert ("execute" in json_conf["chimera_split"] and type(json_conf["chimera_split"]["execute"]) is bool)
    if json_conf["chimera_split"]["execute"] is True:
        assert ("blast_check" in json_conf["chimera_split"] and type(json_conf["chimera_split"]["blast_check"]) is bool)
        if json_conf["chimera_split"]["blast_check"] is True:
            assert "blast_params" in json_conf["chimera_split"] 
            if "evalue" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["evalue"]) in (float,int) or json_conf["chimera_split"]["blast_params"]["evalue"] is None
            else:
                json_conf["chimera_split"]["blast_params"]["evalue"]=None
            if "hsp_evalue" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["hsp_evalue"]) in (float,int) or json_conf["chimera_split"]["blast_params"]["hsp_evalue"] is None
            else:
                json_conf["chimera_split"]["blast_params"]["hsp_evalue"]=None
            if "max_target_seqs" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["max_target_seqs"]) is int or json_conf["chimera_split"]["blast_params"]["max_target_seqs"] is None
            else:
                json_conf["chimera_split"]["blast_params"]["max_target_seqs"]=None
            if "minimal_hsp_overlap" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"]) is float and \
                    0<=json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"]<=1
            else:
                json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"]=0
            if "lenient" not in json_conf["chimera_split"]["blast_params"]:
                json_conf["chimera_split"]["blast_params"]["lenient"]=False
            else:
                assert type(json_conf["chimera_split"]["blast_params"]["lenient"]) is bool

    return json_conf

def check_scoring(json_conf):
    '''Function to check and format the "scoring" section.'''
    parameters_found=set()
    parameters_not_found=[]
    double_parameters=[]
    invalid_filter=set()
    available_metrics = transcript.get_available_metrics()
    if "scoring" not in json_conf or len(json_conf["scoring"].keys())==0:
        raise loci_objects.exceptions.InvalidJson("No parameters specified for scoring!")
    
    for parameter in json_conf["scoring"]:
        if parameter not in available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.add(parameter)
        if "filter" in json_conf["scoring"][parameter]:
            conf=json_conf["scoring"][parameter]["filter"]
            if "operator" not in conf or "value" not in conf:
                invalid_filter.add(parameter)
            elif conf["operator"] not in ("gt","ge","eq","lt","le", "ne","in", "not in"):
                invalid_filter.add(parameter)
        if "rescaling" not in json_conf["scoring"][parameter]:
            raise loci_objects.exceptions.UnrecognizedRescaler("No rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(parameter))
        elif json_conf["scoring"][parameter]["rescaling"] not in ("max","min", "target"):
            raise loci_objects.exceptions.UnrecognizedRescaler("Invalid rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(parameter))
        elif json_conf["scoring"][parameter]["rescaling"]=="target":
            if "value" not in json_conf["scoring"][parameter]:
                raise loci_objects.exceptions.UnrecognizedRescaler("Target rescaling requested for {0}, but no target value specified. Please specify it with the \"value\" keyword.".format(parameter))
            json_conf["scoring"][parameter]["value"]=float(json_conf["scoring"][parameter]["value"])
        
        if "multiplier" not in json_conf["scoring"][parameter]:
            json_conf["scoring"][parameter]["multiplier"]=1
        else:
            if type(json_conf["scoring"][parameter]["multiplier"]) not in (float,int) or json_conf["scoring"][parameter]["multiplier"]==0:
                raise loci_objects.exceptions.InvalidJson("Invalid multiplier: {0}".format(json_conf["scoring"][parameter]["multiplier"]))
            json_conf["scoring"][parameter]["multiplier"]=float(json_conf["scoring"][parameter]["multiplier"])

    if len(parameters_not_found)>0 or len(double_parameters)>0 or len(invalid_filter)>0:
        err_message=''
        if len(parameters_not_found)>0:
            err_message="The following parameters, present in the JSON file, are not available!\n\t{0}\n".format("\n\t".join(parameters_not_found))
        if len(double_parameters)>0:
            err_message+="The following parameters have been specified more than once, please correct:\n\t{0}".format("\n\t".join(list(double_parameters)))
        if len(invalid_filter)>0:
            err_message+="The following parameters have an invalid filter, please correct:\n\t{0}".format("\n\t".join(list(invalid_filter)))
        raise loci_objects.exceptions.InvalidJson(err_message)
    
    return json_conf

def check_requirements(json_conf):
    
    '''Function to check the requirements section.'''
    available_metrics = transcript.get_available_metrics()
    parameters_not_found=[]
    
    if "requirements" in json_conf:
        #Check that the parameters are valid
        if "parameters" not in json_conf["requirements"]:
            raise loci_objects.exceptions.InvalidJson("The requirements field must have a \"parameters\" subfield!")
        for key in json_conf["requirements"]["parameters"]:
            key_name=key.split(".")[0]
            if key_name not in available_metrics:
                parameters_not_found.append(key_name)
            if "operator" not in json_conf["requirements"]["parameters"][key]:
                raise loci_objects.exceptions.InvalidJson("No operator provided for requirement {0}".format(key))
            elif "value" not in json_conf["requirements"]["parameters"][key]:
                raise loci_objects.exceptions.InvalidJson("No value provided for requirement {0}".format(key))
            elif json_conf["requirements"]["parameters"][key]["operator"] not in ("gt","ge","eq","lt","le", "ne","in", "not in"):
                raise loci_objects.exceptions.UnrecognizedOperator("Unrecognized operator: {0}".format(json_conf["requirements"]["parameters"][key]["operator"]))
            json_conf["requirements"]["parameters"][key]["name"]=key_name
        if len(parameters_not_found)>0:
            raise loci_objects.exceptions.InvalidJson("The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                                                                                                                                 "\n\t".join(parameters_not_found)
                                                                                                                                 ))
        if "expression" not in json_conf["requirements"]: #Create automatically a filtering expression
            json_conf["requirements"]["expression"]=" and ".join(list(json_conf["requirements"]["parameters"].keys()))
            keys=json_conf["requirements"]["parameters"].keys()
            newexpr=json_conf["requirements"]["expression"][:]
        else:
            #Parse the filtering expression, verify that it is syntactically correct
            if type(json_conf["requirements"]["expression"]) is list:
                json_conf["requirements"]["expression"]=" ".join(json_conf["requirements"]["expression"])
            newexpr=json_conf["requirements"]["expression"][:]
            keys = list(filter(lambda x: x not in ("and","or", "not", "xor"), re.findall("([^ ()]+)", json_conf["requirements"]["expression"])))
            diff_params=set.difference(set(keys), set(json_conf["requirements"]["parameters"].keys()))
            if len(diff_params)>0:
                raise loci_objects.exceptions.InvalidJson("Expression and required parameters mismatch:\n\t{0}".format("\n\t".join(list(diff_params))))
        for key in keys: #Create the final expression
            newexpr=re.sub(key, "evaluated[\"{0}\"]".format(key), newexpr)
        json_conf["requirements"]["expression"]=newexpr        
                 
    return json_conf

    
def check_blast(json_conf, json_file):
    '''Function to check the BLAST section of the JSON and eventually perform the database indexing.'''

    assert ("execute" in json_conf["blast"] and type(json_conf["blast"]["execute"]) is bool)
    if json_conf["blast"]["execute"] is False:
        return json_conf
    
    if "program" not in json_conf["blast"]:
        raise loci_objects.exceptions.InvalidJson("No BLAST program specified.") 
    elif os.path.basename(json_conf["blast"]["program"]) not in ("blastn","blastx","tblastx"):
        raise loci_objects.exceptions.InvalidJson("""Invalid BLAST program specified: {0}.
        Supported options: blastn, blastx, tblastx.""")
    if os.path.dirname(json_conf["blast"]["program"])=="":
        program=shutil.which(json_conf["blast"]["program"])
    else:
        try:
            program=os.path.abspath(json_conf["blast"]["program"])
        except OSError:
            program=None
    if program is None:
        raise loci_objects.exceptions.InvalidJson("The selected BLAST program {0} has not been found on this system!".format(json_conf["blast"]["program"]))
    json_conf["blast"]["program"]=program
         
    if "evalue" not in json_conf["blast"]:
        json_conf["blast"]["evalue"]=10
    else:
        if type(json_conf["blast"]["evalue"]) not in (float,int) or \
            0>json_conf["blast"]["evalue"]:
                raise loci_objects.exceptions.InvalidJson("Invalid evalue: {0}".format(json_conf["blast"]["evalue"]))
    if "max_target_seqs" in json_conf["blast"]:
        assert type(json_conf["blast"]["max_target_seqs"]) is int
    if "database" not in json_conf["blast"]:
        raise loci_objects.exceptions.InvalidJson("No BLAST database provided!")
    json_conf["blast"]["database"]=os.path.abspath(json_conf["blast"]["database"])
    if not os.path.exists(json_conf["blast"]["database"]):
        db=json_conf["blast"]["database"]
        if os.path.dirname(json_conf["blast"]["database"])=="":
            db=os.path.join(
                                                        os.path.dirname(json_file),
                                                        json_conf["blast"]["database"]
                                                        )
        if not os.path.exists(db):
            raise loci_objects.exceptions.InvalidJson("I need a valid BLAST database! This file does not exist:\n{0}".format(json_conf["blast"]["database"]))
        else:
            json_conf["blast"]["database"]=os.path.abspath(db)
    else:
        json_conf["blast"]["database"]=os.path.abspath(json_conf["blast"]["database"])
    makeblastdb_cmd = os.path.join(os.path.dirname(json_conf["blast"]["program"]), "makeblastdb")
    assert os.path.exists(makeblastdb_cmd)
    retcode=0
    if os.path.basename(json_conf["blast"]["program"])=="blastx" and not os.path.exists("{0}.pog".format(json_conf["blast"]["database"])):
        retcode=subprocess.call("{0} -in {1} -dbtype prot -parse_seqids".format(makeblastdb_cmd,json_conf["blast"]["database"]),
                                shell=True)
    elif os.path.basename(json_conf["blast"]["program"]) in ("blastn","tblastx") and not os.path.exists("{0}.nog".format(json_conf["blast"]["database"])):
        retcode=subprocess.call("{0} -in {1} -dbtype nucl -parse_seqids".format(makeblastdb_cmd, json_conf["blast"]["database"]),
                                shell=True)
    if retcode!=0:
        raise OSError("BLAST indexing failed.")

    return json_conf

def check_orf_loading(json_conf):
    
    if "orf_loading" not in json_conf:
        json_conf["orf_loading"]=dict()
        json_conf["orf_loading"]["strand_specific"]=False
        json_conf["orf_loading"]["minimal_secondary_orf_length"]=0
    else:
        if "strand_specific" not in json_conf:
            json_conf["orf_loading"]["strand_specific"]=False
        else:
            if not type(json_conf["orf_loading"]["strand_specific"]) is bool:
                raise loci_objects.exceptions.InvalidJson("Invalid strand_specific value: {0}".format(json_conf["orf_loading"]["strand_specific"]))
        if "minimal_secondary_orf_length" not in json_conf["orf_loading"]:
            json_conf["orf_loading"]["minimal_secondary_orf_length"]=0
        else:
            if not type(json_conf["orf_loading"]["minimal_secondary_orf_length"]) is int:
                raise loci_objects.exceptions.InvalidJson("Invalid minimal_secondary_orf_length value: {0}".format(json_conf["orf_loading"]["minimal_secondary_orf_length"]))

    return json_conf

def check_json(json_conf, json_file):
    '''Quick function to check that the JSON dictionary is well formed.'''
    
    if "db" not in json_conf:
        raise loci_objects.exceptions.InvalidJson("No database specified.") 
    if "dbtype" not in json_conf:
        raise loci_objects.exceptions.InvalidJson("DB type not specified.")
    if json_conf["dbtype"] not in ("sqlite", "mysql", "psql"):
        raise loci_objects.exceptions.InvalidJson("Invalid DB type: {0}. At the moment we support sqlite, mysql, psql".format(json_conf["dbtype"]))
        
    if "strand_spefific" not in json_conf:
        json_conf["strand_specific"]=False
    else:
        if not type(json_conf["strand_specific"]) is bool:
            raise loci_objects.exceptions.InvalidJson("Invalid boolean value for \"strand_specific\": {0}".format(json_conf["strand_specific"]))
        
    json_conf = check_blast(json_conf, json_file)
    json_conf = check_requirements(json_conf)
    json_conf = check_scoring(json_conf)
    json_conf = check_orf_loading(json_conf)
    json_conf = check_chimera_split(json_conf)
    return json_conf
       
def to_json(string):
    
    '''Function to serialise the JSON for configuration and check its consistency.'''
    
    with open(string) as json_file:
        json_dict = json.load(json_file)
    json_dict=check_json(json_dict, json_file.name)
    return json_dict
