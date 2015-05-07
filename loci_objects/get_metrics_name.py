import os,sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.transcript import transcript

def get_metrics_name(filename=None):
    '''Function to retrieve the available metrics from the metrics.txt
        file, only when the class is called for the first time.
        The method also checks that the requested metrics are defined as
        properties of the transcript class; if that is not the case, it raises
        an exception.'''
    
    if filename is None:
        filename=open(os.path.join( os.path.dirname(__file__),  "metrics.txt"  ))
    else:
        if type(filename) is str:
            assert os.path.exists(filename) and os.path.isfile(filename)
            filename=open(filename)
        elif not hasattr(filename, "mode") or "r" not in filename.mode:
            raise TypeError("Unmanageable object: {0}, type {1}".format(filename, type(filename)) ) 
        
    __available_metrics = [l.rstrip() for l in filename]
    first = ["tid","parent","score"];
    extended=[]
    not_found = []
    not_properties=[]
    for x in __available_metrics:
        if x=='': continue
        if x not in first:
            if not hasattr(transcript, x):
                not_found.append(x)
            elif not type(transcript.__dict__[x]) is property: # @UndefinedVariable
                not_properties.append(x)
            extended.append(x)
    if len(not_found)>0 or len(not_properties)>0:
        err_message=''
        if len(not_found)>0:
            err_message+="The following required metrics are not defined.\n\t{0}\n".format(
                                                                                                             "\n\t".join(not_found)
                                                                                                             )
        if len(not_properties)>0:
            err_message+="""The following metrics are not properties but rather
            methods of the transcript class. I cannot use them properly.
            \t{0}
            """.format("\n\t".join(not_properties))
        raise AttributeError(err_message)
    __available_metrics=first
    __available_metrics.extend(sorted(extended))
    return __available_metrics