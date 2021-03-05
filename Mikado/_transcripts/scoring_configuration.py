import functools

from marshmallow_dataclass import dataclass, List, Dict, Union, Optional
from dataclasses import field, asdict
from marshmallow import validate
from .transcript_base import TranscriptBase
available_metrics = TranscriptBase.get_available_metrics()
from ..exceptions import InvalidConfiguration
import re


key_pattern = re.compile(r"([^ ()]+)")


# DO NOT TOUCH THIS. The compilation at each call, no lazy storing, is necessary to prevent
# multiprocessing madness (eval and pickle do not mix at all). The LRU cache hopefully speeds things up
# in the context of a single process
@functools.lru_cache(maxsize=1000, typed=True)
def compiler(expression):
    try:
        return compile(expression, "<json>", "eval")
    except SyntaxError:
        raise InvalidConfiguration("Invalid expression:\n{}".format(expression))


@dataclass
class SizeFilter:
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["gt", "ge", "lt", "le"])})
    value: float = field(metadata={"required": True})
    metric: Optional[str] = field(metadata={"required": False}, default=None)
    name: Optional[str] = field(default=None)
    source: Optional[str] = field(default=None)


@dataclass
class NumBoolEqualityFilter:
    value: Union[float, bool] = field(metadata={"required": True})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["ne", "eq"])})
    metric: Optional[str] = field(metadata={"required": False}, default=None)
    name: Optional[str] = field(default=None)
    source: Optional[str] = field(default=None)


@dataclass
class InclusionFilter:
    value: list = field(metadata={"required": True})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["in", "not in"])})
    metric: Optional[str] = field(metadata={"required": False}, default=None)
    name: Optional[str] = field(default=None)
    source: Optional[str] = field(default=None)


@dataclass
class RangeFilter:
    class Unique(validate.Validator):
        message = "{input} is Not unique"

        def __init__(self, error=None):
            self.error = error

        def _repr_args(self):
            return ""

        def _format_error(self, value):
            return self.message.format(input=value)

        def __call__(self, value):
            if len(value) != len(set(value)):
                raise validate.ValidationError(self._format_error(value))
            return value

    value: List[float] = field(metadata={
        "required": True,
        "validate": [validate.Length(min=2, max=2), Unique]})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["within", "not within"])})
    metric: Optional[str] = field(metadata={"required": False}, default=None)
    name: Optional[str] = field(default=None)


@dataclass
class MinMaxScore:
    rescaling: str = field(metadata={"required": True, "validate": validate.OneOf(["max", "min"])})
    filter: Optional[Union[SizeFilter, NumBoolEqualityFilter, RangeFilter, InclusionFilter]]
    use_raw: bool = field(default=False)
    multiplier: float = field(default=1, metadata={"required": False})
    default: Optional[Union[float, bool]] = field(default=0, metadata={"required": False})
    rtype: Optional[str] = field(default="float", metadata={"validate": validate.OneOf(["float", "int", "bool"]),
                                                            "required": False})
    percentage: Optional[bool] = field(default=False)


@dataclass
class TargetScore:
    rescaling: str = field(metadata={"required": True, "validate": validate.OneOf(["target"])})
    value: Union[float, bool] = field(metadata={"required": True})
    filter: Optional[Union[SizeFilter, NumBoolEqualityFilter, RangeFilter, InclusionFilter]]
    # Use_raw must be false for target scores
    use_raw: bool = field(default=False, metadata={"validate": validate.OneOf([False])})
    multiplier: float = field(default=1, metadata={"required": False})
    default: Optional[Union[float, bool]] = field(default=0, metadata={"required": False})
    rtype: Optional[str] = field(default="float", metadata={"validate": validate.OneOf(["float", "int", "bool"]),
                                                            "required": False})
    percentage: Optional[bool] = field(default=False)


@dataclass
class Requirements:
    parameters: Dict[str, Union[SizeFilter, InclusionFilter, NumBoolEqualityFilter, RangeFilter]] = field(
        metadata={"validate": validate.Length(min=1), "required": True}
    )
    expression: List[str] = field(default_factory=lambda: [])

    _expression = None

    @staticmethod
    def _create_expression(expression: list, parameters: dict):
        if len(expression) == 0:
            expression = " and ".join(list(parameters.keys()))
            keys = parameters.keys()
        else:
            expression = " ".join(expression)
            keys = set([key for key in key_pattern.findall(expression) if key not in ("and", "or", "not", "xor")])
            diff_params = set.difference(set(keys), set(parameters.keys()))
            if len(diff_params) > 0:
                raise InvalidConfiguration("Expression and required parameters mismatch:\n\t{0}".format(
                    "\n\t".join(list(diff_params))))
        for key in keys:  # Create the final expression
            expression = re.sub(r"\b{}\b".format(key), "evaluated[\"{0}\"]".format(key), expression)
        return expression

    def _check_parameters(self):
        # Check that the parameters are valid
        parameters_not_found = []

        for key in self.parameters:
            key_value = None
            dots = key.split(".")
            if dots[0] == "external" or dots[0] == "attributes":
                if len(dots) == 1 or len(dots) > 3:
                    parameters_not_found.append(dots[0])
                    continue
                elif len(dots) == 2:
                    key_name = ".".join(dots)
                else:
                    key_name = ".".join(dots[:-1])
                    # print(key_name)
                key_value = dots[1]
            else:
                key_name = dots[0]
            if key_name not in available_metrics:
                if key_value is not None:
                    pass
                else:
                    parameters_not_found.append(key_name)
                    continue
            self.parameters[key].name = key_name
            if key_value:
                self.parameters[key].source = key_value

        if len(parameters_not_found) > 0:
            raise InvalidConfiguration(
                "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                    "\n\t".join(parameters_not_found)
                ))

    # If the requirements are changed at run time, this *needs* to be called again
    def _check_my_requirements(self):
        """Check that the provided parameters are as expected"""
        self._expression = None
        assert hasattr(self, "parameters"), (type(self), self.__dict__)
        self._check_parameters()
        self._expression = self._create_expression(self.expression, self.parameters)
        _ = self.compiled

    @property
    def compiled(self):
        return compiler(self._expression)

    def copy(self):
        return self.Schema().load(asdict(self))


@dataclass
class ScoringFile:
    scoring: Dict[str, Union[TargetScore, MinMaxScore]] = field(metadata={"validate": validate.Length(min=1)})
    requirements: Requirements
    not_fragmentary: Optional[Requirements] = field(default=None)
    as_requirements: Optional[Requirements] = field(default=None)
    cds_requirements: Optional[Requirements] = field(default=None)

    def _create_missing_reqs(self, minimal_orf_length):
        """"""

        assert self.requirements.parameters

        for section_name in ["not_fragmentary", "as_requirements", "cds_requirements"]:
            section = getattr(self, section_name)
            if section is not None:
                continue
            if section_name == "cds_requirements":
                # If no CDS requirement section is present, presume that the only req. is that the ORF
                # is equal or longer than the "minimal_orf_length" (taken from pick.orf_loading)
                mock = {
                    "parameters": {"selected_cds_length": {
                        "operator": "ge",
                        "value": minimal_orf_length}},
                    "expression": ["selected_cds_length"]
                }
                self.cds_requirements = Requirements.Schema().load(mock)
                assert self.cds_requirements.parameters
            else:
                setattr(self, section_name, self.requirements.copy())

    def _check_scoring(self):
        """Method to check that the specified metrics are correct"""
        parameters_found = set()
        parameters_not_found = []
        double_parameters = []
        invalid_raw = set()
        invalid_filter = set()

        for parameter in self.scoring:
            if parameter not in available_metrics:
                # Leniency for external_scores
                if "." in parameter and len(parameter.split(".")) == 2 and (parameter.split(".")[0] == "external" or
                                                                            parameter.split(".")[0] == "attributes"):
                    pass
                else:
                    parameters_not_found.append(parameter)

            if parameter in parameters_found:
                double_parameters.append(parameter)

            if (not parameter.startswith("external") and not parameter.startswith("attributes") and
                    getattr(TranscriptBase, parameter).usable_raw is False and
                    self.scoring[parameter].use_raw is True):
                invalid_raw.add(parameter)

            if self.scoring[parameter].filter:
                metric = self.scoring[parameter].filter.metric
                if metric is not None and metric != "" and metric not in available_metrics:
                    parameters_not_found.append(metric)

        if len(parameters_not_found) > 0 or len(double_parameters) > 0 or len(invalid_filter) > 0 or len(
                invalid_raw) > 0:
            err_message = ''
            if len(parameters_not_found) > 0:
                err_message = """The following parameters, present in the JSON file,
                    are not available!\n\t{0}\n""".format(
                    "\n\t".join(parameters_not_found))
            if len(double_parameters) > 0:
                err_message += """The following parameters have been specified more than once,
                    please correct:\n\t{0}""".format("\n\t".join(list(double_parameters)))
            if len(invalid_filter) > 0:
                err_message += """The following parameters have an invalid filter,
                    please correct:
                    \t{0}""".format("\n\t".join(list(invalid_filter)))
            if len(invalid_raw) > 0:
                err_message += """The following parameters cannot be used as raw scores, either
                    because they are not normalized on their own, or because a "target" rescaling has been asked for:
                    \t{0}""".format("\n\t".join(list(invalid_raw)))
            raise InvalidConfiguration(err_message)

    def check(self, minimal_orf_length):
        self.requirements._check_my_requirements()
        self._create_missing_reqs(minimal_orf_length=minimal_orf_length)
        for name, section in zip(["as_requirements", "cds_requirements", "not_fragmentary"],
                                            [self.as_requirements, self.cds_requirements, self.not_fragmentary]):
            assert section is not None, name
            section._check_my_requirements()
        self._check_scoring()
