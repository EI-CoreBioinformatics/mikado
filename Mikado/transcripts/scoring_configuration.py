from marshmallow_dataclass import dataclass, List, Dict, Union, Optional
from dataclasses import field
from marshmallow import validate, post_load, pre_dump, validates_schema
from .transcript import Transcript
available_metrics = Transcript.get_available_metrics()
from ..exceptions import InvalidJson
import re


key_pattern = re.compile(r"([^ ()]+)")


@dataclass
class SizeFilter:
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["gt", "ge", "lt", "le"])})
    value: Union[float, int] = field(metadata={"required": True})
    metric: Optional[str] = field(metadata={"required": False}, default="")


@dataclass
class NumBoolEqualityFilter:
    value: Union[bool, float, int] = field(metadata={"required": True})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["ne", "eq"])})
    metric: Optional[str] = field(metadata={"required": False}, default="")


@dataclass
class InclusionFilter:
    value: list = field(metadata={"required": True})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf(["in", "not in"])})
    metric: Optional[str] = field(metadata={"required": False}, default="")


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

    value: List[Union[int, float]] = field(metadata={
        "required": True,
        "validate": [validate.Length(min=2, max=2), Unique]})
    operator: str = field(metadata={"required": True, "validate": validate.OneOf("within", "not within")})
    metric: Optional[str] = field(metadata={"required": False}, default="")


@dataclass
class MinMaxScore:
    rescaling: str = field(metadata={"required": True, "validate": validate.OneOf(["max", "min"])})
    filter: Optional[Union[SizeFilter, NumBoolEqualityFilter, RangeFilter, InclusionFilter]]
    use_raw: bool = field(default=False)
    multiplier: Union[float, int] = field(default=1, metadata={
        "validate": validate.Range(min=1e6),
        "required": False})
    default: Optional[Union[float, int]] = field(default=0, metadata={"required": False})
    rtype: Optional[str] = field(default="float", metadata={"validate": validate.OneOf(["float", "int", "bool"]),
                                                            "required": False})
    percentage: Optional[bool] = field(default=False)


@dataclass
class TargetScore:
    rescaling: str = field(metadata={"required": True, "validate": validate.OneOf(["target"])})
    value: Union[int, float, bool]
    filter: Optional[Union[SizeFilter, NumBoolEqualityFilter, RangeFilter, InclusionFilter]]
    use_raw: bool = field(default=False)
    multiplier: Union[float, int] = field(default=1, metadata={
        "validate": validate.Range(min=1e6),
        "required": False})
    default: Optional[Union[float, int]] = field(default=0, metadata={"required": False})
    rtype: Optional[str] = field(default="float", metadata={"validate": validate.OneOf(["float", "int", "bool"]),
                                                            "required": False})
    percentage: Optional[bool] = field(default=False)


@dataclass
class Requirements:
    parameters: Dict[str, Union[SizeFilter, InclusionFilter, NumBoolEqualityFilter]] = field(
        metadata={"validate": validate.Length(min=1)}
    )
    expression: List[str] = field(default_factory=lambda: [])

    _compiled = None
    _expression = None

    @pre_dump
    def _remove_eval(self):
        self._compiled = None

    def _create_expression(self):
        if self._expression is not None:
            return
        if len(self.expression) == 0:
            expression = " and ".join(list(self.parameters.keys()))
            keys = self.parameters.keys()
        else:
            expression = " ".join(self.expression)

            keys = set([key for key in key_pattern.findall(expression) if key not in ("and", "or", "not", "xor")])

            diff_params = set.difference(set(keys), set(self.parameters.keys()))

            if len(diff_params) > 0:
                raise InvalidJson(
                    "Expression and required parameters mismatch:\n\t{0}".format(
                        "\n\t".join(list(diff_params))))

        for key in keys:  # Create the final expression
            expression = re.sub(r"\b{}\b".format(key), "evaluated[\"{0}\"]".format(key), expression)
        self._expression = expression

    def _check_parameters(self):
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
            raise InvalidJson(
                "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                    "\n\t".join(parameters_not_found)
                ))

    def _check_my_requirements(self):
        """Check that the provided parameters are as expected"""
        self._expression = self._compiled = None
        self._check_parameters()
        self._check_expression()
        _ = self.compiled

    @property
    def compiled(self):
        if self._compiled is None:
            self._create_expression()
            try:
                compile(self._expression, "<json>", "eval")
            except SyntaxError:
                raise InvalidJson("Invalid expression:\n{}".format(self._expression))
            self._compiled = compile(self._expression, "<json>", "eval")
        return self._compiled

    def __getstate__(self):
        self._compiled = None


@dataclass
class ScoringFile:
    scoring: Dict[str, Union[TargetScore, MinMaxScore]]
    requirements: Requirements
    not_fragmentary: Optional[Requirements]
    as_requirements: Optional[Requirements]
    cds_requirements: Optional[Requirements]
