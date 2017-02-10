"""
This module contains the definitions of the class codes, using a custom class.
"""


from collections import OrderedDict


def _is_digit(value):
    if not ((value is None) or (isinstance(value, (float, int)) and 0 <= value <= 100)):
        raise ValueError("Invalid numeric value: {}, type: {}".format(value, type(value)))
    return True


def _is_boolean(value):
    if value not in (True, False, None):
        raise ValueError("Invalid boolean value: {}, type: {}".format(value, type(value)))
    return True


class _ClassCode:

    """Container for the class codes ."""

    def __init__(self, code):

        self.__code = code
        self.__definition = None
        self.__ref_multi = None
        self.__pred_multi = None
        self.__nucl_f1, self.__nucl_prec, self.__nucl_rec = None, None, None
        self.__junc_f1, self.__junc_prec, self.__junc_rec = None, None, None
        self.__reverse = None
        self.__category = None

    def __eq__(self, other):
        if hasattr(other, "code") and self.code == other.code:
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.code)

    def __str__(self):
        lines = list()
        lines.append("- Code: {}".format(self.code))
        lines.append("- Definition: {}".format(self.definition))
        lines.append("- Reference multiexonic: {}".format(self.ref_multi))
        lines.append("- Prediction multiexonic: {}".format(self.pred_multi))
        lines.append("- Nucleotide recall, precision, F1: {}".format(self.nucl))
        lines.append("- Junction recall, precision, F1: {}".format(self.junc))
        lines.append("- Reverse class code: {}".format(self.reverse))
        lines.append("- Category: {}".format(self.category))
        return "\n".join(lines)

    @property
    def code(self):
        return self.__code

    @property
    def _nucl_prec(self):
        return self.__nucl_prec

    @_nucl_prec.setter
    def _nucl_prec(self, value):
        self.__nucl_prec = value

    @property
    def _nucl_rec(self):
        return self.__nucl_rec

    @_nucl_rec.setter
    def _nucl_rec(self, value):
        self.__nucl_rec = value

    @property
    def _nucl_f1(self):
        return self.__nucl_f1

    @_nucl_f1.setter
    def _nucl_f1(self, value):
        self.__nucl_f1 = value

    @property
    def nucl(self):
        if self._nucl_f1 is None and self._nucl_rec is None and self._nucl_prec is None:
            return "NA"
        else:
            line = []
            for val in (self._nucl_rec, self._nucl_prec, self._nucl_f1):
                if val is not None:
                    line.append("{}%".format(val))
                else:
                    line.append("NA")
            return ", ".join(line)

    @property
    def junc(self):
        if self._junc_f1 is None and self._junc_rec is None and self._junc_prec is None:
            return "NA"
        else:
            line = []
            for val in (self._junc_rec, self._junc_prec, self._junc_f1):
                if val is not None:
                    line.append("{}%".format(val))
                else:
                    line.append("NA")
            return ", ".join(line)

    @property
    def _junc_prec(self):
        return self.__junc_prec

    @_junc_prec.setter
    def _junc_prec(self, value):
        self.__junc_prec = value

    @property
    def _junc_rec(self):
        return self.__junc_rec

    @_junc_rec.setter
    def _junc_rec(self, value):
        self.__junc_rec = value

    @property
    def _junc_f1(self):
        return self.__junc_f1

    @_junc_f1.setter
    def _junc_f1(self, value):
        self.__junc_f1 = value

    @property
    def pred_multi(self):
        if self.__pred_multi is None:
            return "NA"
        else:
            return self.__pred_multi

    @pred_multi.setter
    def pred_multi(self, value):
        if _is_boolean(value):
            self.__pred_multi = value

    @property
    def ref_multi(self):
        if self.__ref_multi is None:
            return "NA"
        else:
            return self.__ref_multi

    @ref_multi.setter
    def ref_multi(self, value):
        if _is_boolean(value):
            self.__ref_multi = value

    @property
    def definition(self):
        if self.__definition is None:
            return "NA"
        else:
            return self.__definition

    @definition.setter
    def definition(self, value):
        if isinstance(value, bytes):
            value = value.decode()
        elif not (isinstance(value, str) or value is None):
            raise ValueError("Invalid value for definition: {}, type {}".format(value, type(value)))
        self.__definition = value

    @property
    def category(self):
        if self.__category is None:
            return "NA"
        else:
            return self.__category

    @category.setter
    def category(self, value):
        if isinstance(value, bytes):
            value = value.decode()
        elif not (isinstance(value, str) or value is None):
            raise ValueError("Invalid value for category: {}, type {}".format(value, type(value)))
        self.__category = value

    @property
    def reverse(self):
        if self.__reverse is None:
            return "NA"
        else:
            return self.__reverse

    @reverse.setter
    def reverse(self, value):
        if isinstance(value, bytes):
            value = value.decode()
        elif not (isinstance(value, str) or value is None):
            raise ValueError("Invalid value for reverse: {}, type {}".format(value, type(value)))

        self.__reverse = value


class Equal:

    _code = _ClassCode("=")
    _code.definition = "Complete intron chain match."
    _code.pred_multi, _code.ref_multi = True, True
    _code._junc_f1, _code._junc_prec, _code._junc_rec = [100] * 3
    _code.reverse = "="
    _code.category = "Match"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class UnderScore:
    
    _code = _ClassCode("_")
    _code.definition = "Complete match between two monoexonic transcripts."
    _code.ref_multi, _code.pred_multi = False, False
    _code._nucl_f1 = ">=80"
    _code.reverse = "_"
    _code.category = "Match"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)
    
    
class CodeN:

    _code = _ClassCode("n")
    _code.definition = """Intron chain extension, ie. both transcripts are multiexonic and
    the prediction has novel splice sites outside of the reference transcript boundaries."""
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = (100, "< 100", "<100")
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = (100, "< 100", "<100")
    _code.reverse = "c"
    _code.category = "Extension"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalJ:
    
    _code = _ClassCode("J")
    _code.definition = """Intron chain extension, ie. both transcripts are multiexonic and
    the prediction has novel splice sites inside of the reference transcript boundaries."""
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = (100, "<= 100", "<100")
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = (100, "< 100", "<100")
    _code.reverse = "C"
    _code.category = "Extension"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeC:
    _code = _ClassCode("c")
    _code.definition = """The prediction is either multiexonic and with its intron chain completely contained
    within that of the reference, or monoexonic and contained within one of the reference exons."""
    _code.pred_multi, _code.ref_multi = None, None
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "< 100", "100", None
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = "< 100", "100", None
    _code.reverse = "n"
    _code.category = "Extension"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalC:

    _code = _ClassCode("C")
    _code.definition = """The prediction intron chain is completely contained within that of the reference
    transcript, but it partially debords either into its introns or outside of the reference boundaries."""
    _code.pred_multi, _code.ref_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "<= 100", "< 100", "< 100"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = "< 100", "100", "< 100"
    _code.reverse = "J or j"
    _code.category = "Extension"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeJ:
    _code = _ClassCode("j")
    _code.definition = """Alternative splicing event."""
    _code.ref_multi, _code.pred_multi = True, True
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = "<= 100", "100", "< 100"
    _code.reverse = "j or C"
    _code.category = "Alternative splicing"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeH:
    _code = _ClassCode("h")
    _code.definition = """Structural match between two models where where no splice site is conserved but at least
    one intron of the reference and one intron of the prediction partially overlap."""
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 0", "> 0", "> 0"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "h"
    _code.category = "Alternative splicing"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeG:

    _code = _ClassCode("g")
    _code.definition = """The monoexonic prediction overlaps one or more exons of the reference
     transcript; the borders of the prediction cannot fall inside the introns of the reference.
     The prediction transcript can bridge multiple exons of the reference model."""
    _code.ref_multi, _code.pred_multi = True, False
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "G"
    _code.category = "Alternative splicing"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalG:
    _code = _ClassCode("G")
    _code.definition = """Generic match of a multiexonic prediction transcript versus a monoexonic reference."""
    _code.ref_multi, _code.pred_multi = False, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "g"
    _code.category = "Alternative splicing"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeO:
    _code = _ClassCode("o")
    _code.definition = """Generic overlap between two multiexonic transcripts,
    which do not share any overlap among their introns."""
    _code.ref_multi, _code.pred_multi = True, True
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "o"
    _code.category = "Overlap"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeE:
    _code = _ClassCode("e")
    _code.definition = """Single exon transcript overlapping one reference exon and at least 10 bps of a
    reference intron, indicating a possible pre-mRNA fragment."""
    _code.ref_multi, _code.pred_multi = True, False
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "G"
    _code.category = "Overlap"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeM:
    _code = _ClassCode("m")
    _code.definition = """Generic match between two monoexonic transcripts."""
    _code.ref_multi, _code.pred_multi = False, False
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = None, None, "< 80"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = None, None, None
    _code.reverse = "m"
    _code.category = "Overlap"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeI:
    _code = _ClassCode("i")
    _code.definition = "Monoexonic prediction completely contained within one intron of the reference transcript."
    _code.ref_multi, _code.pred_multi = True, False
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "ri"
    _code.category = "Intronic"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalI:
    _code = _ClassCode("I")
    _code.definition = "Prediction completely contained within the introns of the reference transcript."
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "rI"
    _code.category = "Intronic"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeRI:
    _code = _ClassCode("ri")
    _code.definition = """Reverse intron transcript - the monoexonic reference is completely contained
    within one intron of the prediction transcript."""
    _code.ref_multi, _code.pred_multi = False, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "i"
    _code.category = "Intronic"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeRCapitalI:
    _code = _ClassCode("rI")
    _code.definition = """Multiexonic reference completely contained within the introns of the prediction transcript."""
    _code.ref_multi, _code.pred_multi = True, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "I"
    _code.category = "Intronic"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeF:
    _code = _ClassCode("f")
    _code.definition = """Fusion - this special code is applied when a prediction intersects more
    than one reference transcript. To be considered for fusions, candidate references must
    **either** share at least one splice junction with the prediction, **or** have at least 10% of
    its bases recalled. If two or more reference transcripts fit these constraints, then the
    prediction model is classified as a fusion."""
    _code.ref_multi, _code.pred_multi = None, None
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = "> 10", 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = "> 0", 0, 0
    _code.reverse = None
    _code.category = "Fusion"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeX:
    _code = _ClassCode("x")
    _code.definition = "Monoexonic match on the **opposite** strand."
    _code.ref_multi, _code.pred_multi = None, False
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = ">0", ">0", ">0"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "x or X"
    _code.category = "Fragment"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalX:
    _code = _ClassCode("X")
    _code.definition = "Multiexonic match on the **opposite** strand."
    _code.ref_multi, _code.pred_multi = None, True
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = ">0", ">0", ">0"
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = None, None, None
    _code.reverse = "x or X"
    _code.category = "Fragment"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeP:
    _code = _ClassCode("p")
    _code.definition = """The prediction is on the same strand of a neighbouring but non-overlapping transcript.
    Probable polymerase run-on"""
    _code.ref_multi, _code.pred_multi = None, None
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "p"
    _code.category = "Fragment"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeCapitalP:
    _code = _ClassCode("P")
    _code.definition = """The prediction is on the opposite strand of a neighbouring but non-overlapping transcript.
    Probable polymerase run-on."""
    _code.ref_multi, _code.pred_multi = None, None
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = "P"
    _code.category = "Fragment"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


class CodeU:
    _code = _ClassCode("u")
    _code.definition = """Unknown - no suitable model has been found near enough the prediction to
    perform a comparison."""
    _code.ref_multi, _code.pred_multi = None, None
    _code._nucl_rec, _code._nucl_prec, _code._nucl_f1 = 0, 0, 0
    _code._junc_rec, _code._junc_prec, _code._junc_f1 = 0, 0, 0
    _code.reverse = None
    _code.category = "Unknown"

    code, definition = _code.code, _code.definition
    pred_multi, ref_multi = _code.pred_multi, _code.ref_multi
    category, reverse = _code.category, _code.reverse
    nucl, junc = _code.nucl, _code.junc

    __doc__ = str(_code)


codes = OrderedDict()
codes["="] = Equal
codes["_"] = UnderScore
codes["n"] = CodeN
codes["J"] = CodeCapitalJ
codes["c"] = CodeC
codes["C"] = CodeCapitalC
codes["j"] = CodeJ
codes["h"] = CodeH
codes["g"] = CodeG
codes["G"] = CodeCapitalG
codes["o"] = CodeO
codes["e"] = CodeE
codes["m"] = CodeM
codes["i"] = CodeI
codes["I"] = CodeCapitalI
codes["ri"] = CodeRI
codes["rI"] = CodeRCapitalI
codes["f"] = CodeF
codes["x"], codes["X"] = CodeX, CodeCapitalX
codes["p"], codes["P"], codes["u"] = CodeP, CodeCapitalP, CodeU

assert len(set(codes.values())) == len(codes), set.difference(set(codes.keys()),
                                                              set([_.code for _ in codes.values()]))
assert all(_ == codes[_].code for _ in codes)
