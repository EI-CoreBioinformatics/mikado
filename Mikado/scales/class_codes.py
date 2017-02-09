"""
This module contains the definitions of the class codes, using a custom class.
"""


from collections import OrderedDict as odict

def _is_digit(value):
    if not ((value is None) or (isinstance(value, (float, int)) and 0 <= value <= 100)):
        raise ValueError("Invalid numeric value: {}, type: {}".format(value, type(value)))
    return True


def _is_boolean(value):
    if value not in (True, False, None):
        raise ValueError("Invalid boolean value: {}, type: {}".format(value, type(value)))
    return True


class ClassCode:

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


def code_equal():
    equal = ClassCode("=")
    equal.definition = "Complete intron chain match."
    equal.pred_multi, equal.ref_multi = True, True
    equal._junc_f1, equal._junc_prec, equal._junc_rec = [100] * 3
    equal.reverse = "="
    equal.category = "Match"
    return equal


def code_underscore():
    underscore = ClassCode("_")
    underscore.definition = "Complete match between two monoexonic transcripts."
    underscore.ref_multi, underscore.pred_multi = False, False
    underscore._nucl_f1 = ">=80"
    underscore.reverse = "_"
    underscore.category = "Match"
    return underscore


def code_n():
    code = ClassCode("n")
    code.definition = """Intron chain extension, ie. both transcripts are multiexonic and
    the prediction has novel splice sites outside of the reference transcript boundaries."""
    code.ref_multi, code.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = (100, "< 100", "<100")
    code._junc_rec, code._junc_prec, code._junc_f1 = (100, "< 100", "<100")
    code.reverse = "c"
    code.category = "Extension"
    return code


def code_capital_j():
    code = ClassCode("J")
    code.definition = """Intron chain extension, ie. both transcripts are multiexonic and
    the prediction has novel splice sites inside of the reference transcript boundaries."""
    code.ref_multi, code.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = (100, "<= 100", "<100")
    code._junc_rec, code._junc_prec, code._junc_f1 = (100, "< 100", "<100")
    code.reverse = "C"
    code.category = "Extension"
    return code


def code_c():
    code = ClassCode("c")
    code.definition = """The prediction is either multiexonic and with its intron chain completely contained
    within that of the reference, or monoexonic and contained within one of the reference exons."""
    code.pred_multi, code.ref_multi = None, None
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "< 100", "100", None
    code._junc_rec, code._junc_prec, code._junc_f1 = "< 100", "100", None
    code.reverse = "n"
    code.category = "Extension"
    return code


def code_capital_c():
    code = ClassCode("C")
    code.definition = """The prediction intron chain is completely contained within that of the reference
    transcript, but it partially debords either into its introns or outside of the reference boundaries."""
    code.pred_multi, code.ref_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "<= 100", "< 100", "< 100"
    code._junc_rec, code._junc_prec, code._junc_f1 = "< 100", "100", "< 100"
    code.reverse = "J or j"
    code.category = "Extension"
    return code


def code_j():
    code = ClassCode("j")
    code.definition = """Alternative splicing event."""
    code.ref_multi, code.pred_multi = True, True
    code._junc_rec, code._junc_prec, code._junc_f1 = "<= 100", "100", "< 100"
    code.reverse = "j or C"
    code.category = "Alternative splicing"
    return code


def code_h():
    code = ClassCode("h")
    code.definition = """Structural match between two models where where no splice site is conserved but at least
    one intron of the reference and one intron of the prediction partially overlap."""
    code.ref_multi, code_n.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 0", "> 0", "> 0"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "h"
    code.category = "Alternative splicing"
    return code


def code_g():
    code = ClassCode("g")
    code.definition = """The monoexonic prediction overlaps one or more exons of the reference transcript;
    the borders of the prediction cannot fall inside the introns of the reference.
    The prediction transcript can bridge multiple exons of the reference model."""
    code.ref_multi, code.pred_multi = True, False
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "G"
    code.category = "Alternative splicing"
    return code


def code_capital_g():
    code = ClassCode("G")
    code.definition = """Generic match of a multiexonic prediction transcript versus a monoexonic reference."""
    code.ref_multi, code.pred_multi = False, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "g"
    code.category = "Alternative splicing"
    return code


def code_o():
    code = ClassCode("o")
    code.definition = """Generic overlap between two multiexonic transcripts,
    which do not share any overlap among their introns."""
    code.ref_multi, code.pred_multi = True, True
    code.ref_multi, code.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "o"
    code.category = "Overlap"
    return code


def code_e():
    code = ClassCode("e")
    code.definition = """Single exon transcript overlapping one reference exon and at least 10 bps of a
    reference intron, indicating a possible pre-mRNA fragment."""
    code.ref_multi, code.pred_multi = True, False
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 0", "> 0", "0% < F1 < 100"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "G"
    code.category = "Overlap"
    return code


def code_m():
    code = ClassCode("m")
    code.definition = """Generic match between two monoexonic transcripts."""
    code.ref_multi, code.pred_multi = False, False
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = None, None, "< 80"
    code._junc_rec, code._junc_prec, code._junc_f1 = None, None, None
    code.reverse = "m"
    code.category = "Overlap"
    return code


def code_i():
    code = ClassCode("i")
    code.definition = "Monoexonic prediction completely contained within one intron of the reference transcript."
    code.ref_multi, code.pred_multi = True, False
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "ri"
    code.category = "Intronic"
    return code


def code_capital_i():
    code = ClassCode("I")
    code.definition = "Prediction completely contained within the introns of the reference transcript."
    code.ref_multi, code.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "rI"
    code.category = "Intronic"
    return code


def code_r_i():
    code = ClassCode("ri")
    code.definition = """Reverse intron transcript - the monoexonic reference is completely contained
    within one intron of the prediction transcript."""
    code.ref_multi, code.pred_multi = False, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "i"
    code.category = "Intronic"
    return code


def code_r_capital_i():
    code = ClassCode("rI")
    code.definition = """Multiexonic reference completely contained within the introns of the prediction transcript."""
    code.ref_multi, code.pred_multi = True, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "I"
    code.category = "Intronic"
    return code


def code_f():
    code = ClassCode("f")
    code.definition = """Fusion - this special code is applied when a prediction intersects more than one
    reference transcript. To be considered for fusions, candidate references must **either** share at least one
    splice junction with the prediction, **or** have at least 10% of its bases recalled.
    If two or more reference transcripts fit these constraints, then the prediction model is classified as a fusion."""
    code.ref_multi, code.pred_multi = None, None
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = "> 10", 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = "> 0", 0, 0
    code.reverse = None
    code.category = "Fusion"
    return code


def code_x():
    code = ClassCode("x")
    code.definition = "Monoexonic match on the **opposite** strand."
    code.ref_multi, code.pred_multi = None, False
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = ">0", ">0", ">0"
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "x or X"
    code.category = "Fragment"
    return code


def code_capital_x():
    code = ClassCode("X")
    code.definition = "Multiexonic match on the **opposite** strand."
    code.ref_multi, code.pred_multi = None, True
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = ">0", ">0", ">0"
    code._junc_rec, code._junc_prec, code._junc_f1 = None, None, None
    code.reverse = "x or X"
    code.category = "Fragment"
    return code


def code_p():
    code = ClassCode("p")
    code.definition = """The prediction is on the same strand of a neighbouring but non-overlapping transcript.
    Probable polymerase run-on"""
    code.ref_multi, code.pred_multi = None, None
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "p"
    code.category = "Fragment"
    return code


def code_capital_p():
    code = ClassCode("P")
    code.definition = """The prediction is on the opposite strand of a neighbouring but non-overlapping transcript.
    Probable polymerase run-on."""
    code.ref_multi, code.pred_multi = None, None
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = "P"
    code.category = "Fragment"
    return code


def code_u():
    code = ClassCode("u")
    code.definition = """Unknown - no suitable model has been found near enough the prediction to
    perform a comparison."""
    code.ref_multi, code.pred_multi = None, None
    code._nucl_rec, code._nucl_prec, code._nucl_f1 = 0, 0, 0
    code._junc_rec, code._junc_prec, code._junc_f1 = 0, 0, 0
    code.reverse = None
    code.category = "Unknown"
    return code


codes = odict()
codes["="] = code_equal()
codes["_"] = code_underscore()
codes["n"] = code_n()
codes["J"] = code_capital_j()
codes["c"] = code_c()
codes["C"] = code_capital_c()
codes["j"] = code_j()
codes["h"] = code_h()
codes["g"] = code_g()
codes["G"] = code_capital_g()
codes["o"] = code_o()
codes["e"] = code_e()
codes["m"] = code_m()
codes["i"] = code_i()
codes["I"] = code_capital_i()
codes["ri"] = code_r_i()
codes["rI"] = code_r_capital_i()
codes["f"] = code_f()
codes["x"], codes["X"] = code_x(), code_capital_x()
codes["p"], codes["P"], codes["u"] = code_p(), code_capital_p(), code_u()

assert len(set(codes.values())) == len(codes), set.difference(set(codes.keys()), set([_.code for _ in codes.values()]))
assert all(_ == codes[_].code for _ in codes)