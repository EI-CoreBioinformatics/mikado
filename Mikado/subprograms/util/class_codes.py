from collections import OrderedDict as odict

codes = odict()

codes["="] = {"Definition": "Complete intron chain match.",
          "ref_multi": True,
          "pred_multi": True,
          "nucl": "NA",
          "junc": "100%, 100%, 100%",
          "reverse": "=",
          "category": "Match"}
codes["_"] = {"Definition": "Complete match between two monoexonic transcripts.",
          "ref_multi": False,
          "pred_multi": False,
          "nucl": "NA",
          "junc": "NA, NA, >=80%",
          "reverse": "_",
          "category": "Match"}
codes["n"] = {"Definition": """Intron chain extension, ie. both transcripts are multiexonic and
the prediction has novel splice sites outside of the reference transcript boundaries.""",
              "ref_multi": True,
              "pred_multi": True,
              "nucl": "100%, < 100%, < 100%",
              "junc": "100%, < 100%, < 100%",
              "reverse": "c",
              "category": "Extension"}
codes["J"] ={"Definition": """Intron chain extension, ie. both transcripts are multiexonic and
the prediction has novel splice sites inside of the reference transcript boundaries.""",
              "ref_multi": True,
              "pred_multi": True,
              "nucl": "100%, <= 100%, < 100%",
              "junc": "100%, < 100%, < 100%",
              "reverse": "C",
              "category": "Extension"}
codes["c"] = {"Definition": """The prediction is either multiexonic and with its intron chain completely contained
within that of the reference, or monoexonic and contained within one of the reference exons.""",
              "ref_multi": "NA",
              "pred_multi": "NA",
              "nucl": "< 100%, 100%, NA",
              "junc": "< 100%, 100%, NA",
              "reverse": "n",
              "category": "Extension"},
codes["C"] = {"Definition": """The prediction intron chain is completely contained within that of the reference
transcript, but it partially debords either into its introns or outside of the reference boundaries.""",
              "ref_multi": True,
              "pred_multi": True,
              "nucl": "<= 100%, < 100%, < 100%",
              "junc": "< 100%, 100%, < 100%",
              "reverse": "J or j",
              "category": "Extension"}
codes["j"] = {"Definition": """Alternative splicing event.""",
              "ref_multi": True,
              "pred_multi": True,
              "nucl": "NA",
              "junc": "<= 100%, 100%, < 100%",
              "reverse": "j",
              "category": "Alternative splicing"}




#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **j**        | Alternative splicing event.  | True         | True          | NA                | <= 100%, < 100%,  | **j**             | **Alternative     |
#     |              |                              |              |               |                   | < 100%            |                   | splicing**        |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **h**        | Structural match between two | True         | True          | > 0%, > 0%, > 0%  | 0%, 0%, 0%        | **h**             | **Alternative     |
#     |              | models where no splice site  |              |               |                   |                   |                   | splicing**        |
#     |              | is conserved but **at least**|              |               |                   |                   |                   |                   |
#     |              | one intron of the reference  |              |               |                   |                   |                   |                   |
#     |              | and one intron of the        |              |               |                   |                   |                   |                   |
#     |              | prediction partially overlap.|              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **g**        | The monoexonic prediction    | True         | False         | > 0%, > 0%,       | 0%                | **G**             | **Alternative     |
#     | ("mo" before | overlaps one or more exons of|              |               | between 0 and 100%|                   |                   | splicing**        |
#     | release 1)   | the reference transcript; the|              |               |                   |                   |                   |                   |
#     |              | borders of the prediction    |              |               |                   |                   |                   |                   |
#     |              | cannot fall inside the       |              |               |                   |                   |                   |                   |
#     |              | introns of the reference.    |              |               |                   |                   |                   |                   |
#     |              | The prediction transcript    |              |               |                   |                   |                   |                   |
#     |              | can bridge multiple exons    |              |               |                   |                   |                   |                   |
#     |              | of the reference model.      |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **G**        | Generic match of a           | False        | True          | > 0%, > 0%, > 0%  | 0%                | **g**             | **Alternative     |
#     | ("O" before  | multiexonic prediction       |              |               |                   |                   |                   | splicing**        |
#     | release 1)   | transcript versus a          |              |               |                   |                   |                   |                   |
#     |              | monoexonic reference.        |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **o**        | Generic overlap between two  | True         | True          | > 0%, > 0%, > 0%  | 0%, 0%, 0%        | **o**             | **Overlap**       |
#     |              | multiexonic transcripts,     |              |               |                   |                   |                   |                   |
#     |              | which do not share **any**   |              |               |                   |                   |                   |                   |
#     |              | overlap among their introns. |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **e**        | Single exon transcript       | True         | False         | > 0%, > 0%,       | 0%                | **G**             | **Overlap**       |
#     |              | overlapping *one* reference  |              |               | between 0 and 100%|                   |                   |                   |
#     |              | exon and at least 10 bps of a|              |               |                   |                   |                   |                   |
#     |              | reference intron, indicating |              |               |                   |                   |                   |                   |
#     |              | a possible pre-mRNA fragment.|              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **m**        | Generic match between two    | False        | False         | NA, NA, **< 80%** | NA                | **m**             | **Overlap**       |
#     |              | monoexonic transcripts.      |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **i**        | Monoexonic prediction        | True         | False         | 0%                | 0%                | **ri**            | **Intronic**      |
#     |              | completely contained within  |              |               |                   |                   |                   |                   |
#     |              | one intron of the reference  |              |               |                   |                   |                   |                   |
#     |              | transcript.                  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **I**        | Prediction completely        | True         | True          | 0%                | 0%                | **rI**            | **Intronic**      |
#     |              | contained within the introns |              |               |                   |                   |                   |                   |
#     |              | of the reference transcript. |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **rI**       | Reference completely         | True         | True          | 0%                | 0%                | **I**             | **Intronic**      |
#     |              | contained within the introns |              |               |                   |                   |                   |                   |
#     |              | of the prediction transcript.|              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **ri**       | Reverse intron transcript -  | False        | True          | 0%                | 0%                | **i**             | **Intronic**      |
#     |              | the monoexonic reference is  |              |               |                   |                   |                   |                   |
#     |              | completely contained within  |              |               |                   |                   |                   |                   |
#     |              | one intron of the prediction |              |               |                   |                   |                   |                   |
#     |              | transcript.                  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **f**        | Fusion - this special code   | NA           | NA            | **> 10%**, NA, NA | **> 0%**, NA, NA  | NA                | **Fusion**        |
#     |              | is applied when a prediction |              |               |                   |                   |                   |                   |
#     |              | intersects more than one     |              |               |                   |                   |                   |                   |
#     |              | reference transcript. To be  |              |               |                   |                   |                   |                   |
#     |              | considered for fusions,      |              |               |                   |                   |                   |                   |
#     |              | candidate references must    |              |               |                   |                   |                   |                   |
#     |              | **either** share at least one|              |               |                   |                   |                   |                   |
#     |              | splice junction with the     |              |               |                   |                   |                   |                   |
#     |              | prediction, **or** have at   |              |               |                   |                   |                   |                   |
#     |              | least 10% of its bases       |              |               |                   |                   |                   |                   |
#     |              | recalled. If two or more     |              |               |                   |                   |                   |                   |
#     |              | reference transcripts fit    |              |               |                   |                   |                   |                   |
#     |              | these constraints, then the  |              |               |                   |                   |                   |                   |
#     |              | prediction model is          |              |               |                   |                   |                   |                   |
#     |              | classified as a **fusion**.  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **x**        | Monoexonic match on the      | NA           | False         | >= 0%             | 0%                | **x** or **X**    | **Fragment**      |
#     |              | *opposite* strand.           |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **X**        | Multiexonic match on the     | NA           | True          | >= 0%             | 0%                | **x** or **X**    | **Fragment**      |
#     |              | *opposite* strand.           |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **p**        | The prediction is on the same| NA           | NA            | 0%                | 0%                | **p**             | **No overlap**    |
#     |              | strand of a neighbouring but |              |               |                   |                   |                   |                   |
#     |              | non-overlapping transcript.  |              |               |                   |                   |                   |                   |
#     |              | Probable polymerase run-on.  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **P**        | The prediction is on the     | NA           | NA            | 0%                | 0%                | **P**             | **No overlap**    |
#     |              | *opposite* strand of a       |              |               |                   |                   |                   |                   |
#     |              | neighbouring but             |              |               |                   |                   |                   |                   |
#     |              | non-overlapping transcript.  |              |               |                   |                   |                   |                   |
#     |              | Probable polymerase run-on.  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
#     | **u**        | Unknown - no suitable model  | NA           | NA            | 0%                | 0%                | NA                | **No overlap**    |
#     |              | has been found near enough   |              |               |                   |                   |                   |                   |
#     |              | the prediction to perform a  |              |               |                   |                   |                   |                   |
#     |              | comparison.                  |              |               |                   |                   |                   |                   |
#     +--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+







