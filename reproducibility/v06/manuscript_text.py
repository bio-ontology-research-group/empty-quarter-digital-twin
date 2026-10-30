"""Normalisation and numeric-literal extraction for the LaTeX sources (used by verify_manuscript.py).

`normalise` turns a .tex file into plain-ish text in which the registry templates are matched; `mask` returns a text
of identical length in which everything that is not a reportable quantity (equation constants, identifiers, dates,
years, DOIs, sub/superscripts, section/table labels) is blanked, so literals found in the mask map back to offsets
in the normalised text.
"""
import re

NUMP = r"-?(?:\d[\d,]*(?:\.\d+)?(?:E-?\d+)?|\.\d+)%?"
NUM = re.compile(r"(?<![\w.\\])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:E-?\d+)?%?(?![A-Za-z\d])|(?<![\w.\\])\.\d+%?(?![A-Za-z\d])")
MONTHS = "January|February|March|April|May|June|July|August|September|October|November|December"
SKIP_ENV = ("equation", "align", "verbatim", "thebibliography")
SKIP_CMD_ARG = ("cite", "label", "ref", "eqref", "includegraphics", "input", "href", "texttt", "bibliography",
                "bibliographystyle", "usepackage", "documentclass", "renewcommand", "newcommand", "section", "subsection",
                "url", "doi", "title")
FORMULA = re.compile(r"/|<|>|\+|\\leq|\\geq|\\log|\\min|\\max|\\lceil|\\rightarrow|\\times|\(|\d(?=[A-Za-z]|\\(?!pm|%|circ))")


def strip_comments(s):
    return re.sub(r"(?<!\\)%.*", "", s)


def normalise(s):
    s = strip_comments(s)
    if "\\begin{document}" in s:
        s = s.split("\\begin{document}", 1)[1]
    for env in SKIP_ENV:
        s = re.sub(r"\\begin\{%s\*?\}.*?\\end\{%s\*?\}" % (env, env), " ", s, flags=re.S)
    s = re.sub(r"\\begin\{(?:tabular|longtable)\*?\}(?:\[[^\]]*\])?\{(?:[^{}]|\{[^{}]*\})*\}", " ", s)
    s = re.sub(r"\\begin\{(figure|table)\*?\}\[[^\]]*\]", r"\\begin{\1}", s)
    for cmd in SKIP_CMD_ARG:
        s = re.sub(r"\\%s\*?(?:\[[^\]]*\])?\{(?:[^{}]|\{[^{}]*\})*\}" % cmd, " ", s)
    s = re.sub(r"\\includegraphics\[[^\]]*\]", " ", s)
    s = re.sub(r"width=\\textwidth", " ", s)
    s = re.sub(r"10\^\{(-?\d+)\}|10\^(-?\d)", lambda m: "1E" + (m.group(1) or m.group(2)), s)
    s = s.replace("\\%", "%").replace("~", " ").replace("\\,", " ").replace("\\;", " ")
    return s


def _blank(m):
    return " " * len(m.group(0))


def _math(m):
    t = m.group(0)
    if FORMULA.search(t) and "\\mathrm{do}" not in t:
        return " " * len(t)
    t = re.sub(r"[_^]\\[A-Za-z]+", _blank, t)
    t = re.sub(r"[_^]\{[^{}]*\}", _blank, t)
    t = re.sub(r"[_^]\w", _blank, t)
    t = re.sub(r"\\mathrm\{[^{}]*\}", _blank, t)
    return t


def mask(s):
    s = re.sub(r"\$[^$]*\$", _math, s)
    s = re.sub(r"(?<![\\\w])v\d+(?:\.\d+)+|(?<![\\\w])[A-Za-z]+\d[\w-]*", _blank, s)
    s = re.sub(r"doi:\S+", _blank, s)
    s = re.sub(r"\b\d{1,2} (?:%s) (?:19|20)\d\d\b" % MONTHS, _blank, s)
    s = re.sub(r"\b(?:%s) (?:19|20)\d\d\b" % MONTHS, _blank, s)
    s = re.sub(r"(?<![\d.,])(?:19|20)\d\d(?![\d.,%])", _blank, s)
    return s


def macros_used(s):
    s = strip_comments(s)
    if "\\begin{document}" in s:
        s = s.split("\\begin{document}", 1)[1]
    return set(re.findall(r"\\([A-Za-z]+)", s))
