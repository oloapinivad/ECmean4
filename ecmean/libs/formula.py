"""
Shared functions for ecmean. Set of formula based tools
"""

#####################
# FORMULA FUNCTIONS #
#####################

import logging
import operator
import re

loggy = logging.getLogger(__name__)


def formula_wrapper(var, face, xfield):
    """
    Wrapper function to extract data-array from a dataset (xfield)
    and to apply the mathematical operation required if it is a derived variable
    """

    if "derived" in face["variables"][var]:
        cmd = face["variables"][var]["derived"]
        return _eval_formula(cmd, xfield)
    return xfield[var]


# Tokenizer and evaluator for the CDO-like formulas of the interface files.
# Supported syntax: variable names, numbers (decimals, scientific notation),
# the four basic operations, unary minus, parentheses and whitespace.

_TOKEN_REGEX = re.compile(
    r"""\s*(?:
        (?P<number>(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)
        |(?P<name>[A-Za-z_][A-Za-z0-9_]*)
        |(?P<op>[-+*/()])
    )""",
    re.VERBOSE,
)


def tokenize_formula(mystring):
    """Split a formula into a list of (kind, text) tokens, kind being
    'number', 'name' or 'op'. Raises ValueError for unsupported syntax."""

    tokens = []
    pos = 0
    mystring = mystring.rstrip()
    while pos < len(mystring):
        match = _TOKEN_REGEX.match(mystring, pos)
        if not match:
            raise ValueError(f"Unsupported syntax in formula '{mystring}' at position {pos}: '{mystring[pos:]}'")
        kind = match.lastgroup
        tokens.append((kind, match.group(kind)))
        pos = match.end()
    if not tokens:
        raise ValueError("Empty formula")
    return tokens


def formula_variables(mystring):
    """Return the unique variable names required by a formula, in order of appearance."""

    return list(dict.fromkeys(text for kind, text in tokenize_formula(mystring) if kind == "name"))


_OPERATIONS = {"/": operator.truediv, "*": operator.mul, "-": operator.sub, "+": operator.add}


def _eval_formula(mystring, xdataset):
    """Evaluate the cmd string provided by the yaml file
    producing a parsing for the derived variables.
    Grammar: expr = term {(+|-) term}; term = factor {(*|/) factor};
    factor = number | name | (expr) | -factor"""

    tokens = tokenize_formula(mystring)
    pos = 0

    def peek():
        return tokens[pos][1] if pos < len(tokens) and tokens[pos][0] == "op" else None

    def advance():
        nonlocal pos
        pos += 1
        return tokens[pos - 1]

    def expr():
        value = term()
        while peek() in ("+", "-"):
            value = _OPERATIONS[advance()[1]](value, term())
        return value

    def term():
        value = factor()
        while peek() in ("*", "/"):
            value = _OPERATIONS[advance()[1]](value, factor())
        return value

    def factor():
        if pos >= len(tokens):
            raise ValueError(f"Unexpected end of formula '{mystring}'")
        kind, text = advance()
        if kind == "number":
            return float(text)
        if kind == "name":
            return xdataset[text]
        if text == "-":
            return -factor()
        if text == "(":
            value = expr()
            if peek() != ")":
                raise ValueError(f"Missing closing parenthesis in formula '{mystring}'")
            advance()
            return value
        raise ValueError(f"Unexpected '{text}' in formula '{mystring}'")

    result = expr()
    if pos != len(tokens):
        raise ValueError(f"Unexpected '{tokens[pos][1]}' in formula '{mystring}'")
    return result
