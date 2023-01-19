#!/usr/bin/env python3
'''
Shared functions for ECmean4. Set of formula based tools
'''

#####################
# FORMULA FUNCTIONS #
#####################

import operator
import re
from ecmean.libs.general import is_number


# this is a tool to parse a CDO-based formula into mathematical operatos
# there might exists something more intelligent such as the pyparsing package

def eval_formula(mystring, xdataset):
    """Evaluate the cmd string provided by the yaml file
    producing a parsing for the derived variables"""

    # Tokenize the original string
    token = [i for i in re.split('(\\W+)', mystring) if i]
    if len(token) > 1:
        # Use order of operations
        out = _operation(token, xdataset)
    else:
        out = xdataset[token[0]]
    return out


def _operation(token, xdataset):
    """Parsing of the CDO-based commands using operator package
    and an ad-hoc dictionary. Could be improved, working with four basic
    operations only."""

    # define math operators: order is important, since defines
    # which operation is done at first!
    ops = {
        '/': operator.truediv,
        "*": operator.mul,
        "-": operator.sub,
        "+": operator.add
    }

    # use a dictionary to store xarray field and call them easily
    dct = {}
    for k in token:
        if k not in ops:
            if not is_number(k):
                dct[k] = xdataset[k]
            else:
                dct[k] = float(k)

    # apply operators to all occurrences, from top priority
    # so far this is not parsing parenthesis
    code = 0
    for p in ops:
        while p in token:
            code += 1
            # print(token)
            x = token.index(p)
            name = 'op' + str(code)
            replacer = ops.get(p)(dct[token[x - 1]], dct[token[x + 1]])
            dct[name] = replacer
            token[x - 1] = name
            del token[x:x + 2]
    return replacer
