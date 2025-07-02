"""
parse_config.py

Read in simulation parameters from a C file with constants and basic expressions

2025.07.01
Anthony Atkinson
"""

import re

class Params(dict):
    
    def __getattr__(self, lhs):
        return self[lhs]

    def __setattr__(self, lhs, rhsue):
        self[lhs] = rhsue

def parse_config(cfname: str, debug: int=0):
    params = Params(cfname=cfname)
    empty_map = {}

    # pattern for C variable assignment
    passign = re.compile(r"^[^/]+\s([a-zA-Z_]\w*)\s*=\s*(.+);$")
    # pattern for C string assignment
    pstring = re.compile(r"^[^/]+char \*([a-zA-Z_]\w*)\s*=\s*(['\"].+['\"]);$")
    # pattern for C define pre-processing directive for non-macros
    pdefine = re.compile(r"^#define\s+([a-zA-Z_]\w*)\s+(\S.*)$")

    # open config file for reading
    with open(cfname, "r") as f:
        # parse each line for simulation parameters
        for line in f:
            
            if debug >= 2:
                print(line)

            # TODO make this nicer or something
            # match C variable assignment
            massign = re.match(passign, line)
            if massign:
                lhs = massign[1]
                rhs = massign[2]
            else:
                # match C string assignment
                mstring = re.match(pstring, line)
                if mstring:
                    lhs = mstring[1]
                    rhs = mstring[2]
                else:
                    # match C define pre-processing directive
                    mdefine = re.match(pdefine, line)
                    if mdefine:
                        lhs = mdefine[1]
                        rhs = mdefine[2]

            # line matches syntax for assignment or define
            if massign or mstring or mdefine: 
                s = f"{lhs} = {rhs}"
                # assign values locally to evaluate rhs expressions
                exec(s, empty_map, params)
                val = eval(rhs, empty_map, params)
                
                # store param and value in dictionary
                params[lhs] = val

                if debug >= 1:
                    print(s, "-->", val)

    if debug >= 2:
        print(params)

    return params


