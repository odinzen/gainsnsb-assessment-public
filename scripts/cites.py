"""Main-text reference numbers for figure legends, read from results/citation_numbers.json."""
import json
import os

_NUM = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                                   "citation_numbers.json")))


def cite(*keys):
    return "[" + ",".join(str(n) for n in sorted(_NUM[k] for k in keys)) + "]"
