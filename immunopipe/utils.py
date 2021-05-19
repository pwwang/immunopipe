"""Utilities"""
from textwrap import dedent

def text_dedent(text):
    splits = text.split('\n', 1)
    if len(splits) > 1:
        first_line, rest = splits
        return f"{first_line}\n{dedent(rest)}"
    return splits[0]
