"""Wrapper to read the lhe files."""

import pylhe


def read_lhe(filename: str):
    """Returns a generator over all the events in the file."""
    lhe_file = pylhe.read_lhe_file(filepath=filename)
    return lhe_file.events
