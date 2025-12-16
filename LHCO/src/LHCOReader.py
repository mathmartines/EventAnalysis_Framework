"""Functions to read the content of .lhco files"""

from typing import List, Dict
from EventAnalysis_Framework.LHCO.src.EventInfo import Event
from EventAnalysis_Framework.LHE.src.read_lhe import read_lhe


def read_LHCO(filename: str) -> List:
    """Yields a single event at time."""
    # Holds all the events
    with open(filename) as lhco_file:
        event_particles = []

        # Searches the event information
        for line in lhco_file:
            # Strip whitespace and skip comments
            current_line = line.strip()
            if current_line.startswith("#"):
                continue

            # Signal a new event
            if current_line.startswith("0"):
                if event_particles:
                    yield Event.from_str_particles_info(event_particles)
                # Reset event for the next particles
                event_particles = []

            else:
                # Remove the first char - info not needed
                event_particles.append(current_line[1:])

        # Add last event if it exists
        if event_particles:
            yield Event.from_str_particles_info(event_particles)


def read_LHCO_all_events(filaname: str):
    """Returns a list with all the events."""
    return [event for event in read_LHCO(filaname)]


def read_LHCO_with_weight(filenames: Dict[str, str]):
    """Reads the events from the .lhco file and the weights from the .lhe files"""
    # Reads the lhe file
    lhe_events = read_lhe(filename=filenames["LHE"])

    # Reads the lhco file
    for event_lhco, event_lhe in zip(read_LHCO(filename=filenames["LHCO"]), lhe_events):
        event_lhco.weights = event_lhe.eventinfo.weight
        yield event_lhco
