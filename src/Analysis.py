"""
    Framework to perform the single event analysis and the iteration over the events.
"""

from typing import List, Callable
from EventAnalysis_Framework.src.Histogram import Histogram
import copy


class EventAnalysis:
    """
    Performs the analysis of a single event.
    Holds information about particle selections and event selection cuts.
    """

    def __init__(self, cuts: List[Callable], particles_selection=None):
        """
        :param particles_selection:
            Returns an event with a list of particles selected for the analysis
        :param cuts:
            List of functions that represent the selection cuts.
            Each function must return True if the event passes the cut and False otherwise.
        """
        self._particles_selections = particles_selection
        self._cuts = cuts

    def launch_analysis(self, event):
        """
        Launches the analysis on the event.
        Returns a tuple, where the first item is a boolean indicating whether the event should be selected,
        and the second is the modified event after all the particle selections.
        This is done in case the modified event is needed by the EventLoop object for histogram booking.
        """
        # Creates an event with only the particles for the analysis.
        if self._particles_selections is not None:
            event = self._particles_selections(event)
        # Applies the event selection cuts
        passed_cuts = all(cut(event) for cut in self._cuts)
        # Returns the boolean and the modified event
        return passed_cuts, event


class EventLoop:
    """
    Iterates over all events in an .lhe file
    and manages histogram booking with the selected events.
    """

    def __init__(self, file_reader: Callable, histogram: Histogram):
        # Function responsible for reading events
        self._file_reader = file_reader
        # Template of the histogram that should be build for each analysis
        self._histogram_template = histogram

    def analyse_events(self, filename: str, event_analysis: EventAnalysis):
        """
        Runs the analysis on events from the .lhe file and returns a histogram
        constructed from the selected events.

        :param filename: Path to the .lhe file storing the events.
        :param event_analysis: performs the analysis of a single event.

        :return: Dict with the booked histogram for each analysis.
        """
        print(f"Reading events from file: {filename}")

        # Count the number of processed events
        evt_number = 0

        # Generates the histogram for the current analysis
        analysis_hist = copy.copy(self._histogram_template)

        # Iterate over events in the file
        for event in self._file_reader(filename):
            if evt_number > 0 and evt_number % 1000 == 0:
                print(f"INFO: Processed {evt_number} events")

            # Runs the analysis on the current event
            select_event, modified_event = event_analysis.launch_analysis(event=event)

            # Updates the histogram generated for the analysis
            if select_event:
                analysis_hist.update_hist(modified_event)

            evt_number += 1

        # Returns the histogram created for the analysis
        return analysis_hist, evt_number
