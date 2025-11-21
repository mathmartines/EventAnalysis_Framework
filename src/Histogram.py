"""Interface to construct a binning histogram."""

from abc import ABC, abstractmethod
import numpy as np
from typing import List, Callable, Dict
import copy


class Histogram(ABC):
    """Defines the interface for a histogram class."""

    @abstractmethod
    def update_hist(self, event):
        """Updates the histogram with a given event."""
        raise RuntimeError("Trying to use a method from an abstract class.")

    @abstractmethod
    def __copy__(self):
        """Clones an empty histogram."""
        pass


class BinIndexFinder:
    """Finds the index of the bin of the histogram that must be updated."""

    def __init__(self, bin_edges):
        self.bin_edges = bin_edges

    def find_bin_index(self, observable_value: float) -> int:
        """Finds the respective bin index for the given value of the observable."""
        # Look for the right bin number
        for bin_index in range(len(self.bin_edges) - 1):
            if self.bin_edges[bin_index] <= observable_value < self.bin_edges[bin_index + 1]:
                return bin_index
        return -1


def unweighted_events(event):
    """No weights for the events."""
    return 1


class ObservableHistogram(Histogram, np.ndarray, BinIndexFinder):
    """
    One-dimensional histogram for a given observable.
    It behaves like a numpy array.
    The histogram is updated using the method:

    def update_hist(event)

    where it takes a single event as the argument.
    """

    def __new__(cls, bin_edges: List[float], observable: Callable, get_weight: Callable = unweighted_events):
        """
        The two necessary parameters for construction are:

        :param bin_edges: The respective bin edges for the histogram.
        :param observable: A function or callable object that computes
                           the observable for a single Event object.
        """
        # Create an empty histogram
        hist = np.zeros(shape=len(bin_edges) - 1).view(cls)
        # Store the bin_edges and observable as attributes
        hist.bin_edges = bin_edges
        hist.observable = observable
        hist.get_weight = get_weight
        # Return the histogram
        return hist

    def __init__(self, bin_edges: List[float], observable: Callable, get_weight: Callable = unweighted_events):
        BinIndexFinder.__init__(self, bin_edges=bin_edges)

    def __array_finalize__(self, hist):
        if hist is None:
            return
        # Add the attributes
        self.observable = getattr(hist, "observable", None)
        self.bin_edges = getattr(hist, "bin_edges", None)
        self.get_weight = getattr(hist, "get_weight", None)

    def update_hist(self, event):
        """Updates the histogram using the Event object."""
        # Calculate the observable
        obs_value = self.observable(event)
        # Find the bin index
        bin_index = self.find_bin_index(observable_value=obs_value)
        # Get the weight
        weight = self.get_weight(event)
        # Update the histogram if the observable is inside the histogram limits
        if 0 <= bin_index < len(self):
            self[bin_index] += weight

    def __copy__(self):
        """Shallow copy of the current histogram."""
        return self.__new__(self.__class__, bin_edges=self.bin_edges, observable=self.observable,
                            get_weight=self.get_weight)


class WeightedHistogramManager(Histogram, BinIndexFinder):
    """Builds one histogram for each reweighted events."""

    def __init__(self, bin_edges: List[float], observale: Callable, get_weights: Callable, hist_names: List[str]):
        BinIndexFinder.__init__(self, bin_edges=bin_edges)
        # Stores the observable needed
        self.observable = observale
        # Returns a dictionary with the different weights for the event
        self.get_weights_func = get_weights
        # One histogram for which reweighted event
        self._hists = {hist_name: np.zeros(len(bin_edges) - 1) for hist_name in hist_names}

    def update_hist(self, event):
        """Updates all the histograms with the current event."""
        observable_value = self.observable(event)
        # Check if the value is inside the limits of the histogram
        bin_index = self.find_bin_index(observable_value=observable_value)
        if 0 <= bin_index < len(self.bin_edges):
            # Get the weights for each of the histograms
            weights = self.get_weights_func(event)
            # Updates each of the histograms
            for hist_name, hist in self._hists.items():
                # In case there's not weight, set it to 1
                hist[bin_index] += weights[hist_name] if hist_name in weights else 1

    def __getitem__(self, hist_name: str):
        """Returns the histogram"""
        if hist_name in self._hists:
            return self._hists[hist_name]
        return np.zeros(len(self.bin_edges) - 1)

    def __copy__(self):
        """Shallow copy of the current hist."""
        return self.__class__(bin_edges=self.bin_edges, observale=self.observable,
                              get_weights=self.get_weights_func, hist_names=list(self._hists.keys()))


class HistogramCompound(Histogram):
    """Stores a set of Histogram objects that must be updated."""

    def __init__(self, histograms: Dict[str, Histogram]):
        # Stores a dictionary whose values represent histogram, and the key is a name used to identify them
        self._hist_dict = histograms

    def update_hist(self, event):
        """Updates all the histograms with the given event."""
        for hist_name in self._hist_dict:
            self._hist_dict[hist_name].update_hist(event=event)

    def get_hist(self, hist_name: str):
        """Returns the Histogram object associated with the key 'hist_name'"""
        if hist_name in self._hist_dict:
            return self._hist_dict[hist_name]

    def __copy__(self):
        """Returns a shallow clone of all histograms."""
        clone_dict = {hist_name: copy.copy(hist) for hist_name, hist in self._hist_dict.items()}
        # Creates a container with now the cloned hists
        return self.__class__(histograms=clone_dict)
