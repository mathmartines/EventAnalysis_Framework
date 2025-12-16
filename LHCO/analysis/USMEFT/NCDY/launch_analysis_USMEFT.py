"""Launches the analysis on the EFT terms for NC DY channel"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHCO.src.LHCOReader import read_LHCO_with_weight
from EventAnalysis_Framework.LHCO.src.Observables import InvariantMass
from EventAnalysis_Framework.LHCO.src.EventInfo import Event
import copy
import math
import json


def get_weight(event: Event):
    """Returns the weight of the event"""
    return event.weights


def number_of_electrons(event: Event):
    """Cut on the number of electrons."""
    return len(event.electrons) > 1


def leptons_kin_cuts(event: Event):
    """Applies the cuts on the transverse momentum of the final state leptons."""
    # Cuts on the pT and pseudo-rapidity
    leptons_kin = all([lepton.pt > 20 and abs(lepton.eta) < 2.5 for lepton in event.electrons])
    # Veto event if needed
    return leptons_kin


def opposite_charge(event: Event):
    """Checks if the electrons have opposite charge"""
    leptons_trk = [lepton.ntrk for lepton in event.electrons]
    return math.prod(leptons_trk) < 0


if __name__ == "__main__":
    # Path to where the files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/USMEFT-HLLHC/NCDY-d8"

    # All the EFT terms we need to launch the analysis
    eft_terms = [
        "SM",
        # Linear d6 terms
        "Cphi1", "Delta4F", "CBW", "C2JB",
        # d6 squared terms
        "Cphi1-Cphi1", "Cphi1-CBW", "Cphi1-C2JB", "Cphi1-Delta4F",
        "Delta4F-Delta4F", "Delta4F-CBW", "Delta4F-C2JB",
        "CBW-CBW", "CBW-C2JB",
        "C2JB-C2JB",
        # Linear d8 terms
        "C2psi4D2", "C3psi4D2", "C7psi4H2"
    ]

    # Total number of bins we simulated
    n_bins = 10

    # Observable we need to use to build the histogram
    dielectron_mass = InvariantMass(particles=["electrons"])
    bin_edges = [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 1000000]
    mee_hist = ObservableHistogram(
        bin_edges=bin_edges,
        observable=dielectron_mass,
        get_weight=get_weight
    )

    # Applies the cuts on a single event
    event_analysis = EventAnalysis(cuts=[number_of_electrons, leptons_kin_cuts, opposite_charge])

    # Iterates over all the events in the file
    event_loop = EventLoop(file_reader=read_LHCO_with_weight, histogram=mee_hist)

    # Books the histogram for each eft term
    eft_hists = {eft_term: copy.copy(mee_hist) for eft_term in eft_terms}

    # Iterates over all EFT terms
    for eft_term in eft_terms:
        # Iterates over the number of bins
        for bin_index in range(1, n_bins + 1):
            # Prefix of the files
            file_prefix = f"{eft_term}-bin-{bin_index}"

            # files needed to perform the analysis
            input_files = {
                "LHE": f"{folderpath}/lhe_files/{file_prefix}.lhe",
                "LHCO": f"{folderpath}/lhco_files/{file_prefix}.lhco"
            }

            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(input_files, event_analysis)
            print(current_hist, number_of_evts)

            # Updates the histogram for the current term (in fb)
            eft_hists[eft_term] += 1e3 * current_hist / number_of_evts

        print(eft_hists[eft_term])

    # Saves the json file
    with open(f"{folderpath}/NCDY.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)
