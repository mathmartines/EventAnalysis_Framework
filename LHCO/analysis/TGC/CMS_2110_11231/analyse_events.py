"""Analysis for the WW production channel for the CMS paper arXiv: 2009.00119"""


from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHCO.src.LHCOReader import read_LHCO
from EventAnalysis_Framework.LHCO.src.Observables import InvariantMass
from EventAnalysis_Framework.LHCO.analysis.TGC.CMS_2110_11231 import selection_cuts
import copy
import json


def construct_eft_term_nams(term: str, indices: str):
    """Constructs the name of the eft term"""
    return "-".join([f"{coef}{indices}" for coef in term.split("-")])


if __name__ == "__main__":
    # Simulated terms
    eft_terms = [
        "SM"
    ]

    # Path where the simulations are store
    folderpath = "/home/martines/work/MG5_aMC_v3_1_1/PhD/TGC/WZ/CMS_2110_11231"

    # Event analysis with all the selection cuts
    event_analysis = EventAnalysis(
        cuts=[
            selection_cuts.number_of_leptons,
            # selection_cuts.check_ossf_pair,
            # selection_cuts.missing_energy_cut,
            # selection_cuts.leptons_cut,
            # selection_cuts.btag_veto,
            # selection_cuts.min_invariant_mass,
            # selection_cuts.invariant_mass_all_leptons
        ]
    )

    # Books the histogram for the analysis
    bin_edges = [100, 200, 300, 400, 500, 1000, 1500, 3000]
    mwz_hist = ObservableHistogram(
        bin_edges=bin_edges,
        observable=InvariantMass(particles=["electrons", "muons", "met"])
    )

    # Constructs the EventLoop object
    event_loop = EventLoop(file_reader=read_LHCO, histogram=mwz_hist)

    # Book the histogram for the signal
    histograms_efts = {
        eft_term: copy.copy(mwz_hist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 7):
            # name of the .lhe file
            file_name_lhe = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            file_name_lhco = f"{folderpath}/lhco_files/{eft_term}-bin-{bin_index}.lhco"
            # Cross-section
            xsection = read_xsection(file_name_lhe)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name_lhco, event_analysis)
            # Updates the histogram for the current term
            histograms_efts[eft_term] += (xsection * 137 * 1000 / number_of_evts) * current_hist
            print(current_hist, xsection)

    print(histograms_efts)

    with open(f"{folderpath}/CMS_2110_11231.json", "w") as file_:
        simulations = {eft_term: dist.tolist() for eft_term, dist in histograms_efts.items()}
        json.dump(simulations, file_, indent=4)

