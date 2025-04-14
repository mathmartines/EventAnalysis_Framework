"""Analysis for the WW production channel for the CMS paper arXiv: 2009.00119"""


from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHCO.src.LHCOReader import read_LHCO
from EventAnalysis_Framework.LHCO.src.Observables import InvariantMass
from EventAnalysis_Framework.LHCO.analysis.TGC.CMS_2009_00119 import selection_cuts
import copy
import json


def construct_eft_term_nams(term: str, indices: str):
    """Constructs the name of the eft term"""
    return "-".join([f"{coef}{indices}" for coef in term.split("-")])


if __name__ == "__main__":
    # Simulated terms
    eft_terms = [
        "C3Hq", "C3Hq-C3Hq"
    ]
    flavor_indices = [
        "11", "12"
    ]

    # Path where the simulations are store
    folderpath = "/home/martines/work/MG5_aMC_v3_1_1/PhD/TGC/WW/CMS_2009_00119"

    # Event analysis with all the selection cuts
    event_analysis = EventAnalysis(
        particles_selection=selection_cuts.select_objects,
        cuts=[
            selection_cuts.opposite_sign_lepton_pair,
            selection_cuts.leptons_pt_cuts,
            selection_cuts.missing_energy_cut,
            selection_cuts.lepton_pair_cuts,
            selection_cuts.btag_veto,
            selection_cuts.number_of_jets,
            selection_cuts.projected_ptmiss_cut
        ]
    )

    # Books the histogram for the analysis
    bin_edges = [100, 200, 300, 400, 500, 600, 700, 750, 800, 850, 1000, 100000000000000]
    mll_hist = ObservableHistogram(
        bin_edges=bin_edges,
        observable=InvariantMass(particles=["electrons", "muons"])
    )

    # Constructs the EventLoop object
    event_loop = EventLoop(file_reader=read_LHCO, histogram=mll_hist)

    # Book the histogram for the signal
    histograms_efts = {
        construct_eft_term_nams(eft_term, indices): copy.copy(mll_hist)
        for eft_term in eft_terms for indices in flavor_indices
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the flavor indices
        for flavor_index in flavor_indices:
            eft_term_name = construct_eft_term_nams(eft_term, flavor_index)
            # Iterates over the simulated bins
            for bin_index in range(1, 11):
                # name of the .lhe file
                file_name_lhe = f"{folderpath}/lhe_files/{eft_term}-{flavor_index}-bin-{bin_index}.lhe"
                file_name_lhco = f"{folderpath}/lhco_files/{eft_term}-{flavor_index}-bin-{bin_index}.lhco"
                # Cross-section
                xsection = read_xsection(file_name_lhe)
                # Run the analysis on the file
                current_hist, number_of_evts = event_loop.analyse_events(file_name_lhco, event_analysis)
                # Updates the histogram for the current term
                histograms_efts[eft_term_name] += (xsection * 35.9 * 1000 / number_of_evts) * current_hist
                print(current_hist, xsection)

    print(histograms_efts)

    with open(f"{folderpath}/CMS_2009_00119.json", "w") as file_:
        simulations = {eft_term: dist.tolist() for eft_term, dist in histograms_efts.items()}
        json.dump(simulations, file_, indent=4)

