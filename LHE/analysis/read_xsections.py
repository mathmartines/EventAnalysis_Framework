"""Reads the cross-sections from a list of .lhe files"""

import numpy as np
import json
from EventAnalysis_Framework.src.Utilities import read_xsection

if __name__ == "__main__":
    # Folder where the files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/HighPT-VH/PartonCrossSection_XCheck/ZH/ddbar/lhe_files"

    # CM energies
    cm_energies = np.arange(300, 10100, 100)

    # List of EFT terms
    eft_terms = ["SM"]

    # Stores the cross-sections (pb)
    xsections = {term: np.empty(len(cm_energies)) for term in eft_terms}

    # Reads the x-section from the files
    for eft_term in eft_terms:
        for index, cm_energy in enumerate(cm_energies):
            lhe_file_name = f"{folderpath}/{eft_term}-sqrt_s_{cm_energy}.lhe"
            xsections[eft_term][index] = read_xsection(path_to_file=lhe_file_name)

    with open(f"{folderpath}/xsections.json", "w") as file_:
        simuations = {term: dist.tolist() for term, dist in xsections.items()}
        json.dump(simuations, file_, indent=4)