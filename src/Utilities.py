

def read_xsection(path_to_file: str, default_line: str = "#  Integrated weight (pb)  :"):
    """Reads the cross-section from a .lhe or banner file"""
    with open(path_to_file) as info_file:
        for line in info_file:
            if line.startswith(default_line):
                _, xsection = line.split(":")
                return float(xsection)
    return None
