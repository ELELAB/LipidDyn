import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import errno

dict_selection={"PSM": "P","POPC":"P", "ERG":"O3" ,"DOPC": "P", "CHL1": "O3"}

u=mda.Universe("heterogeneous_membrane.gro", "heterogeneous_membrane.xtc")


def hg_selection(u, dict_selection):

    '''Parameters:
            u: universe object

            dict_selection: dictionary with lipid type associated

                            to atom heagroups
    '''

    # initialize an empty list to store the atom selection

    selection=list()

    # create a list to parse lipid type names

    lipid_type=[res for res in u.atoms.resnames]

    for key in dict_selection.keys():

        # check the presence of atom names in the dictionary
        if key in lipid_type:

            selection+=list(u.select_atoms(f"resname {key} \
            and name {dict_selection[key]}"))
    if not selection:
        raise UserWarning("No selection has been made")

    return selection

print(hg_selection(u,dict_selection))
