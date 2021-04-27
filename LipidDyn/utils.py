#!/usr/bin/env python


# Copyright (C) 2019, Simone Scrima <simonescrima@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.




dict_selection={"PSM": "P","POPC":"P", "ERG":"O3" ,"DOPC": "P", "CHL1": "O3"}


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

    return selection

print(hg_selection(u,dict_selection))
