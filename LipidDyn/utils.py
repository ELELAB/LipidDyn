#!/usr/bin/env python

# Copyright (C) 2021, Alessia Campo <alessia.campo@studio.unibo.it>,

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



def select_lipid_headgroups(u, dict_selection):

    '''
        Parameters:
            u: universe object
                Universe on which the selection is performed
            dict_selection: dictionary lipid:atom type
                For each lipid type, the atom name of its headgroup
        Returns:
            selections of all the headgroup atoms of all the
    '''

    # initialize an empty list to store the atom selection
    selection=list()

    # create a set to store the lipid residue names from the system
    system_set=set(u.atoms.resnames)

    # create a set to store the lipid residue names from the dict
    keys_set=set(dict_selection.keys())

    # intersection between the two sets
    intersection=system_set.intersection(keys_set)


    # select headgroup atom for each identified type
    for res in intersection:
        selection.append(u.select_atoms(f"resname {res} \
                        and name {dict_selection[res]['headgroup']}"))

    # return single selection including all the atoms
    return sum(selection)
