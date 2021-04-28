#!/usr/bin/env python


# Copyright (C) 2020, Alessia Campo <alessia.campo@studio.unibo.it>

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



def hg_selection(u, dict_selection):

    ''' Parameters:
            u: universe object

            dict_selection: dictionary with lipid types

                            associated to atom heagroups
    '''

    # initialize an empty list to store the atom selection

    selection=list()

    # create a set to store the lipid residue names from the system

    system_set=set(u.atoms.resnames)

    # create a set to store the lipid residue names from the dict

    keys_set=set(dict_selection.keys())

    intersection=system_set.intersection(keys_set)   # intersection between the two sets

    for res in intersection:

        # perform selections

        selection.append(u.select_atoms(f"resname {res} \
                        and name {dict_selection[res]}"))




    return sum(selection)
