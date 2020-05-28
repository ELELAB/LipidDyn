#!/usr/bin/env python

# conftest.py - set up common fixtures for tests
# Copyright (C) 2020 Matteo Tiberti <matteo.tiberti@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import pytest
import MDAnalysis as mda


@pytest.fixture
def root_dir():
    return os.path.dirname(os.path.abspath(__file__))

@pytest.fixture
def data_dir(root_dir):
    return "{0}/data".format(root_dir)

@pytest.fixture
def membrane_sim():
    return "heterogeneous_membrane"

@pytest.fixture
def membrane_prot_sim():
    return "heterogeneous_membrane_protein"

@pytest.fixture
def membrane_uni(data_dir, membrane_sim):

    membrane_sim_file_top = "{0}/{1}.gro".format(data_dir, membrane_sim)
    membrane_sim_file_traj = "{0}/{1}.xtc".format(data_dir, membrane_sim)

    return mda.Universe(membrane_sim_file_top,
                        membrane_sim_file_traj)

@pytest.fixture
def membrane_protein_uni(data_dir, membrane_prot_sim):

    membrane_prot_sim_file_top = "{0}/{1}.gro".format(data_dir, membrane_prot_sim)
    membrane_prot_sim_file_traj = "{0}/{1}.xtc".format(data_dir, membrane_prot_sim)

    return mda.Universe(membrane_prot_sim_file_top,
                        membrane_prot_sim_file_traj)



