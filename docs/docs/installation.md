# Installing LypidDyn

LypidDyn is written in Python, specifically Python 3.6.
The following instruction will set up a virtual environment with copy of the pipeline up and running on your machine.


## Requirements

Softwares and libraries required to work :

    python3.6
    Gromacs5.X
    FatSlim
    GromacsWrapper
    MDAnalysis

# Installing

Step by step guide to setup the virtual environment for Python 3.6:

On your terminal run :

    virtualenv "LypidDyn" -p /usr/bin/python3.6

Activate the virtual environment and install all the packages into the environment

    source LypidDyn/bin/activate

Once the environment is activated we can continue with the different libraries using pip.

    pip install fatslim 
    pip install GromacsWrapper
    pip install MDAnalysis

To deactivate the virtual environment:

    deactivate



