Test Set Wrapper
================

Python script to run testsets in bulk with simple input files.


Requirements
------------

- Python 3.5 or higher
- DFTB+ (any version)


For calculations with DTNN:

- ASE >= 3.19
- torch >= 1.5
- SchNetPack >= 0.3


Installation
------------

Get the repository with

    git clone git@github.com:MKubillus/testset_wrapper.git

Or download the zip file (green "Code" button on the top right) and extract the archive with

    unzip testset_wrapper-master.zip


Usage
-----

Run the main executable with

    python3 run_testsets.py

in an empty folder to get a default input file "testsets_config.yml" and edit it to fit your system. If you run the program then again, it will read the input file automatically.

If you want to learn how to use the test set wrapper please check out the example folder for some working examples. Note that you might have to adjust the dftb_in.hsd to your system (Slater-Koster file locations) and add the full path to your DFTB+ executable if it is not in your system path.


License
-------

This software is licensed under the GNU General Public License v3.0 (GPL3). See LICENSE file for more information.
