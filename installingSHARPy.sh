#!/bin/bash
# Installing SHARPy onto a fresh Ubuntu 16.04
# r.simpson11@imperial.ac.uk
# 15/08/16
# Prerequisites: eclipse neon with CDT version 9.x, pydev verison 5.1.2

# First, get all the software into your local repo with a one-liner from GitHub
cd Software/git/
git clone https://github.com/RJSSimpson/SHARPy.git

# We're going to make the beam libraries and use f2py to auto-genrate a lot of beam
# functions for Python (3).
# This will require compiling the fortran code (using gfortran), and using the Python 
# headers that come in the python3-dev package.
sudo apt-get install gfortran
sudo apt-get install python3-f2py
sudo apt-get install python3-dev
cd BeamLib/src/wrapper/
make
cd f2py/
./runf2py.sh

# Now we're going to open the UVLMLib as a c/c++ project in eclipse and compile it.
# As stated in the README we will need Eigen and Boost C++ libraries to do so.
cd ~/Software/Downloads/
sudo apt install mercurial
hg clone https://bitbucket.org/eigen/eigen/
wget -O boost_1_61_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz/download
tar zxvf boost_1_61_0.tar.gz

# Move to UVLMLib directory, start eclipse and open project.
cd ../../../../UVLMLib
~/eclipse/cpp-neon/eclipse/./eclipse &
echo "INSTRUCTIONS for opening UVLMLib: Whilst in the C/C++ developer perspective go to File->Open Projects from File System, select the current directory (UVLMLib) and eclipse should find a project. As it is linked to the SHARPy project (not yet opened) it may ask you to configure the Python interpreter - choose the advanced autoconfig and select the Python 3 version you will use and leave all the other boxes checked.\n"

# Configure eclipse to include Eigen and boost in the compilation
echo "INSTRUCTIONS for compiling UVLMLib.so: Right-click on project->Preferences, then open C/C++ Build from side menu -> Settings, then open GCC C++ compiler - > Includes, and edit or add the paths to the 'Eigen' and 'boost_1_61_0' directories that we just downloaded. After confirming these changes right-click on the project -> Build Project.\n"

# For good measure at this point we can also ass the BeamLib as a project in eclipse.
echo "INSTRUCTIONS for opening BeamLib: Whilst in the C/C++ developer perspective go to File->Open Projects from File System, select the BeamLib directory (SHARPy/BeamLib) and eclipse should find a project. We compiled this above using make -- if you want to test the project configration then right click on the project -> Clean Project, then right click again -> Build Project.\n"

# install python packages required for SHARPy
sudo apt-get install python3-numpy
sudo apt-get install python3-scipy
sudo apt-get install python3-matplotlib

# Open SHARPy and modify the path to reflect your system.
echo "INSTRUCTIONS for opening SHARPy: Whilst in PyDev perspective go to File->Open Projects from File System, select the SHARPy directory and eclipse should find a project. Uncheck all the auto-discover options so eclipse imports this as a project (not a c/c++ project). Modify [SharPyProjectDir]/src/Main/SharPySettings.py so that SharPyProjectDir = [SharPyProjectDir] on your system.\n"

# finished
echo "FINISHED: Now run [SharPyProjectDir]/src/Main/PyCoupled/test/test_Coupled_NlnStatic.py to test the structure, aerodynamics and interface for nonlinear static problems."
