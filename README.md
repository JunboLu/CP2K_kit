# CP2K_kit - A helper of CP2K
---------

CP2K_kit is a code to:  
(1) Generate CP2K input file  
(2) Analyze CP2K trajectory  
(3) Optimize GTH pseudopotentials and basis sets  
(4) Do thermodynamic integration  
(5) Optimize neural network force field  

The code is mainly written in Python and Fortran.  

Email: lujunbo15@gmail.com
  
# Installation for CP2K_kit

* Prerequisites
   - Python 3.5 or higher
   - Numpy 1.8.0 or higher
   - psutil

   Suggestion: as we will use deepmd-kit software, we recommend users could install deepmd-kit  
   at first. Then add the environmental variable of deepmd-kit. It will include python and numpy.  
   Then intall psutil by using pip3 in deepmd-kit/bin directory:  
   pip3 install psutil  

* Install GNU parallel

    Download GNU parallel source code from https://www.gnu.org/software/parallel/  
    tar -jxvf parallel-latest.tar.bz2  
    cd parallel-latest  
    ./configure --prefix=parallel_install_path  
    make  
    make install  

* Compile core module
  
    git clone https://github.com/JunboLu/CP2K_kit.git  
    !Caution: When downloading CP2K_kit through zip version, please change the name "CP2K_kit-main"  
    to "CP2K_kit" after you unzip it.  
    cd CP2K_kit_directory/lib  
    change directory of f2py in Makefile  
    !Caution: When the gcc version is low, f2py cannot compile core code successfully, please update your  
    gcc up to 6.3.  
    make  

* Environmental variable

    export PYTHONPATH=../CP2K_kit_directory:$PYTHONPATH  
    change python_exe and CP2K_kit directory in CP2K_kit_directory/bin/CP2K_kit file  
    export PATH=CP2K_kit_directory/bin:$PATH  

# How to use 
* CP2K_kit is an user-friendly code.  

  The input files are in example directory.  
  Users just need to run:  
  CP2K_kit input.inp  
