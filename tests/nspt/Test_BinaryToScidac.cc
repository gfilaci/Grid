/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_runstaggered.cc

Copyright (C) 2015-2017

Author: Gianluca Filaci <g.filaci@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace QCD;
using namespace QCDpt;

/*
 Filenames must have the structure
 
 configname.number
 rngname.number
 
 and outputs are
 
 scidac_configname.number
 scidac_rngname.number
 
 The lattice size must be given with --grid
 and Grid must be compiled with the correct Np
 */
/***************************************/
typedef TwistedGimpl_ptR gimpl;
string configname = "nameofthecfg";
string rngname    = "nameofrnd";
int number = 0;
/***************************************/

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);
    
    GridParallelRNG pRNG(&Grid);
    GridSerialRNG   sRNG;
    
    typename gimpl::GaugeField U(&Grid);
    
    CheckpointerParameters BCPparams;
    CheckpointerParameters SCPparams;
    
    BCPparams.config_prefix = configname;
    BCPparams.rng_prefix = rngname;
    BCPparams.saveInterval = 1;
    BCPparams.format = "IEEE64BIG";
    
    SCPparams.config_prefix = "scidac_" + configname;
    SCPparams.rng_prefix = "scidac_" + rngname;
    SCPparams.saveInterval = 1;
    SCPparams.format = "IEEE64BIG";
    
    BinaryHmcCheckpointer<gimpl> BCP(BCPparams);
    ScidacHmcCheckpointer<gimpl,ActionParameters> SCP(SCPparams);
    
    BCP.CheckpointRestore(number,U,sRNG,pRNG);
    SCP.TrajectoryComplete(number,U,sRNG,pRNG);
    
    Grid_finalize();
    return EXIT_SUCCESS;
}
