/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/langevin/LangevinRun.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef QCD_LANGEVIN_RUN_H
#define QCD_LANGEVIN_RUN_H

namespace Grid {
namespace QCD {
namespace QCDpt {

class LangevinParams {
public:
    double tau;
    double alpha;
    double gfprecision;
    
    int sweeps;
    int save_every;
    
    bool rk;
    
    LangevinParams(int argc, char **argv) {
    
        std::string arg;
        
        if( GridCmdOptionExists(argv,argv+argc,"--tau") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--tau");
            std::stringstream ss(arg);
            ss>>tau;
        } else{
            std::cout << GridLogError << "Use --tau to set the time step for the Langevin process " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--alpha") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--alpha");
            std::stringstream ss(arg);
            ss>>alpha;
        } else alpha = -0.5*tau;
        
        if( GridCmdOptionExists(argv,argv+argc,"--sweeps") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--sweeps");
            std::stringstream ss(arg);
            ss>>sweeps;
        } else{
            std::cout << GridLogError << "Use --sweeps to set the number of sweeps to be performed " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--save-every") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--save-every");
            std::stringstream ss(arg);
            ss>>save_every;
        } else save_every = 0;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-gf") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--enable-gf");
            std::stringstream ss(arg);
            ss>>gfprecision;
        } else gfprecision = 0.;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-rk") ){
            rk = true;
        } else rk = false;
    }
    
};

template <class gimpl>
class LangevinRun {

private:
    
    LangevinParams Params;
    
    GridBase* grid;
    GridParallelRNG *pRNG;
    
    CheckpointerParameters CPparams;
    BinaryHmcCheckpointer<gimpl> *CP;
    
protected:
    
    PertLangevin<gimpl> L;
    
public:
    //ActionSet<typename Gimpl::GaugeField,PNoHirep> TheActions;
    
    
    
    LangevinRun(GridCartesian* grid_, GridParallelRNG *pRNG_, LangevinParams Params_):
    Params(Params_),
    L(grid_,pRNG_,Params_.tau,Params_.alpha)/*
    pRNG(pRNG_)*/
    {
        
        CPparams.config_prefix = "NSPTckpoint_lat";
        CPparams.rng_prefix = "NSPTckpoint_rng";
        CPparams.saveInterval = Params.save_every;
        CPparams.format = "IEEE64BIG";
        
        CP = new BinaryHmcCheckpointer<gimpl>(CPparams);
        
    };
    
    template<class Action>
    void push_back(Action A){
        L.TheActions.push_back(A);
    }
    
    
};

    
}
}
}
#endif
