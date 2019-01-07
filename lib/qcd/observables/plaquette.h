/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef HMC_PLAQUETTE_H
#define HMC_PLAQUETTE_H

namespace Grid {
namespace QCD {

// this is only defined for a gauge theory
template <class Impl>
class PlaquetteLogger : public HmcObservable<typename Impl::Field> {
 public:
  // here forces the Impl to be of gauge fields
  // if not the compiler will complain
  INHERIT_GIMPL_TYPES(Impl);

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    RealD plaq = WilsonLoops<Impl>::avgPlaquette(U);

    int def_prec = std::cout.precision();

    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "Plaquette: [ " << traj << " ] "<< plaq << std::endl;

    std::cout.precision(def_prec);

  }
};

struct WilsonLoopParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonLoopParameters,
                                    int, interval,
                                    bool, do_smearing,
                                    std::vector<std::vector<int>>, sizes);
    
    WilsonLoopParameters(std::vector<std::vector<int>> sizes = {},
                         int interval = 1,
                         bool do_smearing = false):
    sizes(sizes), interval(interval), do_smearing(do_smearing){}
    
    template <class ReaderClass >
    WilsonLoopParameters(Reader<ReaderClass>& Reader){
        read(Reader, "WilsonLoopMeasurement", *this);
    }
    
    void AddLoopTxN(std::vector<int> sizes_){
        sizes.push_back(sizes_);
    }
};

// this is only defined for a gauge theory
template <class Impl>
class WilsonLoopLogger : public HmcObservable<typename Impl::Field> {
    WilsonLoopParameters Pars;
    
public:
    // here forces the Impl to be of gauge fields
    // if not the compiler will complain
    INHERIT_GIMPL_TYPES(Impl);
    
    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;
    
    WilsonLoopLogger(std::vector<std::vector<int>> sizes, int interval = 1, bool do_smearing = false):
    Pars(sizes, interval, do_smearing){}
    
    WilsonLoopLogger(WilsonLoopParameters P):Pars(P){
        std::cout << GridLogDebug << "Creating WilsonLoop " << std::endl;
    }
    
    void TrajectoryComplete(int traj,
                            Field &U,
                            GridSerialRNG &sRNG,
                            GridParallelRNG &pRNG) {
        
        if (traj%Pars.interval == 0){
            
            Field Usmear = U;
            int def_prec = std::cout.precision();
            
            // set smear parameters: number of steps, values of weights...
            if (Pars.do_smearing){
                Field Utmp = U;
                Smear_APE<Impl> APEsmearing;
                for(int i=0; i<10; i++){
                    APEsmearing.smear(Utmp, Usmear);
                    Usmear = ProjectOnGroup(Usmear+Utmp);
                }
            }
            
            for(int i=0; i<Pars.sizes.size(); i++){
                Real plaq = WilsonLoops<Impl>::avgWilsonRectLoop(Usmear,Pars.sizes[i][0],Pars.sizes[i][1]);
                
                if(Pars.do_smearing) std::cout << GridLogMessage << "YesSmearing ";
                else std::cout << GridLogMessage << "NoSmearing ";
                std::cout << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                << "WilsonLoop " << Pars.sizes[i][0] << "x" << Pars.sizes[i][1] << ": [ " << traj << " ] " << plaq << std::endl;
                std::cout.precision(def_prec);
            }
        }
    }
};
    
}  // namespace QCD
}  // namespace Grid

#endif  // HMC_PLAQUETTE_H
