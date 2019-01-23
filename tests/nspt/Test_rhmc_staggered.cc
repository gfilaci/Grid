    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_rhmc_staggered.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

namespace Grid {
    class GaugeActionParameters : Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeActionParameters,
                                        double, beta);
        template <class ReaderClass >
        GaugeActionParameters(Reader<ReaderClass>& Reader){
            read(Reader, "GaugeAction", *this);
        }
    };
    class FermionActionParameters : Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FermionActionParameters,
                                        double, mass);
        template <class ReaderClass >
        FermionActionParameters(Reader<ReaderClass>& Reader){
            read(Reader, "FermionAction", *this);
        }
    };
}

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  typedef Grid::JSONReader       Serialiser;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  typedef QCDpt::NaiveStaggeredImplD FermionImplPolicy;
  typedef StaggeredFermion<QCDpt::NaiveStaggeredImplD> FermionAction;
  typedef typename FermionAction::FermionField FermionField;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCWrapper TheHMC;

  TheHMC.ReadCommandLine(argc, argv);
  if (TheHMC.ParameterFile.empty()){
      std::cout << "Input file not specified. Use --ParameterFile option in the command line.\nAborting" << std::endl;
      exit(1);
  }

  Serialiser Reader(TheHMC.ParameterFile);

  TheHMC.Resources.AddFourDimGrid("gauge");

  CheckpointerParameters CPparams(Reader);
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  typedef WilsonLoopMod<HMCWrapper::ImplPolicy> WLoop;
  WilsonLoopParameters WLParams(Reader);
  TheHMC.Resources.AddObservable<WLoop>(WLParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  GaugeActionParameters GPar(Reader);
  WilsonGaugeActionR Waction(GPar.beta);

  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  LatticeGaugeField U(GridPtr);
  FermionActionParameters FPar(Reader);
  FermionAction FermOp(U, *GridPtr, *GridRBPtr, FPar.mass);

  // 2 staggered flavours (parameters taken from the one-flavour Wilson fermion action)
  OneFlavourRationalParams Params(1.0e-4, 64.0, 2000, 1.0e-6);
  TwoFlavourRationalStaggeredPseudoFermionAction<FermionImplPolicy> WilsonNf2(FermOp,Params);
  WilsonNf2.is_smeared = false;

  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&WilsonNf2);
  ActionLevel<HMCWrapper::Field> Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.initialize(Reader);
  TheHMC.Run();

  Grid_finalize();

} // main
