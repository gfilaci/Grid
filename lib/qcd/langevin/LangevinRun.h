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

// by default, a configuration with one order less
// (in beta) than the current one is usually loaded.
// change if needed
const int load_Nplow = Np - 10;

struct ActionParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ActionParameters,
				    double, beta)

    ActionParameters(): beta(0.) {};

  };

template <class PFermionImpl>
class LangevinParams {
public:
    double tau;
    double alpha;
    double gfprecision;
    
    int sweeps;
    int save_every;
    int pog_every;
    int StartTrajectory;
    
    int Nf;
    PRealD mass;
    std::vector<Complex> boundary_phases;
    std::vector<Complex> prop_phases1;
    std::vector<Complex> prop_phases2;
    std::vector<Complex> prop_phases3;
    
    bool rk;
    bool measureprop;
    bool fagf;
    bool enablenorm;
    bool pog_at_startup;
    std::string StartingType;
    std::vector<int> rngseed;
    
    std::string basename, nameprefix = "";
    CheckpointerParameters CPparams;
    
    StochasticFermionAction<PFermionImpl> *FA;
    
    // in the quenched case, the parameters Nf, mass and boundary_phases can be omitted
    LangevinParams(int argc, char **argv, int Nf_=0, PRealD mass_=zero, std::vector<Complex> boundary_phases_={0.}, StochasticFermionAction<PFermionImpl> *FA_ = NULL):
    FA(FA_)
    {
    
        std::string arg;
        
        Nf = Nf_;
        if(Nf_!=0){
            mass = mass_;
            boundary_phases = boundary_phases_;
            prop_phases1 = boundary_phases_;
            prop_phases2 = boundary_phases_;
            prop_phases3 = boundary_phases_;
            
            Complex im(0.,1.);
            Real phase_exp;
            // reference momenta for L=32, L=24, L=20
            phase_exp = M_PI*((double)GridDefaultLatt()[0])/32.;
            prop_phases1[0] = exp(im*phase_exp);
            phase_exp = M_PI*((double)GridDefaultLatt()[0])/24.;
            prop_phases2[0] = exp(im*phase_exp);
            phase_exp = M_PI*((double)GridDefaultLatt()[0])/20.;
            prop_phases3[0] = exp(im*phase_exp);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--tau") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--tau");
            std::stringstream ss(arg);
            ss>>tau;
        } else{
            std::cout << GridLogError << "Use --tau to set the time step for the Langevin process" << std::endl;
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
            std::cout << GridLogError << "Use --sweeps to set the number of sweeps to be performed" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--save-every") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--save-every");
            std::stringstream ss(arg);
            ss>>save_every;
        } else save_every = 0;
        
        if( GridCmdOptionExists(argv,argv+argc,"--pog-every") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--pog-every");
            std::stringstream ss(arg);
            ss>>pog_every;
        } else pog_every = 0;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-gf") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--enable-gf");
            std::stringstream ss(arg);
            ss>>gfprecision;
        } else gfprecision = 0.;
        
        if( GridCmdOptionExists(argv,argv+argc,"--fourier-acceleration") ){
            fagf = true;
        } else fagf = false;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-rk") ){
            rk = true;
        } else rk = false;
        
        if( GridCmdOptionExists(argv,argv+argc,"--pog-at-startup") ){
            pog_at_startup = true;
        } else pog_at_startup = false;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-propagator") ){
            measureprop = true;
        } else measureprop = false;
        
        if( GridCmdOptionExists(argv,argv+argc,"--enable-norm") ){
            enablenorm = true;
        } else enablenorm = false;
        
        if( GridCmdOptionExists(argv,argv+argc,"--start-from") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--start-from");
            std::stringstream ss(arg);
            ss>>StartTrajectory;
        } else StartTrajectory = -1;
        
        
        if (GridCmdOptionExists(argv,argv+argc,"--StartingType")) {
            arg = GridCmdOptionPayload(argv,argv+argc,"--StartingType");
            if (arg != "ColdStart" && arg != "CheckpointStart" && arg != "LowerOrderStart") {
                std::cout << GridLogError << "Unknown option in --StartingType" << std::endl;
                std::cout << GridLogError << "Valid [ColdStart, CheckpointStart, LowerOrderStart]" << std::endl;
                exit(EXIT_FAILURE);
            }
            StartingType = arg;
        } else{
            std::cout << GridLogError << "Use --StartingType to specify how to initialise the gauge field" << std::endl;
            std::cout << GridLogError << "Valid [ColdStart, CheckpointStart, LowerOrderStart]" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--rng") ){
            arg = GridCmdOptionPayload(argv,argv+argc,"--rng");
            GridCmdOptionIntVector(arg,rngseed);
        } else if(StartingType=="ColdStart"){
            std::cout << GridLogError << "This is a ColdStart: use --rng to initialise the random number generator" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if(StartTrajectory==-1 && (StartingType=="CheckpointStart"||StartingType=="LowerOrderStart")){
            std::cout << GridLogError << "This is a " << StartingType << ": use --start-from to specify which configuration has to be loaded" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if(StartTrajectory!=-1 && StartingType=="ColdStart"){
            std::cout << GridLogError << "This is a ColdStart: the option --start-from must not be used" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if( GridCmdOptionExists(argv,argv+argc,"--prefix") ){
            nameprefix = GridCmdOptionPayload(argv,argv+argc,"--prefix");
        }
        
        basename = initialiseName();
        
        CPparams.config_prefix = nameprefix + basename + "_lat";
        CPparams.rng_prefix = nameprefix + basename + "_rng";
        CPparams.saveInterval = 1;
        CPparams.format = "IEEE64BIG";
        
    }
    
    std::string initialiseName(int Norder = Np){
    
        std::stringstream ss;
        std::string _basename;
        
        _basename = "SU";
        ss << Nc;
        _basename += ss.str() + "_";
        ss.str(std::string());
        
        if(Nf!=0) {
            ss << Nf;
            _basename += "Nf" + ss.str() + "_";
            ss.str(std::string());
        }
        
        ss << Norder;
        _basename += "Np" + ss.str() + "_";
        ss.str(std::string());
        
        for (int d=0; d<Nd; d++) {
            ss << GridDefaultLatt()[d];
            _basename += ss.str();
            if(d!=Nd-1) _basename += "x";
            else _basename += "_";
            ss.str(std::string());
        }
        
        if(rk==true) {
            _basename += "RK";
            ss.str(std::string());
        }
        
        ss<<tau;
        _basename += "tau" + ss.str();
        ss.str(std::string());
        
        return _basename;
    }
    
};

template <class gimpl, class PFermionImpl = PWilsonSmellImplR>
class LangevinRun {

private:
    
    std::ofstream log, plaqfile, liveplaqfile, propfile, normfile;
    
    PRealD plaq;
    std::vector<PRealD> plaqbuffer;
    PComplexD norm;
    
    GridSerialRNG sRNG;
    
    TwistValencePropagator<PFermionImpl> *TVP = NULL;
    
protected:
    
    GridCartesian *grid;
    GridParallelRNG *pRNG;
    
    typename gimpl::GaugeField U;
    
    LangevinParams<PFermionImpl> Params;
    PertLangevin<gimpl> L;
    
    ScidacHmcCheckpointer<gimpl,ActionParameters> CP;
    
public:
    
    LangevinRun(GridCartesian *grid_, GridParallelRNG *pRNG_, LangevinParams<PFermionImpl> Params_):
    grid(grid_),
    pRNG(pRNG_),
    U(grid_),
    Params(Params_),
    L(grid_,pRNG_,Params_.tau,Params_.alpha),
    CP(Params_.CPparams)
    {
        if(Params.save_every!=0) plaqbuffer.resize(Params.save_every);
        
        if(Params.StartingType=="LowerOrderStart"){
            CheckpointerParameters CPparams_low;
            CPparams_low.config_prefix = Params.nameprefix + Params.initialiseName(load_Nplow) + "_lat";
            CPparams_low.rng_prefix = Params.nameprefix + Params.initialiseName(load_Nplow) + "_rng";
            CPparams_low.saveInterval = 1;
            CPparams_low.format = "IEEE64BIG";
            ScidacHmcCheckpointer<TwistedGaugeImpl<GaugeImplTypes_pt<vComplex,Nc,load_Nplow>>,ActionParameters> CPlow(CPparams_low);
            LoadLowerOrder(&CPlow);
            Params.StartTrajectory = -1;
        } else if(Params.StartingType=="ColdStart"){
            pRNG->SeedFixedIntegers(Params.rngseed);
            PertVacuum(U);
        } else if(Params.StartingType=="CheckpointStart"){
            CP.CheckpointRestore(Params.StartTrajectory,U,sRNG,*pRNG);
            if(pog_at_startup){
                std::cout << GridLogMessage << "Projecting checkpoint configuration on group..." << std::endl;
                U = ProjectOnGroup(U);
                std::cout << GridLogMessage << "done." << std::endl;
            }
        }
        
        if(grid->_processor==0){
            openLog();
            openlivePlaq();
            if(Params.save_every!=0) openPlaq();
            if(Params.measureprop) openProp();
            if(Params.enablenorm) openNorm();
        }
        
        if(Params.Nf!=0 && Params.measureprop) TVP = new TwistValencePropagator<PFermionImpl>(grid,Params.FA,Params.prop_phases1,Params.prop_phases2,Params.prop_phases3);
    };
    
    template<class Action>
    void push_back(Action A){
        L.TheActions.push_back(A);
    }
    
    void openLog(){
    
        std::string logname = Params.basename + ".log";
        
        if(Params.StartingType=="CheckpointStart"){
            if(exists(logname))
                log.open(logname, std::ios::app);
            else log.open(logname.c_str());
        }
        else if(!exists(logname)){
            log.open(logname.c_str());
        } else{
            std::cout << GridLogError << "There is already a log file with same parameters: use CheckpointStart to start from a previous checkpoint or move/delete/rename the log for a new ColdStart" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if (!log.is_open()){
            std::cout << GridLogError << "Unable to open log file" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        
        log << std::endl << "/**********************************************/" << std::endl;
        log << now() << std::endl << std::endl;
        log << "ORDERS:       " << Np << std::endl;
        log << "GAUGE GROUP:  SU(" << Nc << ")" << std::endl;
        log << "FERMIONS:     ";
        if(Params.Nf!=0) {
            log << "Nf=" << Params.Nf << " Wilson fermions with mass" << std::endl;
            log << "              " << Params.mass <<std::endl;
            log << "              and boundary phases" << std::endl;
            log << "              " << Params.boundary_phases <<std::endl;
        } else{
            log << "no fermions" << std::endl;
        }
        
        log<<"TWISTED BC:   directions ";
        for (int d=0; d<Nd; d++) {
            if(istwisted(d)) log << d << " ";
        }
        log << std::endl;
        
        log << "LATTICE:      " << GridCmdVectorIntToString(GridDefaultLatt()) << std::endl;
        
        log << std::endl;
        
        log <<"OpenMP threads: "<<GridThread::GetThreads()<<std::endl;
        log <<"MPI tasks:      "<<GridCmdVectorIntToString(GridDefaultMpi())<<std::endl;
        log <<"SIMD layout:    "<<GridCmdVectorIntToString(grid->_simd_layout)<<std::endl;
        
        log << std::endl;
        
        if(Params.rk==true) log << "Integrator:            Runge-Kutta" << std::endl;
        else log << "Integrator:            Euler" << std::endl;
        if(Params.fagf==true) log << "Fourier acceleration:  enabled" << std::endl;
        else log << "Fourier acceleration:  disabled" << std::endl;
        if(Params.measureprop==true) log << "Propagator measure:    enabled" << std::endl;
        else log << "Propagator measure:    disabled" << std::endl;
        log << "StartingType:          " << Params.StartingType << std::endl;
        if(Params.StartingType=="ColdStart") log << "rng seeds:             " << GridCmdVectorIntToString(Params.rngseed) << std::endl;
        if(Params.StartingType=="LowerOrderStart") log << "lower order:               " << load_Nplow << std::endl;
        if(Params.StartingType=="CheckpointStart") log << "starting from:             " << Params.StartTrajectory << std::endl;
        if(Params.gfprecision!=0)
            log << "Saving configurations in Landau gauge with precision " << Params.gfprecision << std::endl;
        log << "POG at startup:        ";
        if(pog_at_startup) log << "YES" << std::endl;
        else log << "NO" << std::endl;
        
        log << std::endl;
        
        log << "tau        = " << Params.tau << std::endl;
        log << "alpha      = " << Params.alpha << std::endl;
        log << "sweeps     = " << Params.sweeps << std::endl;
        log << "save_every = " << Params.save_every << std::endl;
        log << "pog_every  = " << Params.pog_every << std::endl;
        log << std::endl;
        
        log << std::endl;
        
    }
    
    void openPlaq(){
        
        std::string plaqname = Params.basename + ".plaq";
        
        if(Params.StartingType=="CheckpointStart"){
            if(exists(plaqname))
                plaqfile.open(plaqname, std::ios::app);
            else plaqfile.open(plaqname.c_str());
        }
        else if(!exists(plaqname)){
            plaqfile.open(plaqname.c_str());
        } else{
            std::cout << GridLogError << "There is already a plaquette file with same parameters: use CheckpointStart to start from a previous checkpoint or move/delete/rename the plaquette file for a new ColdStart" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if (!plaqfile.is_open()){
            std::cout << GridLogError << "Unable to open plaquette file" << std::endl;
            exit(EXIT_FAILURE);
        }
        plaqfile << std::setprecision(std::numeric_limits<Real>::digits10) << std::scientific;
    }
    
    void openlivePlaq(){
        std::string liveplaqname = Params.basename + ".plaq.live";
        liveplaqfile.open(liveplaqname.c_str());
        if (!liveplaqfile.is_open()){
            std::cout << GridLogError << "Unable to open live plaquette file" << std::endl;
            exit(EXIT_FAILURE);
        }
        liveplaqfile << std::setprecision(std::numeric_limits<Real>::digits10) << std::scientific;
    }
    
    void openNorm(){
        
        std::string normname = Params.basename + ".norm";
        
        if(Params.StartingType=="CheckpointStart"){
            if(exists(normname))
                normfile.open(normname, std::ios::app);
            else normfile.open(normname.c_str());
        }
        else if(!exists(normname)){
            normfile.open(normname.c_str());
        } else{
            std::cout << GridLogError << "There is already a norm file with same parameters: use CheckpointStart to start from a previous checkpoint or move/delete/rename the norm file for a new ColdStart" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if (!normfile.is_open()){
            std::cout << GridLogError << "Unable to open norm file" << std::endl;
            exit(EXIT_FAILURE);
        }
        normfile << std::setprecision(std::numeric_limits<Real>::digits10) << std::scientific;
    }
    
    void openProp(){
        
        std::string propname = Params.basename + ".prop";
        
        if(Params.StartingType=="CheckpointStart"){
            if(exists(propname))
                propfile.open(propname, std::ios::app | std::ios::binary);
            else propfile.open(propname.c_str(), std::ios::binary);
        }
        else if(!exists(propname)){
            propfile.open(propname.c_str(), std::ios::binary);
        } else{
            std::cout << GridLogError << "There is already a propagator file with same parameters: use CheckpointStart to start from a previous checkpoint or move/delete/rename the propagator file for a new ColdStart" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if (!propfile.is_open()){
            std::cout << GridLogError << "Unable to open propagator file" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        propfile.precision(30);
        propfile << std::scientific;
        
    }
    
    void FinaliseRun(){
        log << GridLogMessage << "RUN COMPLETED" << std::endl;
        log << now() << std::endl;
        log.flush();
        log.close();
        if(Params.save_every!=0) plaqfile.close();
        liveplaqfile.close();
        if(Params.measureprop) propfile.close();
        if(Params.enablenorm) normfile.close();
    }
    
    bool exists (const std::string &name) {
        std::ifstream f(name.c_str());
        return f.good();
    }
    
    std::string now( const char* format = "%c" ){
        time_t t = time(0) ;
        char cstr[128] ;
        strftime( cstr, sizeof(cstr), format, localtime(&t) ) ;
        return cstr ;
    }
    
    void Run(){
        
        // print status from time to time
        // e.g. 1000 times per run
        int pp = Params.sweeps/1000;
        if(pp==0) pp = 1;
        
        int plaqcount = 0;
        
        log << std::endl << now() << std::endl;
        log << GridLogMessage << "RUN STARTED" << std::endl;
        log.flush();
        
        for (int i=0; i<Params.sweeps; i++) {
            
            if(i%pp==0) log << GridLogMessage << "Starting sweep number "<< i << " (" << now() << ")" << std::endl;
            log.flush();
            
            if(Params.rk==false) L.EulerStep(U);
            else L.RKStep(U);
            
            L.StochasticGF(U);
            
            if(Params.pog_every!=0 && i%Params.pog_every==(Params.pog_every-1)){
                std::cout << GridLogMessage << "Projecting on group..." << std::endl;
                U = ProjectOnGroup(U);
                std::cout << GridLogMessage << "done." << std::endl;
            }
            
            plaq = WilsonLoops<gimpl>::avgPlaquette(U);
            for (int k=0; k<Np; k++) liveplaqfile << plaq(k) << std::endl;
            if(Params.save_every!=0) plaqbuffer[plaqcount++] = plaq;
            
            if(Params.enablenorm){
                norm = Pnorm2(U);
                for (int k=0; k<Np; k++) normfile << norm(k).real() << std::endl;
            }
            
            if(Params.save_every!=0 && i%Params.save_every==(Params.save_every-1)){
                for (int c=0; c<plaqbuffer.size(); c++) {
                    for (int k=0; k<Np; k++) plaqfile << plaqbuffer[c](k) << std::endl;
                }
                plaqcount = 0;
                
                if(Params.gfprecision!=0){
                    log << GridLogMessage << "Landau gauge fixing ..." <<std::endl;
                    log.flush();
                    if(Params.fagf==true) L.LandauGF_FA(U, Params.gfprecision);
                    else L.LandauGF(U, Params.gfprecision);
                    log << GridLogMessage << "... completed" <<std::endl;
                    log.flush();
                }
                
                if(Params.measureprop) TVP->measure(propfile, U);
                
                CP.TrajectoryComplete(Params.StartTrajectory+i+1,U,sRNG,*pRNG);
                log << GridLogMessage << "Configuration " << Params.StartTrajectory+i+1 << " saved" <<std::endl;
                log.flush();
            }
            
        }
        
        if(grid->_processor==0) FinaliseRun();
    }
    
    template <class impl>
    void LoadLowerOrder(ScidacHmcCheckpointer<impl,ActionParameters> *CPlow){
        Lattice<iVector<iScalar<iPert<iMatrix<vComplex,Nc>,load_Nplow>>,Nd>> Ulow(grid);
        iVector<iScalar<iPert<iMatrix<Complex,Nc>,load_Nplow>>,Nd> tmplow = zero;
        
        CPlow->CheckpointRestore(Params.StartTrajectory,Ulow,sRNG,*pRNG);

        QCDpt::LorentzPertColourMatrix tmp = zero;
        
        for (int x0=0; x0<grid->_fdimensions[0]; x0++) {
            for (int x1=0; x1<grid->_fdimensions[1]; x1++) {
                for (int x2=0; x2<grid->_fdimensions[2]; x2++) {
                    for (int x3=0; x3<grid->_fdimensions[3]; x3++) {
                        peekSite(tmplow,Ulow,std::vector<int>({x0,x1,x2,x3}));
                        for (int mu=0; mu<Nd; mu++) {
                            for (int k=0; k<load_Nplow; k++) {
                                tmp(mu)()(k) = tmplow(mu)()(k);
                            }
                        }
                        pokeSite(tmp,U,std::vector<int>({x0,x1,x2,x3}));
                    }
                }
            }
        }
        U = ProjectOnGroup(U);
    }
    
};

}
}
}
#endif
