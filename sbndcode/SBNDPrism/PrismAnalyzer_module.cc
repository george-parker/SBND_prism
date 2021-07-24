////////////////////////////////////////////////////////////////////////
// Class:       PrismAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        PrismAnalyzer_module.cc
//
// Generated at Tue Dec  1 11:40:17 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"

class PrismAnalyzer;


class PrismAnalyzer : public art::EDAnalyzer {
public:
  explicit PrismAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PrismAnalyzer(PrismAnalyzer const&) = delete;
  PrismAnalyzer(PrismAnalyzer&&) = delete;
  PrismAnalyzer& operator=(PrismAnalyzer const&) = delete;
  PrismAnalyzer& operator=(PrismAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(art::SubRun const& sr) override;
  virtual void endSubRun(art::SubRun const& sr) override;
  void endJob();

private:
  /// Returns the off-axis angle in radians
  float GetOffAxisAngle(float x, float y, float z);

  /// Returns the reco energy using QE formula
  float GetEnergyQE(float muon_energy, float muon_px, float muon_py, float muon_pz);

  std::string _mctruth_producer = "generator";
  float _beam_origin_x = 73.78;
  float _beam_origin_y = 0.0;
  float _beam_origin_z = 11000.0;
  // std::vector<float> _beam_origin; // = {-0.457 * 100, 0 * 100, 110 * 100}; // cm

  const TDatabasePDG *_pdg_db = TDatabasePDG::Instance();

  TTree* _tree;

  int _run, _subrun, _event;
  
  //trial
  float _nu_e; ///< Neutrino energy
  long int _nu_pdg; ///< Neutrino PDG code
  int _nu_ccnc; ///< 0: CC, 1: NC
  int _nu_mode; ///< Neutrino interaction mode
  int _nu_int_type; ///< Neutrino interaction type
  float _nu_vtx_x; ///< Neutrino vertex X
  float _nu_vtx_y; ///< Neutrino vertex Y
  float _nu_vtx_z; ///< Neutrino vertex Z
  float _nu_vtx_t; ///< Neutrino vertex T
  float _nu_px; ///< Neutrino momentum along X
  float _nu_py; ///< Neutrino momentum along Y
  float _nu_pz; ///< Neutrino momentum along Z
  float _nu_oaa; ///< Off-Axis Angle (angle w.r.t. beam)
  float _nu_roaa; ///< Real Off-Axis Angle (angle between neutrino parent p and neutrino p)
  float _nu_l; ///< L [cm]: Distance the neutrino travelled from production to interaction point
  int _nu_decay; ///< Neutrino parent decay type (see dkproc_t in dk2nu)
 
  float _nu_lepton_e; ///< Final state lepton energy
  float _nu_lepton_px; ///< Final state lepton px
  float _nu_lepton_py; ///< Final state lepton py
  float _nu_lepton_pz; ///< Final state lepton pz

  float _nu_electron_e; ///< Final state electron energy
  float _nu_electron_px; ///< Final state electron px
  float _nu_electron_py; ///< Final state electron py
  float _nu_electron_pz; ///< Final state electron pz

  float _nu_positron_e; ///< Final state positron energy
  float _nu_positron_px; ///< Final state positron px
  float _nu_positron_py; ///< Final state positron py
  float _nu_positron_pz; ///< Final state positron pz

  float _nu_photon1_e; ///< Final state photon energy
  float _nu_photon1_px; ///< Final state photon px
  float _nu_photon1_py; ///< Final state photon py
  float _nu_photon1_pz; ///< Final state photon pz

  float _nu_photon2_e; ///< Final state photon energy
  float _nu_photon2_px; ///< Final state photon px
  float _nu_photon2_py; ///< Final state photon py
  float _nu_photon2_pz; ///< Final state photon pz

  float _opening_angle; ///<Opening angle of electron-positron

  float _nu_e_reco; ///< Neutrino reconstructe energy using QE formula
  std::vector<long int> _pars_pdg; ///< All other particles produced - pdg code
  std::vector<float> _pars_e; ///< All other particles produced - energy

  std::vector<long int> pdg_num; ///< particles produced - pdg code
  std::vector<long int> pdg_count; ///<  particles produced - pdg count
  std::vector<long int> pdg_row; ///< particles produced 
  std::vector<long int> pdg_single; ///<  particles produced
  std::vector<long int> pdg_double; ///< particles produced
  std::vector<long int> pdg_other; ///<  particles produced
  int counter;

  float _nu_prod_vtx_x; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_y; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_z; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_x_beam; ///< Neutrino production vertex in beamline coordinates
  float _nu_prod_vtx_y_beam; ///< Neutrino production vertex in beamline coordinates
  float _nu_prod_vtx_z_beam; ///< Neutrino production vertex in beamline coordinates

  int _p_type; ///< Neutrino parent PDG code
  float _p_dpx; ///< Neutrino parent px at neutrino production vertex
  float _p_dpy; ///< Neutrino parent py at neutrino production vertex
  float _p_dpz; ///< Neutrino parent pz at neutrino production vertex

  int _nu_pip_mult; ///< Pi0 multiplicity
  int _nu_pi0_mult; ///< Pi plus multiplicity
  int _nu_p_mult; ///< Proton multiplicity
  
  int _lept_pairs;
  int _phot_pairs;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun

};


PrismAnalyzer::PrismAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  _phot_pairs = 0;
  _lept_pairs = 0;
  counter = 0;

  _beam_origin_x = p.get<float>("BeamCenterX");
  _beam_origin_y = p.get<float>("BeamCenterY");
  _beam_origin_z = p.get<float>("BeamCenterZ");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree","");

  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");

  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
  _tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
  _tree->Branch("nu_int_type", &_nu_int_type, "nu_int_type/I");
  _tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
  _tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
  _tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
  _tree->Branch("nu_vtx_t", &_nu_vtx_t, "nu_vtx_t/F");
  _tree->Branch("nu_px", &_nu_px, "nu_px/F");
  _tree->Branch("nu_py", &_nu_py, "nu_py/F");
  _tree->Branch("nu_pz", &_nu_pz, "nu_pz/F");
  _tree->Branch("nu_oaa", &_nu_oaa, "nu_oaa/F");
  _tree->Branch("nu_roaa", &_nu_roaa, "nu_roaa/F");
  _tree->Branch("nu_l", &_nu_l, "nu_l/F");
  _tree->Branch("nu_decay", &_nu_decay, "nu_decay/I");
  _tree->Branch("nu_lepton_e", &_nu_lepton_e, "nu_lepton_e/F");
  _tree->Branch("nu_lepton_px", &_nu_lepton_px, "nu_lepton_px/F");
  _tree->Branch("nu_lepton_py", &_nu_lepton_py, "nu_lepton_py/F");
  _tree->Branch("nu_lepton_pz", &_nu_lepton_pz, "nu_lepton_pz/F");

  _tree->Branch("opening_angle", &_opening_angle, "opening_angle/F");

  _tree->Branch("nu_electron_e", &_nu_electron_e, "nu_electron_e/F");
  _tree->Branch("nu_electron_px", &_nu_electron_px, "nu_electron_px/F");
  _tree->Branch("nu_electron_py", &_nu_electron_py, "nu_electron_py/F");
  _tree->Branch("nu_electron_pz", &_nu_electron_pz, "nu_electron_pz/F");

  _tree->Branch("nu_positron_e", &_nu_positron_e, "nu_positron_e/F");
  _tree->Branch("nu_positron_px", &_nu_positron_px, "nu_positron_px/F");
  _tree->Branch("nu_positron_py", &_nu_positron_py, "nu_positron_py/F");
  _tree->Branch("nu_positron_pz", &_nu_positron_pz, "nu_positron_pz/F");

  _tree->Branch("nu_photon1_e", &_nu_photon1_e, "nu_photon1_e/F");
  _tree->Branch("nu_photon1_px", &_nu_photon1_px, "nu_photon1_px/F");
  _tree->Branch("nu_photon1_py", &_nu_photon1_py, "nu_photon1_py/F");
  _tree->Branch("nu_photon1_pz", &_nu_photon1_pz, "nu_photon1_pz/F");

  _tree->Branch("nu_photon2_e", &_nu_photon2_e, "nu_photon2_e/F");
  _tree->Branch("nu_photon2_px", &_nu_photon2_px, "nu_photon2_px/F");
  _tree->Branch("nu_photon2_py", &_nu_photon2_py, "nu_photon2_py/F");
  _tree->Branch("nu_photon2_pz", &_nu_photon2_pz, "nu_photon2_pz/F");

  _tree->Branch("nu_e_reco", &_nu_e_reco, "nu_e_reco/F");

  _tree->Branch("nu_prod_vtx_x", &_nu_prod_vtx_x, "nu_prod_vtx_x/F");
  _tree->Branch("nu_prod_vtx_y", &_nu_prod_vtx_y, "nu_prod_vtx_y/F");
  _tree->Branch("nu_prod_vtx_z", &_nu_prod_vtx_z, "nu_prod_vtx_z/F");
  _tree->Branch("nu_prod_vtx_x_beam", &_nu_prod_vtx_x_beam, "nu_prod_vtx_x_beam/F");
  _tree->Branch("nu_prod_vtx_y_beam", &_nu_prod_vtx_y_beam, "nu_prod_vtx_y_beam/F");
  _tree->Branch("nu_prod_vtx_z_beam", &_nu_prod_vtx_z_beam, "nu_prod_vtx_z_beam/F");

  _tree->Branch("pars_pdg", "std::vector<long int>", &_pars_pdg);
  _tree->Branch("pars_e", "std::vector<float>", &_pars_e);

  _tree->Branch("p_type", &_p_type, "p_type/I");
  _tree->Branch("p_dpx", &_p_dpx, "p_dpx/F");
  _tree->Branch("p_dpy", &_p_dpy, "p_dpy/F");
  _tree->Branch("p_dpz", &_p_dpz, "p_dpz/F");

  _tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
  _tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
  _tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void PrismAnalyzer::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  //
  // Get the associated MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_producer, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product " << _mctruth_producer << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

  //
  // Get the associated MCFlux
  //
  art::FindManyP<simb::MCFlux> mct_to_mcf (mct_h, e, _mctruth_producer);

  // _pars_pdg.clear();
  //_pars_e.clear();
  _nu_pi0_mult = _nu_pip_mult = _nu_p_mult = 0;
  //counter = 0;

  //
  // Loop over the neutrino interactions in this event
  //
  for (size_t i = 0; i < mct_v.size(); i++) {
    if (mct_v.at(i)->Origin() != simb::kBeamNeutrino) {
      std::cout << "[PrismAnalyzer] MCTruth from generator does not have neutrino origin?!" << std::endl;
    }

    _nu_e = mct_v[i]->GetNeutrino().Nu().E();
    _nu_pdg = mct_v[i]->GetNeutrino().Nu().PdgCode();
    _nu_ccnc = mct_v[i]->GetNeutrino().CCNC();
    _nu_mode = mct_v[i]->GetNeutrino().Mode();
    _nu_int_type = mct_v[i]->GetNeutrino().InteractionType();
    _nu_vtx_x = mct_v[i]->GetNeutrino().Nu().Vx();
    _nu_vtx_y = mct_v[i]->GetNeutrino().Nu().Vy();
    _nu_vtx_z = mct_v[i]->GetNeutrino().Nu().Vz();
    _nu_vtx_t = mct_v[i]->GetNeutrino().Nu().T();
    _nu_px = mct_v[i]->GetNeutrino().Nu().Px();
    _nu_py = mct_v[i]->GetNeutrino().Nu().Py();
    _nu_pz = mct_v[i]->GetNeutrino().Nu().Pz();

    _nu_oaa = GetOffAxisAngle(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z);

    // Get the associated MCFlux
    std::vector<art::Ptr<simb::MCFlux>> mcf_v = mct_to_mcf.at(i);
    assert(mcf_v.size() == 1);
    auto mcf = mcf_v[0];
    assert(mcf->fntype == _nu_pdg);

    // Save the neutrino production vertex, this is in beam coordinates
    _nu_prod_vtx_x_beam = mcf->fvx;
    _nu_prod_vtx_y_beam = mcf->fvy;
    _nu_prod_vtx_z_beam = mcf->fvz;

    // Also convert it to detector coordinates
    _nu_prod_vtx_x = mcf->fvx - _beam_origin_x;
    _nu_prod_vtx_y = mcf->fvy - _beam_origin_y;
    _nu_prod_vtx_z = mcf->fvz - _beam_origin_z;

    auto diff_x = _nu_vtx_x - _nu_prod_vtx_x;
    auto diff_y = _nu_vtx_y - _nu_prod_vtx_y;
    auto diff_z = _nu_vtx_z - _nu_prod_vtx_z;
    _nu_l = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);

    _nu_decay = mcf->fndecay;
    _p_type = mcf->fptype;
    _p_dpx = mcf->fpdpx;
    _p_dpy = mcf->fpdpy;
    _p_dpz = mcf->fpdpz;
    _nu_roaa = TVector3(_nu_px, _nu_py, _nu_pz).Angle(TVector3(_p_dpx, _p_dpy, _p_dpz));

    _nu_lepton_e = -9999;
    _nu_lepton_px = -9999;
    _nu_lepton_py = -9999;
    _nu_lepton_pz = -9999;
    _nu_e_reco = -9999;

    _nu_electron_e = -9999;
    _nu_electron_px = -9999;
    _nu_electron_py = -9999;
    _nu_electron_pz = -9999;

    _nu_positron_e = -9999;
    _nu_positron_px = -9999;
    _nu_positron_py = -9999;
    _nu_positron_pz = -9999;

    _nu_photon1_e = -9999;
    _nu_photon1_px = -9999;
    _nu_photon1_py = -9999;
    _nu_photon1_pz = -9999;

    _nu_photon2_e = -9999;
    _nu_photon2_px = -9999;
    _nu_photon2_py = -9999;
    _nu_photon2_pz = -9999;

    float electron_e = -9999;
    float electron_px = -9999;
    float electron_py = -9999;
    float electron_pz = -9999;

    float positron_e = -9999;
    float positron_px = -9999;
    float positron_py = -9999;
    float positron_pz = -9999;

    float photon1_e = -9999;
    float photon1_px = -9999;
    float photon1_py = -9999;
    float photon1_pz = -9999;

    float photon2_e = -9999;
    float photon2_px = -9999;
    float photon2_py = -9999;
    float photon2_pz = -9999;
    
    int _positron = 0;
    int _electron = 0;

    for (int p = 0; p < mct_v[i]->NParticles(); p++) {
      auto const & mcp = mct_v[i]->GetParticle(p);
      
      if (mcp.StatusCode() != 1) continue;

      if (std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode()) ==pdg_num.end()) {
	// someName not in name, add it
	pdg_num.push_back(mcp.PdgCode());
	pdg_count.push_back(0);
	pdg_row.push_back(0);
	pdg_single.push_back(0);
	pdg_double.push_back(0);
	pdg_other.push_back(0);
      }

      //std::vector<long int>::iterator it = std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode());       
      //std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode());
      //std::cout << "The index is at:" << std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode()) - pdg_num.begin() << std::endl;
      pdg_count[std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode()) - pdg_num.begin()]++;
      pdg_row[std::find(pdg_num.begin(), pdg_num.end(), mcp.PdgCode()) - pdg_num.begin()]++;
      //pdg_num.push_back(mcp.PdgCode());

      _pars_pdg.push_back(mcp.PdgCode());
      _pars_e.push_back(mcp.E());

      if (mcp.PdgCode() == 111) {
        _nu_pi0_mult++;
	
      } else if (std::abs(mcp.PdgCode()) == 211) {
        _nu_pip_mult++;
      }
      else if (mcp.PdgCode() == 22) {
	//_phot_pairs++;
	//std::cout << "Photon:" << _phot_pairs << std::endl;
	
        if (photon1_e == -9999) {
	  photon1_e = mcp.E();
          photon1_px = mcp.Px();
	  photon1_py = mcp.Py();
	  photon1_pz = mcp.Pz(); }
	else {
	  photon2_e = mcp.E();
          photon2_px = mcp.Px();
          photon2_py = mcp.Py();
          photon2_pz = mcp.Pz(); 

	  _nu_photon1_e = photon1_e;
	  _nu_photon1_px = photon1_px;
	  _nu_photon1_py = photon1_py;
	  _nu_photon1_pz = photon1_pz;

	  _nu_photon2_e = photon2_e;
	  _nu_photon2_px = photon2_px;
	  _nu_photon2_py = photon2_py;
	  _nu_photon2_pz = photon2_pz;
	  
	  _phot_pairs++;
	  //std::cout << "Photon pairs is:" << _phot_pairs << std::endl;
	}
      } 


      else if (std::abs(mcp.PdgCode()) == 2112) {
	_nu_p_mult++;
      } else if (std::abs(mcp.PdgCode()) == 13 || std::abs(mcp.PdgCode()) == 11) {
        _nu_lepton_e = mcp.E();
        _nu_lepton_px = mcp.Px();
        _nu_lepton_py = mcp.Py();
        _nu_lepton_pz = mcp.Pz();
        _nu_e_reco = GetEnergyQE(mcp.E(), mcp.Px(), mcp.Py(), mcp.Pz());

	if (mcp.PdgCode()  == 11) {
	  _electron = 1;
  
	  electron_e = mcp.E();
	  electron_px = mcp.Px();
	  electron_py = mcp.Py();
	  electron_pz = mcp.Pz();
	  
	  if (_positron == 1) {

	    _nu_electron_e = electron_e;
	    _nu_electron_px = electron_px;
	    _nu_electron_py = electron_py;
	    _nu_electron_pz = electron_pz;

	    _nu_positron_e = positron_e;
	    _nu_positron_px = positron_px;
	    _nu_positron_py = positron_py;
	    _nu_positron_pz = positron_pz;
	    
	    _lept_pairs++;
	    //std::cout << "Number of positrons:" << _lept_pairs << std::endl;
	    _opening_angle = std::acos((_nu_electron_px * _nu_positron_px + _nu_electron_py * _nu_positron_py + _nu_electron_pz * _nu_positron_pz)/(std::sqrt(_nu_electron_px * _nu_electron_px + _nu_electron_py * _nu_electron_py + _nu_electron_pz * _nu_electron_pz)*std::sqrt(_nu_positron_px * _nu_positron_px + _nu_positron_py * _nu_positron_py + _nu_positron_pz * _nu_positron_pz)));
	  }
        }

        if (mcp.PdgCode()  == -11) {
	  _positron = 1;

	  //std::cout << "Number of positrons:" << _lept_pairs << std::endl;

          positron_e = mcp.E();
          positron_px = mcp.Px();
          positron_py = mcp.Py();
          positron_pz = mcp.Pz();
	  
	  if (_electron == 1) {

            _nu_electron_e = electron_e;
            _nu_electron_px = electron_px;
            _nu_electron_py = electron_py;
            _nu_electron_pz = electron_pz;

            _nu_positron_e = positron_e;
            _nu_positron_px = positron_px;
            _nu_positron_py = positron_py;
            _nu_positron_pz = positron_pz;

	    _lept_pairs++;
	    //std::cout << "Number of e+e- pairs is:" << _lept_pairs << std::endl;
	    _opening_angle = std::acos((_nu_electron_px * _nu_positron_px + _nu_electron_py * _nu_positron_py + _nu_electron_pz * _nu_positron_pz)/(std::sqrt(_nu_electron_px * _nu_electron_px + _nu_electron_py * _nu_electron_py + _nu_electron_pz * _nu_electron_pz)*std::sqrt(_nu_positron_px * _nu_positron_px + _nu_positron_py * _nu_positron_py + _nu_positron_pz * _nu_positron_pz)));

          }
	  //_positron.push_back(p);
        }
      }

      //for (int x = 0; x < (int)_electron.size(); x++) {
      //auto const & mcp_e = mct_v[i]->GetParticle(_electron[x]);
      //for (int y = 0; y < (int)_positron.size(); y++) {
      // auto const & mcp_p = mct_v[i]->GetParticle(_positron[y]);
      // if (mcp_e->Mother() == mcp_p->Mother()) {
      // _nu_electron_e = mcp_e.E();
      //_nu_positron_e = mcp_p.E();}}}
    }
 

    // _n_pi0 = 0;
    // for (int p = 0; p < mct_v[i]->NParticles(); p++) {
    //   auto const & mcp = mct_v[i]->GetParticle(p);
    //   if (mcp.StatusCode() != 1) continue;
    //   if (mcp.PdgCode() != 111) continue;
    //   _n_pi0++;
    // }
    // for (int p = 0; p < mct_v[i]->NParticles(); p++) {
    //   auto const & mcp = mct_v[i]->GetParticle(p);
    //   if (mcp.StatusCode() != 1) continue;
    //   if (mcp.PdgCode() != 2112) continue;
    //   std::cout << "Got neutron, mother is " << mcp.Mother() << ", energy is " << mcp.E() << std::endl;
    // }

    _tree->Fill();

    for(std::vector<long int>::size_type y=0; y<pdg_row.size(); ++y) {
      if (pdg_row[y] == 1) {
        pdg_single[y]++;}
      else if (pdg_row[y] == 2) {
	pdg_double[y]++;}
      else if (pdg_row[y] > 2) {
        pdg_other[y]++;}
    }
    counter += std::accumulate(pdg_row.begin(), pdg_row.end(), 0);
    std::fill(pdg_row.begin(), pdg_row.end(), 0);
  }
  

  TDatabasePDG *db   = TDatabasePDG::Instance();

   std::ofstream myfile;
   myfile.open("pdg.txt");
   myfile << std::setw(10) << "Particle:" << std::setw(15) << "PDG:" << std::setw(10) << "Total:" << std::setw(10) << "Single:" << std::setw(10) << "Double:" << std::setw(10) << "More:" << " \n" ;
   for(std::vector<long int>::size_type z=0; z<pdg_num.size(); ++z) {
     if (pdg_num[z] < 10000) {
       TParticlePDG * p = db->GetParticle(pdg_num[z]);
       myfile << std::setw(10) << p->GetName() << std::setw(15) << pdg_num[z] << std::setw(10) << pdg_count[z] << std::setw(10) << pdg_single[z] << std::setw(10) << pdg_double[z] << std::setw(10) << pdg_other[z] << "\n" ;}
     else {
       myfile << std::setw(10) << "Nuclei"  << std::setw(15) << pdg_num[z] << std::setw(10) << pdg_count[z] << std::setw(10) << pdg_single[z] << std::setw(10) << pdg_double[z] << std::setw(10) << pdg_other[z] << "\n" ;}
   }  
   myfile << "\n" << "Total number of particles         :" << counter << "\n";
   myfile << "Total number of particles (check) :" << std::accumulate(pdg_count.begin(), pdg_count.end(), 0) << "\n\n";
   myfile << "Total number of events w/ at least:"  <<  "\n";
   myfile << "2 photons: " << _phot_pairs << "\n";
   myfile << "e+ and e-: " << _lept_pairs << "\n";
   myfile.close();

  // art::Handle<std::vector<simb::MCParticle>> mcp_h;
  // e.getByLabel("largeant", mcp_h);
  // if(!mcp_h.isValid()){
  //   std::cout << "MCParticle product not found..." << std::endl;
  //   throw std::exception();
  // }
  // std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  // art::fill_ptr_vector(mcp_v, mcp_h);
  // for (size_t p = 0; p < mcp_v.size(); p++) {
  //   auto mcp = mcp_v[p];
  //   if (mcp->StatusCode() != 1) continue;
  //   if (mcp->Mother() != 1) continue;
  //   std::cout << "MCParticle " << p
  //             << ", trackid " << mcp->TrackId()
  //             << ", pdg " << mcp->PdgCode()
  //             << ", E " << mcp->E()
  //             << ", mother " << mcp->Mother()
  //             << ", process " << mcp->Process() << std::endl;
  // }
}

float PrismAnalyzer::GetOffAxisAngle(float x, float y, float z) {
  TVector3 nu_vtx(x, y, z);
  TVector3 beam_center(_beam_origin_x, _beam_origin_y, _beam_origin_z);
  TVector3 nu_vtx_beam = nu_vtx + beam_center;
  TVector3 beam(0, 0, 1);

  return beam.Angle(nu_vtx_beam);
}

float PrismAnalyzer::GetEnergyQE(float muon_energy, float px, float py, float pz) {
  float p_mass = _pdg_db->GetParticle(2212)->Mass(); // proton mass
  float n_mass = _pdg_db->GetParticle(2112)->Mass(); // neutron mass
  float m_mass = _pdg_db->GetParticle(13)->Mass(); // muon mass
  float Eb = 0.030; // binding energy

  float p = std::sqrt(px * px + py * py + pz * pz);
  // float pxy = TMath::Sqrt( px*px + py*py )
  float theta = TMath::ACos(pz / p);

  float num = p_mass * p_mass;
  num -= m_mass * m_mass;
  num -= (n_mass - Eb) * (n_mass - Eb);
  num += 2 * (n_mass - Eb) * muon_energy;

  float den = 2. * (n_mass - Eb - muon_energy + p * std::cos(theta));

  return num / den;
}


void PrismAnalyzer::beginSubRun(art::SubRun const& sr) {

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(_mctruth_producer, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }
  std::cout << "POT for this subrun: " << _sr_pot << std::endl;

  _sr_tree->Fill();

}

void PrismAnalyzer::endJob() {

}

void PrismAnalyzer::endSubRun(art::SubRun const& sr) {

}

DEFINE_ART_MODULE(PrismAnalyzer)

















