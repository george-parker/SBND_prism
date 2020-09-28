////////////////////////////////////////////////////////////////////////
//// Class:       SBNDFlashFinder
//// Module Type: producer
//// File:        SBNDFlashFinder_module.cc
////
//// Adapted form ICARUSFlashFinder by Kazuhiro Terao
//////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>
#include <string>
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderManager.h"
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderFMWKInterface.h"
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/PECalib.h"

namespace opdet{

  class SBNDFlashFinder;

  class SBNDFlashFinder : public art::EDProducer {
  public:
    explicit SBNDFlashFinder(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SBNDFlashFinder(SBNDFlashFinder const &) = delete;
    SBNDFlashFinder(SBNDFlashFinder &&) = delete;
    SBNDFlashFinder & operator = (SBNDFlashFinder const &) = delete;
    SBNDFlashFinder & operator = (SBNDFlashFinder &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

  private:

    ::lightana::FlashFinderManager _mgr;
    ::lightana::PECalib _pecalib;
    std::vector<std::string> _hit_producers;

    void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);

  };

  SBNDFlashFinder::SBNDFlashFinder(lightana::Config_t const & p)
  : EDProducer{p}
  // Initialize member data here.
  {
    _hit_producers = p.get<std::vector<std::string>>("OpHitProducers");

    auto const flash_algo  = p.get<std::string>("FlashFinderAlgo");
    auto const flash_pset = p.get<lightana::Config_t>("AlgoConfig");
    auto algo_ptr = ::lightana::FlashAlgoFactory::get().create(flash_algo,flash_algo);
    algo_ptr->Configure(flash_pset);
    _mgr.SetFlashAlgo(algo_ptr);
    _pecalib.Configure(p.get<lightana::Config_t>("PECalib"));

    produces< std::vector<recob::OpFlash>   >();
    produces< art::Assns <recob::OpHit, recob::OpFlash> >();
  }

  void SBNDFlashFinder::produce(art::Event & e)
  {

    // produce OpFlash data-product to be filled within module
    std::unique_ptr< std::vector<recob::OpFlash> > opflashes(new std::vector<recob::OpFlash>);
    std::unique_ptr< art::Assns <recob::OpHit, recob::OpFlash> > flash2hit_assn_v
      (new art::Assns<recob::OpHit, recob::OpFlash>);

    std::vector<art::Ptr<recob::OpHit>> ophit_v;

    ::lightana::LiteOpHitArray_t ophits;
    double trigger_time=1.1e20;

    for (auto producer : _hit_producers) {

      // load OpHits previously created
      art::Handle<std::vector<recob::OpHit> > ophit_h;
      e.getByLabel(producer, ophit_h);

      // make sure hits look good
      if(!ophit_h.isValid()) {
        std::cerr << "\033[93m[ERROR]\033[00m ... could not locate OpHit!" << std::endl;
        throw std::exception();
      }

      std::vector<art::Ptr<recob::OpHit>> temp_v;
      art::fill_ptr_vector(temp_v, ophit_h);
      ophit_v.insert(ophit_v.end(), temp_v.begin(), temp_v.end());
    }

    for(auto const oph : ophit_v) {
      ::lightana::LiteOpHit_t loph;
      if(trigger_time > 1.e20) trigger_time = oph->PeakTimeAbs() - oph->PeakTime();
      loph.peak_time = oph->PeakTime();

      size_t opdet = ::lightana::OpDetFromOpChannel(oph->OpChannel());
      loph.pe = _pecalib.Calibrate(opdet,oph->Area());
      loph.channel = oph->OpChannel();
      ophits.emplace_back(std::move(loph));
    }

    auto const flash_v = _mgr.RecoFlash(ophits);

    for(const auto& lflash :  flash_v) {

      double Ycenter, Zcenter, Ywidth, Zwidth;
      GetFlashLocation(lflash.channel_pe, Ycenter, Zcenter, Ywidth, Zwidth);
      recob::OpFlash flash(lflash.time, lflash.time_err, trigger_time + lflash.time,
                         (trigger_time + lflash.time) / 1600., lflash.channel_pe,
                         0, 0, 1, // this are just default values
                         Ycenter, Ywidth, Zcenter, Zwidth);
      opflashes->emplace_back(std::move(flash));


      for(auto const& hitidx : lflash.asshit_idx) {
        const art::Ptr<recob::OpHit> hit_ptr(ophit_v.at(hitidx));
        util::CreateAssn(*this, e, *opflashes, hit_ptr, *flash2hit_assn_v);
      }
    }

    e.put(std::move(opflashes));
    e.put(std::move(flash2hit_assn_v));
  }

  void SBNDFlashFinder::GetFlashLocation(std::vector<double> pePerOpChannel,
                                         double& Ycenter,
                                         double& Zcenter,
                                         double& Ywidth,
                                         double& Zwidth)
  {

    // Reset variables
    Ycenter = Zcenter = 0.;
    Ywidth  = Zwidth  = -999.;
    double totalPE = 0.;
    double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
    for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {
      if (opch > 31 && opch < 100){  //TO ADAPT FOR SBND USING THE MAP
        //  std::cout << "Ignoring channel " << opch << " as it's not a real channel" << std::endl;
        continue;
      }
      // Get physical detector location for this opChannel
      double PMTxyz[3];
      ::lightana::OpDetCenterFromOpChannel(opch, PMTxyz);
      // Add up the position, weighting with PEs
      sumy    += pePerOpChannel[opch]*PMTxyz[1];
      sumy2   += pePerOpChannel[opch]*PMTxyz[1]*PMTxyz[1];
      sumz    += pePerOpChannel[opch]*PMTxyz[2];
      sumz2   += pePerOpChannel[opch]*PMTxyz[2]*PMTxyz[2];
      totalPE += pePerOpChannel[opch];
    }

    Ycenter = sumy/totalPE;
    Zcenter = sumz/totalPE;

    // This is just sqrt(<x^2> - <x>^2)
    if ( (sumy2*totalPE - sumy*sumy) > 0. )
      Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

    if ( (sumz2*totalPE - sumz*sumz) > 0. )
      Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
  }

  DEFINE_ART_MODULE(SBNDFlashFinder)

}//end namespace
