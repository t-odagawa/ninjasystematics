#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>

#include <B2Reader.hh>
#include <B2Enum.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>

#include <McsClass.hpp>

namespace logging = boost::log;

int main(int argc, char* argv[]) {

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0] << " <input file> <reweight file> <dial> <sigma> <outputfile>";
    return 1;
  }

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  try {

    std::string imomchname = std::string(argv[1]);
    BOOST_LOG_TRIVIAL(info) << "Input file name : " << imomchname;
    auto ev_vec = Momentum_recon::ReadEventInformationBin(imomchname);
    BOOST_LOG_TRIVIAL(info) << "Event size : " << ev_vec.size();

    B2Reader reader((std::string)argv[2]);

    std::string dial_name = (std::string)argv[3];
    std::string dial_sigma = (std::string)argv[4];
    double sign = 0.;
    if ( dial_sigma == "plus" ) sign = 1.;
    else if ( dial_sigma == "minus" ) sign = -1.;

    for ( auto &ev : ev_vec ) {
      
      int entry = ev.entry_in_daily_file;
      
      reader.ReadSpill(entry);
      auto &spill_summary = reader.GetSpillSummary();
      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      auto interaction_type = event->GetInteractionType();

      double rwfactor = 1.;

      if ( dial_name == "CC_DIS_Norm_Nu" &&
	   interaction_type == B2InteractionMode::MODE_CC_DIS ) {
	rwfactor += (sign * 0.1);
      }
      else if ( dial_name == "CC_MultiPi_Norm_Nu" &&
		interaction_type == B2InteractionMode::MODE_CC_MULTI_PI ) {
	rwfactor += (sign * 0.1);
      }
      else if ( dial_name == "CC_DIS_MultiPi_Norm_Nu" &&
		( interaction_type == B2InteractionMode::MODE_CC_DIS ||
		  interaction_type == B2InteractionMode::MODE_CC_MULTI_PI ) ) {
	rwfactor += (sign * 0.1);
      }

      ev.weight *= rwfactor;      

    }

    std::string ofilename = argv[5];
    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    return 1;
  }

  return 0;

}
