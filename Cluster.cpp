#include "Cluster.h"

#include <algorithm>
#include <iostream>

constexpr int RIGHT_BX{-380};

Cluster::Cluster(std::vector<TriggerPrimitive> const& tps, double x_cut, int wh, int sec, int st)
    : wheel{wh}, sector{sec}, station{st} {
  // in a given wheel-station-sector finds the in-time highest-quality TP
  // and looks for other TP in a |x_cut| interval

  std::vector<TriggerPrimitive> tps_in_chamber;
  // prendo tutte quelle nella stessa wheel-station-sector
  std::copy_if(tps.begin(), tps.end(), std::back_inserter(tps_in_chamber),
               [=](auto& tp) { return tp.wheel == wh && tp.station == st && tp.sector == sec; });

  auto n_tps_in_chamber{tps_in_chamber.size()};

  //riordino il vettore e metto prima le cose in tempo poi i oot ghost
  auto first_oot = std::partition(tps_in_chamber.begin(), tps_in_chamber.end(),
                                  [=](auto& tp) { return tp.BX == RIGHT_BX; });

  // sposto i ghost fuori tempo e poi li cancello da tps_in_chamber
  std::move(first_oot, tps_in_chamber.end(), std::back_inserter(_ootGhosts));
  tps_in_chamber.erase(first_oot, tps_in_chamber.end());

  // se mi è rimasto più di una tp in tempo-> cerco la qualità più alta e la metto in _bestTP

  if (tps_in_chamber.size() > 1){
    auto compareQuality = [](auto& tp1 , auto& tp2) { return tp1.quality < tp2.quality; } ; 
    auto best_tp = ( std::max_element(tps_in_chamber.begin(), tps_in_chamber.end(), compareQuality));
    _bestTP = *best_tp;
  }
  // se ne ho solo una me la tengo buona in _bestTp
  else if (tps_in_chamber.size() == 1) _bestTP = tps_in_chamber[0];


  if (bestTPQuality() >= 0) {
    // se ho la bestTP la cancello dal vettore dei ghost it
    tps_in_chamber.erase(std::remove(tps_in_chamber.begin(), tps_in_chamber.end(), _bestTP));
  }

  _itGhosts = std::move(tps_in_chamber);

}

int Cluster::ootSize() const { return _ootGhosts.size(); };
const std::vector<TriggerPrimitive> & Cluster::ootGhosts() const { return _ootGhosts; };

int Cluster::itSize() const { return _itGhosts.size(); };
const std::vector<TriggerPrimitive> & Cluster::itGhosts() const { return _itGhosts; };


int Cluster::bestTPIndex() const { return _bestTP.index; };
int Cluster::bestTPQuality() const { return _bestTP.quality; };
const TriggerPrimitive & Cluster::bestTP() const { return _bestTP; };

int Cluster::itCountIf(std::function<bool(TriggerPrimitive const&)> f) const {
  return std::count_if(_itGhosts.begin(), _itGhosts.end(), f);
}

int Cluster::ootCountIf(std::function<bool(TriggerPrimitive const&)> f) const {
  return std::count_if(_ootGhosts.begin(), _ootGhosts.end(), f);
}

void Cluster::MatchSegment( int MuWh, int MuStat, int MuSec,  double MuXedge, double MuYedge, double MuX, int MuIndex, int nMu ) {
  if (MuWh == wheel && MuStat == station && MuSec == sector && MuXedge < -5  && MuYedge < -5){
  // if the extrapolated segment is within 5 cm from the HQ it matches
    if (  std::abs( _bestTP.xLoc - MuX ) < 10 ){
        MuMatchedIndex[0] = nMu;
        MuMatchedIndex[1] = MuIndex;
        MuMatched = true;
        //std::cout << " matched" << std::endl;
    }
  }
}