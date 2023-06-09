#include "Cluster.h"

#include <algorithm>
#include <iostream>


constexpr int RIGHT_BX{-380};


Cluster::Cluster(std::vector<TriggerPrimitive> const& tps, double x_cut, int wh, int sec, int st) 
    : wheel{wh}, sector{sec}, station{st} {
  // in a given wheel-station-sector finds the in-time highest-quality TP
  // and looks for other TP in a |x_cut| interval

  std::vector<TriggerPrimitive> tps_in_chamber;

    // Put all the TP in the right wheel-station-sector in the vector 
  std::copy_if(tps.begin(), tps.end(), std::back_inserter(tps_in_chamber),
               [=](auto& tp) { return tp.wheel == wh && tp.station == st && tp.sector == sec; });

  auto n_tps_in_chamber{tps_in_chamber.size()};
    
    // put at the beginning of the vector the tps with the right BX, mi returna l'indice della prima col BX sbagliato
  auto first_oot = std::partition(tps_in_chamber.begin(), tps_in_chamber.end(), [=](auto& tp) { return tp.BX == RIGHT_BX; });

    //sposto le ultime nel vettore dei oot ghost
  std::move(first_oot, tps_in_chamber.end(), std::back_inserter(_ootGhosts));

    // cancello gli oot ghost dal vettore originale in cui ho it ghost e best quality
  tps_in_chamber.erase(first_oot, tps_in_chamber.end());

    //bool compareBX = [](auto tp1 , auto tp2) { return tp1.quality < tp2.quality; } ; 
   // std::max_element(tps_in_chamber.begin(), tps_in_chamber.end(), compareBX);
    //std::cout << _bestTP.quality << endl;
 // find "best" TP in wh st sec CB std::max
for (const auto& tp : tps_in_chamber) {
    if (tp.BX == RIGHT_BX && tp.quality > bestTPQuality()) {
      _bestTP = tp;
    }
  }

  if (bestTPQuality() >= 0) {
    tps_in_chamber.erase(std::remove(tps_in_chamber.begin(), tps_in_chamber.end(), _bestTP));
  }

  _itGhosts = std::move(tps_in_chamber);

};

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


