#include "Cluster.h"
#include "Segment.h"

#include <algorithm>
#include <iostream>

constexpr int RIGHT_BX{-380};

Cluster::Cluster(std::vector<TriggerPrimitive> const& tps, std::vector<Segment> const& seg, double xCut, int wh, int sec, int st)
    : wheel{wh}, sector{sec}, station{st} {
  // in a given wheel-station-sector finds the in-time highest-quality TP
  // and looks for other TP in a |xCut| interval 
  // CB TODO: actually xCut not used yet, making 1 cluster per chamber
  // it's cluster even if no TPs but there's a segment

  std::vector<TriggerPrimitive> tps_in_chamber;
  // select primitives from a given chamber: wheel-station-sector
  std::copy_if(tps.begin(), tps.end(), std::back_inserter(tps_in_chamber),
               [=](auto& tp) { return tp.wheel == wh && tp.station == st && tp.sector == sec; });

  auto n_tps_in_chamber{tps_in_chamber.size()};

  //sorting (actually partitioning) the vector:
  // - first part is in-time TPs
  // - second part is out-of-time TPs
  auto first_oot = std::partition(tps_in_chamber.begin(), tps_in_chamber.end(),
                                  [=](auto& tp) { return tp.BX == RIGHT_BX; });

  // moving OOT TPs to dedicated std::vector and removing them from tps_in_chamber
  std::move(first_oot, tps_in_chamber.end(), std::back_inserter(_ootGhosts));
  tps_in_chamber.erase(first_oot, tps_in_chamber.end());

  // If I have more than one in-time TP -> look for the highest quality and assign to _bestTP
  if (tps_in_chamber.size() > 1){
    auto compareQuality = [](auto& tp1 , auto& tp2) { return tp1.quality < tp2.quality; } ; 
    auto best_tp = ( std::max_element(tps_in_chamber.begin(), tps_in_chamber.end(), compareQuality));
    _bestTP = *best_tp;
    foundTP = true;
  }
  // If I have one in-time TP, assign to _bestTP
  else if (tps_in_chamber.size() == 1) {
    _bestTP = tps_in_chamber[0];
    foundTP = true;
  }


  if (bestTPQuality() >= 0) {
    // If I've found _bestTP I remove it from tps_in_chamber
    tps_in_chamber.erase(std::remove(tps_in_chamber.begin(), tps_in_chamber.end(), _bestTP));
  }

  _itGhosts = std::move(tps_in_chamber);

  // segments in same chamber

  std::vector<Segment> segments_in_chamber;  
  std::copy_if(seg.begin(), seg.end(), std::back_inserter(segments_in_chamber),
               [=](auto& segm) { return segm.wheel == wh && segm.station == st && segm.sector == sec; });
  
  if (segments_in_chamber.size() > 0) {
    foundSeg = true;
    if (segments_in_chamber.size() == 1){
      _bestSeg = segments_in_chamber[0];
    }
    
    else{
    auto compareQuality = [](Segment& s1, Segment& s2){return s1.nPhiHits < s2.nPhiHits; };
    auto best_seg = std::max_element(segments_in_chamber.begin(), segments_in_chamber.end(), compareQuality);
    _bestSeg = *best_seg;
    }    
  }

  // match bestSeg e bestTP
  if (foundTP && foundSeg && ( std::abs( _bestSeg.xLoc - _bestTP.xLoc ) < xCut ) ){
    segMatched = true;
    _matchedSeg = _bestSeg;
  }
 

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

void Cluster::MatchMu( int muWh, int muStat, int muSec,  double muXedge, double muYedge, double muX, int muIndex, int iMu ) {
  if (muWh == wheel && muStat == station && muSec == sector && muXedge < -5  && muYedge < -5){
    // if the extrapolated segment is within 10 cm from _bestTP
    if (segMatched && std::abs( _bestSeg.xLoc - muX ) < 10) {
      muMatchedIndex[0] = iMu;
      muMatchedIndex[1] = muIndex;
      muMatched = true;
    }
    if ( !segMatched && std::abs( _bestTP.xLoc - muX ) < 10 ){
      muMatchedIndex[0] = iMu;
      muMatchedIndex[1] = muIndex;
      muMatched = true;
    }
  }
}

int Cluster::matchedSegIndex() const {return _matchedSeg.index; };
int Cluster::matchedSegPhiHits() const {return _matchedSeg.nPhiHits; };
const Segment& Cluster::matchedSeg() const {return _matchedSeg; };

int Cluster::bestSegPhiHits() const {return _bestSeg.nPhiHits; };
const Segment& Cluster::bestSeg() const {return _bestSeg; };

/*void Cluster::MatchSegment(Segment segment, double xCut){
  if (segment.wheel == wheel && segment.station == station && segment.sector == sector){
    if (segMatched == false){
      if ( std::abs(segment.xLoc - _bestTP.xLoc ) < xCut ){
        segMatched = true;
        segIndex = segment.index;
        segHits = segment.nPhiHits;
      }
    }
    else {
      if ( std::abs(segment.xLoc - _bestTP.xLoc ) < xCut && segment.nPhiHits > segHits) {
        segIndex = segment.index;
        segHits = segment.nPhiHits;
      } 
    }
  }
}*/