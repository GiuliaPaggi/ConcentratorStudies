#include "Cluster.h"
#include "Segment.h"

#include <algorithm>
#include <iostream>

constexpr int RIGHT_BX{-380};

Cluster::Cluster(std::vector<TriggerPrimitive> & tps, std::vector<Segment> & seg, std::vector<Digi> & digis, double xCut, double digiCut, int wh, int sec, int st)
    : wheel{wh}, sector{sec}, station{st} {
  // in a given wheel-station-sector finds the in-time highest-quality TP
  // and looks for other TP in a |xCut| interval 
  // CB TODO: actually xCut not used yet, making 1 cluster per chamber
  // it's cluster even if no TPs but there's a segment, or more than 10 digis in a SL


  // ########## TPs cluster #############
  std::vector<TriggerPrimitive> tps_in_chamber;
  // select primitives from a given chamber: wheel-station-sector
  std::copy_if(tps.begin(), tps.end(), std::back_inserter(tps_in_chamber),
               [=](auto& tp) { return tp.wheel == wh && tp.station == st && tp.sector == sec && tp.inCluster == false; });

  auto n_tps_in_chamber{tps_in_chamber.size()};
  
  if (n_tps_in_chamber > 0 ) {
    // cluster closeby TPS from tps in chamber 
    std::vector<TriggerPrimitive> cluster;
    
    std::copy_if( tps_in_chamber.begin(), tps_in_chamber.end(), std::back_inserter(cluster), 
                  [=] (auto &tp) {return tp.wheel == wh && tp.station == st && tp.sector == sec && std::abs(tp.xLoc-tps_in_chamber[0].xLoc) < xCut;} );
    for (auto &tp : tps){
      for (auto &Cl_tp :cluster){
        if (Cl_tp.index == tp.index) tp.inCluster = true;
      }
    }

    //sorting (actually partitioning) the vector:
    // - first part is in-time TPs
    // - second part is out-of-time TPs
    auto first_oot = std::partition(cluster.begin(), cluster.end(),
                                    [=](auto& tp) { return tp.BX == RIGHT_BX; });

    // moving OOT TPs to dedicated std::vector and removing them from tps_in_chamber
    std::move(first_oot, cluster.end(), std::back_inserter(_ootGhosts));
    cluster.erase(first_oot, cluster.end());

    // If I have more than one in-time TP -> look for the highest quality and assign to _bestTP
    if (cluster.size() > 1){
      auto compareQuality = [](auto& tp1 , auto& tp2) { return tp1.quality < tp2.quality; } ; 
      auto best_tp = ( std::max_element(cluster.begin(), cluster.end(), compareQuality));
      _bestTP = *best_tp;
      foundTP = true;
    }
    // If I have one in-time TP, assign to _bestTP
    else if (cluster.size() == 1) {
      _bestTP = tps_in_chamber[0];
      foundTP = true;
    }

    if (bestTPQuality() > 0) {
      // If I've found _bestTP I remove it from cluster so I only have the ItGhost in it
      cluster.erase(std::remove(cluster.begin(), cluster.end(), _bestTP));
    }

    _itGhosts = std::move(cluster);

  }

  // ########## Segments cluster #############

  // if I found TP cluster-> look for segment around it
  if (foundTP){
    std::copy_if(seg.begin(), seg.end(), std::back_inserter(_SegmentCluster),
                [=](auto& segm){ return  segm.wheel == wh && segm.station == st && segm.sector == sec  && segm.inCluster == false && std::abs(segm.xLoc - _bestTP.xLoc) < xCut; });
    for (auto &s : seg){
      for (auto &Cl_seg : _SegmentCluster){
        if (Cl_seg.index == s.index) {
          s.inCluster = true;
        }
      }
    }
    
    // if I find any look for best qual
    if (_SegmentCluster.size() > 0) {
      foundSeg = true;
      segMatched = true;

      if (_SegmentCluster.size() == 1){
        _bestSeg = _SegmentCluster[0];
      }
      // if there's more than a segment find the higher quality one
      else{
      auto compareQuality = [](Segment& s1, Segment& s2){return s1.nPhiHits < s2.nPhiHits; };
      auto best_seg = std::max_element(_SegmentCluster.begin(), _SegmentCluster.end(), compareQuality);
      _bestSeg = *best_seg;
      }   

      _matchedSeg = _bestSeg; 
    }
  }

  // if I don't have a TPs cluster I try to make a segment cluster
  else {
    if ( seg.size() == 1) {
      _bestSeg = seg[0];
      seg.erase(seg.begin(), seg.end());
      seg[0].inCluster = true;
    }

    else{
      // find segments in same chamber
      std::vector<Segment> segments_in_chamber;  
      std::copy_if(seg.begin(), seg.end(), std::back_inserter(segments_in_chamber),
                [=](auto& segm) { return segm.wheel == wh && segm.station == st && segm.sector == sec && segm.inCluster == false; });

      std::copy_if(segments_in_chamber.begin(), segments_in_chamber.end(), std::back_inserter(_SegmentCluster),
          [=](auto& segm) { return std::abs(segm.xLoc -segments_in_chamber[0].xLoc) < xCut; });

      for (auto &s : seg){
        for (auto &Cl_seg : _SegmentCluster){
          if (Cl_seg.index == s.index) s.inCluster = true;
        }
      }

      // find best quality for segment -> higher n of Phihits
      if (_SegmentCluster.size() > 0) {
        foundSeg = true;
        if (_SegmentCluster.size() == 1){
          _bestSeg = _SegmentCluster[0];
        }
        
        else{
        auto compareQuality = [](Segment& s1, Segment& s2){return s1.nPhiHits < s2.nPhiHits; };
        auto best_seg = std::max_element(_SegmentCluster.begin(), _SegmentCluster.end(), compareQuality);
        _bestSeg = *best_seg;
        }    
      }
    }
  }

  // ########## Digi cluster #############
  std::vector<Digi> DigiClusterSl1;
  std::vector<Digi> DigiClusterSl3;
  
  // if I found TP cluster-> look for digis around it -> se c'è pure segmento guardo intorno al segmento
  if (foundSeg){
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(DigiClusterSl1),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 1 && dig.inCluster == false && std::abs(dig.xLoc - _bestSeg.xLoc) < digiCut;});
    
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(DigiClusterSl3),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 3 && dig.inCluster == false && std::abs(dig.xLoc - _bestSeg.xLoc) < digiCut;});

  }
  else if (foundTP){
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(DigiClusterSl1),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 1 && dig.inCluster == false && std::abs(dig.xLoc - _bestTP.xLoc) < digiCut;});
    
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(DigiClusterSl3),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 3 && dig.inCluster == false && std::abs(dig.xLoc - _bestTP.xLoc) < digiCut;});

  }

  else {
    std::vector<Digi> digis_in_chamber_sl1;
    std::vector<Digi> digis_in_chamber_sl3;

    std::copy_if(digis.begin(), digis.end(), std::back_inserter(digis_in_chamber_sl1),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 1 && dig.inCluster == false;});
    
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(digis_in_chamber_sl3),
                [=](auto& dig) {return dig.wheel == wh && dig.station == st && dig.sector == sec && dig.superlayer == 3 && dig.inCluster == false;});

    if (digis_in_chamber_sl1.size() > 10 || digis_in_chamber_sl3.size() > 10){
      DigiClusterSl1 = digis_in_chamber_sl1[0].FindCluster(digis_in_chamber_sl1, digiCut);
      DigiClusterSl3 = digis_in_chamber_sl3[0].FindCluster(digis_in_chamber_sl3, digiCut);
    }
  }  

  // at least 10 digi per superlayer to build a digi cluster
  if (DigiClusterSl1.size() > 10 ){
    _DigiCluster = DigiClusterSl1;
    foundDigi = true;
    SL1Cluster = true;
  }
  if (DigiClusterSl3.size() > 10 ){
    
    if (foundDigi == false) {
      _DigiCluster = DigiClusterSl3;
      foundDigi = true;
    }
    else _DigiCluster.insert(_DigiCluster.end(), DigiClusterSl3.begin(), DigiClusterSl3.end());
    SL3Cluster = true;
  }

  for (auto &d : digis){
    for (auto Cl_d : _DigiCluster) d.inCluster = true;
  }
  
}

int Cluster::ootSize() const { return _ootGhosts.size(); };
const std::vector<TriggerPrimitive> & Cluster::ootGhosts() const { return _ootGhosts; };

int Cluster::itSize() const { return _itGhosts.size(); };
const std::vector<TriggerPrimitive> & Cluster::itGhosts() const { return _itGhosts; };

int Cluster::tpClusterSize() const { return _ootGhosts.size() + _itGhosts.size() + 1; };


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
    else if ( !segMatched && std::abs( _bestTP.xLoc - muX ) < 10 ){
      muMatchedIndex[0] = iMu;
      muMatchedIndex[1] = muIndex;
      muMatched = true;
    }
  }
}

int Cluster::matchedSegIndex() const {return _matchedSeg.index; };
int Cluster::matchedSegPhiHits() const {return _matchedSeg.nPhiHits; };
const Segment& Cluster::matchedSeg() const {return _matchedSeg; };
int Cluster::segClusterSize() const {return _SegmentCluster.size(); };
const std::vector<Segment> Cluster::segCluster() const {return _SegmentCluster;};

int Cluster::bestSegPhiHits() const {return _bestSeg.nPhiHits; };
const Segment& Cluster::bestSeg() const {return _bestSeg; };

void Cluster::MatchDigi(std::vector<Digi> const& digis, double xCut){
  Segment seg = _bestSeg;

  for (auto const &digi: digis) {
    if (digi.sector != seg.sector && digi.wheel != seg.wheel && digi.station!= seg.station) continue;
    for (int i = 0; i < seg.nPhiHits; ++i){
      if (std::abs(seg.xLoc - digi.xLoc ) < xCut ) {
        digiMatched = true;
        _matchedDigis.emplace_back(digi);
      }
    }
  }
};

const std::vector<Digi> Cluster::matchedDigi() const{ return _matchedDigis; };

const int Cluster::GetNDigi() const {return _matchedDigis.size(); };

const int Cluster::WhichSL() const{
  if (SL1Cluster && SL3Cluster) return 5;
  else if (SL1Cluster && !SL3Cluster) return 1;
  else if (!SL1Cluster && SL3Cluster) return 3;
  else return 0;
};

