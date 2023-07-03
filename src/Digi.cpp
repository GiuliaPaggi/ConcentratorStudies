#include "include/Digi.h"

#include <algorithm>
Digi::Digi(const Geometry& geom, int i, int wh, int sec, int stat, int SL, int L, int w, double t)
    : wheel{wh}, sector{sec}, station{stat}, superlayer{SL}, layer{L}, wire{w}, time{t} {
  index = i;

  if (sector == 13)
    sector = 4;
  else if (sector == 14)
    sector = 10;

  const auto firstWire{geom.firstWire(wheel, station, sector, superlayer, layer)};
  const auto xFirstWire{geom.xFirstWire(wheel, station, sector, superlayer, layer)};

  xLoc = xFirstWire + (wire - firstWire) * geom.CELL_WIDTH;
}

std::vector<Digi> Digi::findCluster(std::vector<Digi> digis, double cut) const {
  std::vector<Digi> digiCluster;

  // find digis close to the one that calls the function
  std::copy_if(digis.begin(), digis.end(), std::back_inserter(digiCluster),
               [=](auto& d) { return std::abs(xLoc - d.xLoc) < cut; });

  // if there are more digis in the chamber-> find mean value of xLoc in the cluster we have for now and use that to
  // look for close digis
  if (digis.size() > digiCluster.size()) {
    double mean_xLoc = xLoc;
    for (auto d : digis) mean_xLoc += d.xLoc;  // CB std::accumulate
    mean_xLoc = mean_xLoc / digiCluster.size();
    // clear cluster to fill it again with new center
    digiCluster.clear();
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(digiCluster),
                 [=](auto& d) { return std::abs(mean_xLoc - d.xLoc) < cut; });
  }

  return digiCluster;
}