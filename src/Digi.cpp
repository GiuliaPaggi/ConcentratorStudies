#include "include/Digi.h"

#include <algorithm>
#include <numeric>

#include "include/Geometry.h"

Digi::Digi(int i, int wh, int sec, int stat, int SL, int L, int w, double t)
    : wheel{wh}, sector{GEOM.sectorInRange(sec)}, station{stat}, superlayer{SL}, layer{L}, wire{w}, time{t} {
  index = i;

  const auto firstWire{GEOM.firstWire(wheel, station, sector, superlayer, layer)};
  const auto xFirstWire{GEOM.xFirstWire(wheel, station, sector, superlayer, layer)};

  xLoc = xFirstWire + (wire - firstWire) * GEOM.CELL_WIDTH;
}

std::vector<Digi> Digi::findCluster(const std::vector<Digi>& digis, double cut) const {
  std::vector<Digi> cluster;

  // find digis close to the one that calls the function
  std::copy_if(digis.begin(), digis.end(), std::back_inserter(cluster),
               [=](auto& d) { return std::abs(xLoc - d.xLoc) < cut; });

  // if there are more digis in the chamber than in the cluster, compute mean_xLoc
  // as the mean value of in the cluster just built and use it to look for close
  // digis again
  if (digis.size() > cluster.size()) {
    auto mean_xLoc = std::accumulate(cluster.begin(), cluster.end(), 0.0,
                                     [](const auto sum, const auto digi) { return sum + digi.xLoc; }) /
                     cluster.size();

    cluster.clear();
    std::copy_if(digis.begin(), digis.end(), std::back_inserter(cluster),
                 [=](auto& d) { return std::abs(mean_xLoc - d.xLoc) < cut; });
  }

  return cluster;
}