#ifndef Analyser_h
#define Analyser_h

#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>

#include <map>
#include <vector>

#include "include/AnalyserBase.h"

class Analyser : public AnalyserBase {
 private:
  std::map<std::string, TH1 *> m_plots;
  std::map<std::string, TH2 *> m_2Dplots;
  std::map<std::string, TEfficiency *> m_effs;
  std::vector<std::string> tags;

  std::map<std::string, double> m_counters;
  double nClustersGhosts{};
  double ooTHQCount{};
  double nClusters{};

  void ClusterAnalisis(const std::vector<Cluster> &, const std::string &, const std::vector<Segment> &);
  void DefinePlot();

 public:
  Analyser(std::string file) : AnalyserBase{file} {}
  ~Analyser() final {}

  void Loop() final;
};

#endif
