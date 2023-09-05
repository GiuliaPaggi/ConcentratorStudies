#ifndef Analyser_h
#define Analyser_h

#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <map>
#include <vector>

#include "include/AnalyserBase.h"

class Analyser : public AnalyserBase {
 private:
  std::map<std::string, TH1 *> m_plots;
  std::map<std::string, TH2 *> m_2Dplots; // CB in theory not needed
  std::map<std::string, TEfficiency *> m_effs;
  std::map<std::string, TProfile2D *> m_2Dprofiles; 
  std::vector<std::string> tags;

  std::map<std::string, double> m_counters;
  // double nClustersGhosts{};
  // double ooTHQCount{}; // CB not used
  // double nClusters{};

  void ClusterAnalisis(const std::vector<Cluster> &, const std::string &, const std::vector<Segment> &);
  void DefinePlot();
  void FillEfficiency(const std::string &, const std::string &, const int &, const int &, const double &, const Cluster &);
  void FillEfficiency2D(const std::string &, const std::string &, const int &, const int &, const double &, const double &, const Cluster &);
  void FillGhostRatio(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster);
  void FillGhostProfile(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster);
  void FillBackground(const std::string &typeStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster);

 public:
  Analyser(std::string file) : AnalyserBase{file} {}
  ~Analyser() final {}

  void Loop() final;
};

#endif
