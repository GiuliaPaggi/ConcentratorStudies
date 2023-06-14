#ifndef TestAnalyser_h
#define TestAnalyser_h

#include "include/Analyser.h"
#include <TH2.h>
#include <TEfficiency.h>


class TestAnalyser : public Analyser {
private:
   std::map<std::string, TH1*> m_plots;
   std::map<std::string, TH2*> m_2Dplots;
   std::map<std::string, TEfficiency*> m_effs;
   std::vector<std::string> tags;

   std::map<std::string, double> m_counters;
   double nClustersGhosts{};
   double ooTHQCount{};
   double nClusters{};

   void ClusterAnalisis(std::vector<Cluster> CLtoAnalize, std::string CLtype, std::vector<Segment> Segments);
   void DefinePlot();
public:
  TestAnalyser(std::string file) : Analyser{file} {}
  void Loop() final;
};

#endif
