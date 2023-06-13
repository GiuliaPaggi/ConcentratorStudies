#ifndef TestAnalyser_h
#define TestAnalyser_h

#include "include/Analyser.h"

class TestAnalyser : public Analyser {
public:
  TestAnalyser(std::string file) : Analyser{file} {}
  void Loop() final;
};

#endif
