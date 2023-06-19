
#include "include/Analyser.h"

#include <iostream>
#include <filesystem>

int main(int argc, char* argv[] ) {
  
  if(argc != 2 || !std::filesystem::is_regular_file(argv[1])) {
    std::cout << "[" << argv[0] << "] Usage:" << argv[0] << " <ROOT_INPUT_DTTREE>\n";
    return EXIT_FAILURE;
  }

  Analyser analysis{argv[1]};
  analysis.Loop();
};
