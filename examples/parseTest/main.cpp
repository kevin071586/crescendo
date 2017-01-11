#include <iostream>
#include <string>
#include "Parser.h"

int main(int argc, char** argv) {

  if (argc == 1) {
    std::cout << "Error: Must supply an input file" << std::endl;
    return 0;
  }

  std::string fileName(argv[1]);
  std::cout << "Input file: " << fileName << std::endl;

  Parser parser;
  parser.setFile(fileName);
  parser.parse();

  return 0;
}
