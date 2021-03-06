#include <iostream>
#include <string>
#include "Parser.h"
#include <ParserUtil.h>
#include <ParserDictionary.h>
#include <ParserCmdBlock.h>

int main(int argc, char** argv) {

  if (argc == 1) {
    std::cout << "Error: Must supply an input file" << std::endl;
    return 0;
  }

  std::string fileName(argv[1]);
  std::cout << "Input file: " << fileName << std::endl;

  Parser parser(fileName);
  try {
    //parser.setFile(fileName);
    parser.parse();
  }
  catch (int e) {
    std::cout << "A parsing error occurred." << std::endl;
    return 1;
  }

  // std::cout << "Input Values" << std::endl;
  // std::cout << " Param1:    " << parser.getFieldInt("Param1") << std::endl;
  // std::cout << " Param2:    " << parser.getFieldDouble("Param2") << std::endl;
  // std::cout << " Param3:    " << parser.getFieldString("param3") << std::endl;
  // std::cout << " FakeParam:    " << parser.getFieldString("fakeparam") << std::endl;

  // Try to access a parameter that wasnt defined
  // std::cout << "Testings that should fail: " << std::endl;
  // parser.getFieldInt("doesntExist");
  // parser.getFieldDouble("doesntExist");
  // parser.getFieldString("doesntExist");
  // std::cout << std::endl;

  // Build a ParserUtil object and use a method to make sure it doesnt fail
  ParserUtil m_parserUtil;
  std::string testString = "  Trim left and right  ";
  m_parserUtil.trim(testString);

  // Dump the available output
  ParserDictionary dict;
  dict.print();

  // Print out some values from the input deck
  ParserCmdBlock cmdBlock;
  cmdBlock = parser.getCmdBlock("finite element model");
  double val = cmdBlock.getFieldDouble("youngs modulus");
  std::cout << "  youngs modulus: " << val << std::endl;

  return 0;
}
