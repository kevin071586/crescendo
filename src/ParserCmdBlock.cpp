#include <ParserCmdBlock.h>

// Constructor
// ==================================================================
ParserCmdBlock::ParserCmdBlock() 
{
}

ParserCmdBlock::ParserCmdBlock(std::string type, std::string name) 
{
  m_name = name;
  m_type = type;
}

// Add a key value pair 
// ==================================================================
void ParserCmdBlock::addKeyValuePair(std::string key, std::string value)
{
  // 1. Check, is this key defined in block of this type?
  // 2. If so, determine the type of key this is (INT, STR, DBL)
  // 3. Use that information to assign the key to the correct map.
  //    - Usage in the code requires prior knowledge of how the
  //      key type (INT, STR, DBL) is defined, but that is fine..
  //      can't expect the code to magically figure that out.
  //    - Proposed usage something like:
  //        cmdBlock = parser.getCmdBlock("block type")
  //        cmdBlock.getParamInt("keyName") maybe?


  // Check if key exists in this block
  // if (std::find(v.begin(), v.end(), "abc") != v.end()) {
  //   std::cout << "Block name: " << m_name << std::endl; 
  //   std::cout << "Block type: " << m_type << std::endl; 
  // }
  // else {
  // }
  
  m_paramString[key] = value;

  std::cout << "Adding to block " << m_name << ": " 
            << "Key: " << key << ", Type: " << m_type << std::endl;

  return;
}
