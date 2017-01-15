#include <ParserCmdBlock.h>
#include <ParserUtil.h>

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

  m_paramVal[key] = value;

  std::cout << "Adding to block " << m_name << ": " 
            << "Key: " << key << ", Value: " << value << std::endl;
  return;
}

// Access methods, return fields by type 
// ==================================================================
int ParserCmdBlock::getFieldInt(std::string key)
{
  ParserUtil::toLower(key);
  try {
    return std::stoi(m_paramVal[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return 0;
}

double ParserCmdBlock::getFieldDouble(std::string key)
{
  ParserUtil::toLower(key);
  try {
    return std::stod(m_paramVal[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return 0.0;
}

std::string ParserCmdBlock::getFieldString(std::string key)
{
  ParserUtil::toLower(key);
  return m_paramVal[key];
}
