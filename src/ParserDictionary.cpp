#include <iostream>
#include <ParserDictionary.h>
#include <ParserUtil.h>

// ParserDictionary, Constructors
// ==================================================================

ParserDictionary::ParserDictionary() 
{
  addCmdBlock("finite element model", true, 1);
  addCmdBlockKey("database",          keyType::KEY_STRING, true);
  addCmdBlockKey("youngs modulus",    keyType::KEY_DOUBLE, true);
  addCmdBlockKey("density",           keyType::KEY_DOUBLE, true);
  addCmdBlockKey("poissons ratio",    keyType::KEY_DOUBLE, true);

  addCmdBlock("solver parameters", false);
}


// ParamBlockMetadata, Constructors
// ==================================================================
ParamKeyMetadata::ParamKeyMetadata()
{
}

// ParamBlockMetadata, Constructors
// ==================================================================
ParamBlockMetadata::ParamBlockMetadata()
{
}

// Add a command block
// ==================================================================
void ParserDictionary::addCmdBlock(std::string name,
                                   bool isRequired, int maxOccur)
{
  ParserUtil::trim(name);

  ParamBlockMetadata pbmd;
  pbmd.setName(name);
  pbmd.setIsRequired(isRequired);
  pbmd.setMaxOccur(maxOccur);

  // Add block metadata to vector
  m_blockDefs.push_back(pbmd);

  return;
}

// Add a command block, unlimited occurences
// ==================================================================
void ParserDictionary::addCmdBlock(std::string name, bool isRequired)
{
  ParserUtil::trim(name);

  ParamBlockMetadata pbmd;
  pbmd.setName(name);
  pbmd.setIsRequired(isRequired);

  // Add block metadata to vector
  m_blockDefs.push_back(pbmd);

  return;
}


// Add a command block, unlimited occurences
// ==================================================================
void ParserDictionary::addCmdBlockKey(std::string name, keyType type, 
                                      bool isRequired)
{
  ParserUtil::trim(name);

  // Last command block that was added
  ParamBlockMetadata& lastBlock = m_blockDefs.back();

  ParamKeyMetadata key;
  key.setName(name);
  key.setIsRequired(isRequired);
  key.setType(type);

  switch (type) {
    case keyType::KEY_INT:
      lastBlock.addKeyInt(key);
      break;

    case keyType::KEY_DOUBLE:
      lastBlock.addKeyDouble(key);
      break;

    case keyType::KEY_STRING:
      lastBlock.addKeyString(key);
      break;

    // no default case
  }
}

// Print out available command syntax
// ==================================================================
void ParserDictionary::Print()
{
  for (int i=0; i < m_blockDefs.size(); ++i) {
    ParamBlockMetadata blk = m_blockDefs[i];
    std::cout << "begin " << blk.getName() << std::endl;

    // print int keys
    std::vector<ParamKeyMetadata> keys;
    keys = blk.getKeyInt();
    for (int j=0; j < keys.size(); ++j) {
      std::cout << "  " << keys[j].getName() << " = <int>"
                << std::endl;
    }

    // print double keys
    keys = blk.getKeyDouble();
    for (int j=0; j < keys.size(); ++j) {
      std::cout << "  " << keys[j].getName() << " = <real>"
                << std::endl;
    }

    // print string keys
    keys = blk.getKeyString();
    for (int j=0; j < keys.size(); ++j) {
      std::cout << "  " << keys[j].getName() << " = <string>"
                << std::endl;
    }

    // end command block
    std::cout << "end\n" << std::endl;
  }
  return; 
}
