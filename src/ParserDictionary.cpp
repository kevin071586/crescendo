#include <iostream>
#include <algorithm>
#include <ParserDictionary.h>
#include <ParserUtil.h>

// ParserDictionary, Constructors
// ==================================================================

ParserDictionary::ParserDictionary() 
{
  addCmdBlock("finite element model", true, 1);
  addCmdBlockKey("database name",          keyType::KEY_STRING, true);
  addCmdBlockKey("youngs modulus",    keyType::KEY_DOUBLE, true);
  addCmdBlockKey("density",           keyType::KEY_DOUBLE, true);
  addCmdBlockKey("poissons ratio",    keyType::KEY_DOUBLE, true);

  addCmdBlock("solver parameters", false);

  addCmdBlock("results output", false);
  addCmdBlockKey("database name",          keyType::KEY_STRING, true);
}


// ParamKeyMetadata, Constructors
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
  t_blockMap::iterator it;
  it = m_blockDefs.insert(std::make_pair(name, pbmd));
  m_lastBlockAdded = &(it->second);

  return;
}

// Add a command block
// ==================================================================
void ParserDictionary::addCmdBlock(std::string name, bool isRequired)
{
  ParserUtil::trim(name);

  ParamBlockMetadata pbmd;
  pbmd.setName(name);
  pbmd.setIsRequired(isRequired);

  // Add block metadata to vector
  t_blockMap::iterator it;
  it = m_blockDefs.insert(std::make_pair(name, pbmd));
  m_lastBlockAdded = &(it->second);

  return;
}


// Add a command block
// ==================================================================
void ParserDictionary::addCmdBlockKey(std::string name, keyType type, 
                                      bool isRequired)
{
  ParserUtil::trim(name);

  // Last command block that was added
  ParamBlockMetadata& lastBlock = *m_lastBlockAdded;

  ParamKeyMetadata key;
  key.setName(name);
  key.setIsRequired(isRequired);
  key.setType(type);
  lastBlock.addKey(key);
}

// Get block metadata from dictionary, by block type name
// ==================================================================
ParamBlockMetadata ParserDictionary::getBlock(std::string name)
{
  t_blockMap blockMap = getBlocks();
  t_blockMap::iterator it = blockMap.find(name);
  return it->second;
}

// Get block key metadata from dictionary, by block type name
// ==================================================================
t_keyMap ParserDictionary::getBlockKeys(std::string name)
{
  ParamBlockMetadata block = getBlock(name);
  return block.getKeys();
}

// Check if this block is defined in the dictionary 
// ==================================================================
bool ParserDictionary::isValidBlock(std::string name)
{
  t_blockMap blockMap = getBlocks();
  t_blockMap::iterator it = blockMap.find(name);
  return (it != blockMap.end()) ? true : false;
}


// Check if a block-key pair is defined in the dictionary
// ==================================================================
bool ParserDictionary::isValidBlockKey(std::string blockName, std::string key)
{
  if (isValidBlock(blockName)) {
    t_keyMap keys = getBlockKeys(blockName);
    t_keyMap::iterator it = keys.find(key);
    return (it != keys.end()) ? true : false;

  }
  else {
    return false;
  }
}



// Print out available command syntax
// ==================================================================
void ParserDictionary::print()
{
  t_blockMap::iterator blk_it;
  for (blk_it=m_blockDefs.begin(); blk_it != m_blockDefs.end(); ++blk_it) {
    ParamBlockMetadata blk = blk_it->second;
    std::cout << "begin " << blk.getName() << std::endl;
    
    t_keyMap keyMap = blk.getKeys();
    t_keyMap::iterator it;

    for (it = keyMap.begin(); it != keyMap.end(); ++it) {
      ParamKeyMetadata key = it->second;
      std::string keyTypeStr;

      switch (key.getKeyType()) {
        case keyType::KEY_INT:
          keyTypeStr = "<int>";
          break;
        case keyType::KEY_DOUBLE:
          keyTypeStr = "<real>";
          break;
        case keyType::KEY_STRING:
          keyTypeStr = "<string>";
          break;
      }
      std::cout << "  " << key.getName() << " = "
                << keyTypeStr
                << std::endl;
    }

    // end command block
    std::cout << "end\n" << std::endl;
  }
  return; 
}
