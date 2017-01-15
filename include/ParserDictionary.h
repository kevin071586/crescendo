#ifndef PARSERDICTIONARY_H
#define PARSERDICTIONARY_H

#include <string>
#include <map>
#include <vector>

// Forward declarations for typedefs
class ParamBlockMetadata;
class ParamKeyMetadata;
typedef std::multimap<std::string, ParamBlockMetadata> t_blockMap;
typedef std::multimap<std::string, ParamKeyMetadata> t_keyMap;

// Command Block-Key Metadata
// ==================================================================
class ParamKeyMetadata
{
  public:
    enum keyType
    {
      KEY_INT,
      KEY_DOUBLE,
      KEY_STRING
    };

    ParamKeyMetadata();
    void setName(std::string name) { m_name = name; }
    void setIsRequired(bool isReq) { m_isRequired = isReq; }
    void setType(keyType type) { m_keyType = type; }

    std::string getName() { return m_name; }
    keyType getKeyType() { return m_keyType; }
    
  private:
    std::string m_name;
    keyType m_keyType;
    bool m_isRequired = true;
};


// Command Block Metadata
// ==================================================================
class ParamBlockMetadata
{
  public:
    ParamBlockMetadata();


    // Set methods
    void setName(std::string name) { m_name = name; }
    void setIsRequired(bool isReq) { m_isRequired = isReq; }
    void setMaxOccur(int maxOccur) { m_maxOccur = maxOccur; }

    // Get methods
    std::string getName() { return m_name; }
    t_keyMap getKeys() { return m_keyMap; }

    // Add keys
    void addKey(ParamKeyMetadata key) { 
      m_keyMap.insert(std::make_pair(key.getName(), key)); 
    }

  private:
    std::string m_name = "undefined_command_block_name";
    bool m_isRequired = true;
    int m_maxOccur = 999; // essentially unlimited
    t_keyMap m_keyMap;
};

// Parser Dictionary
// ==================================================================
class ParserDictionary 
{
  public:
    ParserDictionary();
    void print();
    t_blockMap getBlocks() { return m_blockDefs; }
    t_keyMap getBlockKeys(std::string name);
    ParamBlockMetadata getBlock(std::string name);

    bool isValidBlock(std::string name);
    bool isValidBlockKey(std::string blockName, std::string key);
  
  private:
    t_blockMap m_blockDefs;
    ParamBlockMetadata* m_lastBlockAdded;

    typedef ParamKeyMetadata::keyType keyType;

    void addCmdBlock(std::string name, bool isRequired, int maxOccur);
    void addCmdBlock(std::string name, bool isRequired);
    void addCmdBlockKey(std::string name, keyType type, bool isRequired);

    int getBlockMetadataIndex(std::string name);
};


#endif // PARSERDICTIONARY_H
