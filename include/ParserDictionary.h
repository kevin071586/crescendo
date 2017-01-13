#ifndef PARSERDICTIONARY_H
#define PARSERDICTIONARY_H

#include <string>
#include <map>
#include <vector>


// Command Block-Key Metadata
// ==================================================================
class ParamKeyMetadata
{
  public:
    ParamKeyMetadata();
    void setName(std::string name) { m_name = name; }
    void setIsRequired(bool isReq) { m_isRequired = isReq; }

    std::string getName() { return m_name; }

  private:
    std::string m_name;
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
    std::vector<ParamKeyMetadata> getKeyInt() { return keyInt; }
    std::vector<ParamKeyMetadata> getKeyDouble() { return keyDouble; }
    std::vector<ParamKeyMetadata> getKeyString() { return keyString; }

    void addKeyInt(ParamKeyMetadata key) { keyInt.push_back(key); }
    void addKeyDouble(ParamKeyMetadata key) { keyDouble.push_back(key); }
    void addKeyString(ParamKeyMetadata key) { keyString.push_back(key); }

  private:
    std::string m_name;
    bool m_isRequired = true;
    int m_maxOccur = 999; // essentially unlimited
    
    std::vector<ParamKeyMetadata> keyInt;
    std::vector<ParamKeyMetadata> keyDouble;
    std::vector<ParamKeyMetadata> keyString;

};

// Parser Dictionary
// ==================================================================
class ParserDictionary 
{
  public:
    ParserDictionary();
    void Print();
  
  private:
    std::vector<ParamBlockMetadata> m_blockDefs;

    enum keyType
    {
      KEY_INT,
      KEY_DOUBLE,
      KEY_STRING
    };

    void addCmdBlock(std::string name, bool isRequired, int maxOccur);
    void addCmdBlock(std::string name, bool isRequired);
    void addCmdBlockKey(std::string name, keyType type, bool isRequired);
};


#endif // PARSERDICTIONARY_H
