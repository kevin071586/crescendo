#ifndef PARSERCMDBLOCK_H 
#define PARSERCMDBLOCK_H

#include <string>
#include <map>
#include <iostream>

class ParserCmdBlock 
{
  public:
    ParserCmdBlock();
    ParserCmdBlock(std::string type, std::string name);

    // Get methods
    std::string getName() {return m_name;}
    std::string getType() {return m_type;}

    // Data access methods
    int getFieldInt( std::string key );
    double getFieldDouble( std::string key );
    std::string getFieldString( std::string key );

    // Data assignment methods
    void addKeyValuePair(std::string key, std::string value);

    std::map<std::string, std::string> m_paramVal;
  private:
    std::string m_name = "DEFAULT_CMD_BLOCK";
    std::string m_type = "DEFAULT_TYPE";
};

#endif // PARSERCMDBLOCK_H
