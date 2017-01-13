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

    void addKeyValuePair(std::string key, std::string value);

  private:
    std::string m_name = "DEFAULT_CMD_BLOCK";
    std::string m_type = "DEFAULT_TYPE";

    std::map<std::string, int> m_paramInt;
    std::map<std::string, double> m_paramDouble;
    std::map<std::string, std::string> m_paramString;
};

#endif // PARSERCMDBLOCK_H
