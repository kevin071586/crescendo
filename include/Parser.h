#ifndef PARSER_H 
#define PARSER_H

#include <string>
#include <map>
#include <vector>
#include <ParserCmdBlock.h>
#include <ParserUtil.h>
#include <ParserDictionary.h>

class Parser 
{
  public:
    // Constructors
    Parser(); 
    Parser( std::string filename ); 

    // Basic use methods
    void setFile(std::string filename);
    void parse();
    ParserCmdBlock getCmdBlock(std::string type);

  private:
    std::string m_filename;
    ParserDictionary m_dictionary;

    ParserUtil m_parserUtil;
    std::vector<ParserCmdBlock> m_cmdBlocks;
    
    bool isKeyValuePair( std::string line );
    bool isCmdBlock( std::string line );
    void getKeyValuePair( std::string line, std::string& key, std::string& value );

};

#endif // PARSER_H
