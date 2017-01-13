#ifndef PARSER_H 
#define PARSER_H

#include <string>
#include <map>
#include <ParserCmdBlock.h>
#include <ParserUtil.h>

class Parser 
{
  public:
    // Constructors
    Parser(); 
    Parser( std::string filename ); 

    // Basic use methods
    void setFile( std::string filename );
    void parse();

    // Access methods
    int getFieldInt( std::string key );
    int getFieldInt( std::string key, int defaultValue );
    double getFieldDouble( std::string key );
    double getFieldDouble( std::string key, double defaultValue );
    std::string getFieldString( std::string key );
    std::string getFieldString( std::string key, std::string defaultValue );

  private:
    std::string m_filename;
    std::map<std::string, std::string> m_inputParam;

    ParserUtil m_parserUtil;
    ParserCmdBlock m_mainBlock;
    
    bool isKeyValuePair( std::string line );
    bool isCmdBlock( std::string line );
    void getKeyValuePair( std::string line, std::string& key, std::string& value );

};

#endif // PARSER_H
