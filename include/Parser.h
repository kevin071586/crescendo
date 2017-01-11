#ifndef PARSER_H 
#define PARSER_H

#include <string>

class Parser 
{
  public:
    Parser(); 
    void setFile( std::string m_filename );
    void parse();

  private:
    std::string m_filename;
    
    bool isKeyValuePair( std::string line );
    void getKeyValuePair( std::string line, std::string& key, std::string& value );

    // std::string trimming methods
    static inline std::string& ltrim(std::string &s);
    static inline std::string& rtrim(std::string &s);
    static inline std::string& trim(std::string &s);
};

#endif // PARSER_H
