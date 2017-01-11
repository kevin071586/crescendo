#ifndef PARSER_H 
#define PARSER_H

#include <string>
#include <map>

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
    
    bool isKeyValuePair( std::string line );
    void getKeyValuePair( std::string line, std::string& key, std::string& value );

    // std::string trimming methods
    static inline std::string& ltrim(std::string &s);
    static inline std::string& rtrim(std::string &s);
    static inline std::string& trim(std::string &s);

    // std::string cleaning, purge trailing comments
    static inline std::string& stripComment( std::string &s );
    static inline std::string& toLower( std::string &s );
};

#endif // PARSER_H
