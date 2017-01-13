#ifndef PARSERUTIL_H 
#define PARSERUTIL_H

#include <string>

class ParserUtil
{
  public:
    ParserUtil();

    // std::string trimming methods
    static std::string& ltrim(std::string &s);
    static std::string& rtrim(std::string &s);
    static std::string& trim(std::string &s);

    // std::string cleaning, purge trailing comments
    static std::string& stripComment( std::string &s );
    static std::string& toLower( std::string &s );

};

#endif // PARSERUTIL
