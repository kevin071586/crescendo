#include <ParserUtil.h>

#include <string>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

// Default constructor
// ==================================================================

ParserUtil::ParserUtil() 
{
}


// Trim from start
// ==================================================================
std::string& ParserUtil::ltrim( std::string &s )
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
          std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}


// Trim from end
// ==================================================================
std::string& ParserUtil::rtrim( std::string &s )
{
  s.erase(std::find_if(s.rbegin(), s.rend(),
          std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}


// Trim from both ends
// ==================================================================
std::string& ParserUtil::trim( std::string &s )
{
  return ltrim(rtrim(s));
}


// Strip trailing comments 
// ==================================================================
std::string& ParserUtil::stripComment( std::string &s )
{
  int pos = s.find("#");
  if (pos == -1) {
    return s;
  }
  else {
    return s = s.substr(0, pos);
  }
}

 
// Set string to lowercase for comparisons 
// ==================================================================
std::string& ParserUtil::toLower( std::string &s ) {
  std::transform(s.begin(), s.end(), s.begin(), tolower);
  return s;
}
