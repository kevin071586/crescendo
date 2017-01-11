#include <Parser.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

// Includes for trim algorithms
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

// Default constructor
// ==================================================================
Parser::Parser() 
{
}

// Set input file
// ==================================================================
void Parser::setFile( std::string filename )
{
  m_filename = filename;
  return;
}



// Read the file, line by line.
// ==================================================================
void Parser::parse()
{
  std::fstream inFile( m_filename );
  std::string line;

  while (std::getline(inFile, line)) {
    std::cout << trim(line) << std::endl;

    if (isKeyValuePair(line)) {
      std::string key, value;
      getKeyValuePair( line, key, value );
      std::cout << "Key: " << key << ", Value: " << value << std::endl;
    }
  }
  return;
}



// Check if string is in key-value pair format
// ==================================================================
bool Parser::isKeyValuePair( std::string line )
{
  if (line.find("=") == -1) {
    return false;
  }
  return true;
}



// Get key-value pair
// ==================================================================
void Parser::getKeyValuePair( std::string line, std::string& key, std::string& value )
{
  int pos = line.find("=");
  key = line.substr(0, pos);
  value = line.substr(pos+1, line.length());
      //TODO: Make sure pos+1 isnt the end of the string.

  trim(key);
  trim(value);
  return;
}


// Trim from start
// ==================================================================
std::string& Parser::ltrim( std::string &s )
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
          std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}



// Trim from end
// ==================================================================
std::string& Parser::rtrim( std::string &s )
{
  s.erase(std::find_if(s.rbegin(), s.rend(),
          std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}



// Trim from both ends
// ==================================================================
std::string& Parser::trim( std::string &s )
{
  return ltrim(rtrim(s));
}


