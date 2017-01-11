#include <Parser.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

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

Parser::Parser( std::string filename ) 
{
  m_filename = filename;
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
  
  // Read each line from the input file
  while (std::getline(inFile, line)) {
    if (isKeyValuePair(line)) {
      std::string key, value;
      getKeyValuePair( line, key, value );
      m_inputParam[key] = value;
    }
  }
  return;
}


// Access methods, return fields by type 
// ==================================================================
int Parser::getFieldInt( std::string key )
{
  toLower(key);
  try {
    return std::stoi(m_inputParam[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return 0;
}

int Parser::getFieldInt( std::string key, int defaultValue )
{
  toLower(key);
  try {
    return std::stoi(m_inputParam[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return defaultValue;
}

double Parser::getFieldDouble( std::string key )
{
  toLower(key);
  try {
    return std::stod(m_inputParam[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return 0.0;
}

double Parser::getFieldDouble( std::string key, double defaultValue )
{
  toLower(key);
  try {
    return std::stod(m_inputParam[key]);
  }
  catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument for parameter " << key
              << ":" << ia.what() << '\n';
  }
  return defaultValue;
}

std::string Parser::getFieldString( std::string key )
{
  toLower(key);
  return m_inputParam[key];
}
std::string Parser::getFieldString( std::string key, std::string defaultValue )
{
  toLower(key);
  if ( m_inputParam.find(key) == m_inputParam.end() ) {
    return defaultValue;
  } 
  else {
    return m_inputParam[key];
  }
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
  stripComment(key);
  toLower(key);

  trim(value);
  stripComment(value);
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


// Strip trailing comments 
// ==================================================================
std::string& Parser::stripComment( std::string &s )
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
std::string& Parser::toLower( std::string &s ) {
  std::transform(s.begin(), s.end(), s.begin(), tolower);
  return s;
}
