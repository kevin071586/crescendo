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

#include <Parser.h>
#include <ParserCmdBlock.h>

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
  
  int nestingCount = 0;
  std::vector<ParserCmdBlock> cmdBlocks;

  // Read each line from the input file
  while (std::getline(inFile, line)) {

    // Read key-value pair
    if (isKeyValuePair(line)) {
      std::string key, value;
      getKeyValuePair( line, key, value );
      m_inputParam[key] = value;
    }

    // Read command blocks
    if (isCmdBlock(line)) {
      ++nestingCount;
      ParserCmdBlock cmdBlock;
      cmdBlocks.push_back(cmdBlock);

      int pos = line.find_last_of(" ");
      std::string blockName = line.substr( pos+1, line.length() );
      std::cout << "Found command block: " << blockName << std::endl;
    }

    // Leaving a command block
    if (m_parserUtil.toLower(line).find("end") == 0) {
      --nestingCount;
      std::cout << "Leaving a command block." << std::endl;
    }

  }

  std::cout << "Total number of command blocks: " << cmdBlocks.size() 
    << std::endl;
  return;
}


// Access methods, return fields by type 
// ==================================================================
int Parser::getFieldInt( std::string key )
{
  m_parserUtil.toLower(key);
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
  m_parserUtil.toLower(key);
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
  m_parserUtil.toLower(key);
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
  m_parserUtil.toLower(key);
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
  m_parserUtil.toLower(key);
  return m_inputParam[key];
}
std::string Parser::getFieldString( std::string key, std::string defaultValue )
{
  m_parserUtil.toLower(key);
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


// Check if string is in begin <block_type> <name> format
// ==================================================================
bool Parser::isCmdBlock( std::string line )
{
  m_parserUtil.ltrim(line);
  if (line.find("begin") == 0){
    return true;
  }
  return false;
}



// Get key-value pair
// ==================================================================
void Parser::getKeyValuePair( std::string line, std::string& key, std::string& value )
{
  int pos = line.find("=");
  key = line.substr(0, pos);
  value = line.substr(pos+1, line.length());
      //TODO: Make sure pos+1 isnt the end of the string.

  m_parserUtil.trim(key);
  m_parserUtil.stripComment(key);
  m_parserUtil.toLower(key);

  m_parserUtil.trim(value);
  m_parserUtil.stripComment(value);
  return;
}


