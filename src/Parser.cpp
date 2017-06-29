#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <algorithm>

#include <Parser.h>
#include <ParserCmdBlock.h>
#include <ParserDictionary.h>

// Default constructor
// ==================================================================
Parser::Parser() 
{
  ParserCmdBlock main_block;
  m_cmdBlocks.push_back(main_block);
}

Parser::Parser( std::string filename ) 
{
  ParserCmdBlock main_block;
  m_cmdBlocks.push_back(main_block);
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

  // Read each line from the input file
  while (std::getline(inFile, line)) {

    // Read key-value pair
    if (isKeyValuePair(line)) {
      std::string key, value;
      getKeyValuePair(line, key, value);

      // Add key-value pairs to current block
      ParserCmdBlock& currBlock = m_cmdBlocks.back();
      if (m_dictionary.isValidBlockKey(currBlock.getType(), key)) {
        currBlock.addKeyValuePair(key,value);
      }
    }

    // Read command blocks
    if (isCmdBlock(line)) {
      ++nestingCount;

      int posBlkStart = line.find_first_of(" ") + 1;
      int posBlkEnd = line.find_last_of(" ");
      std::string blockType = line.substr(posBlkStart, posBlkEnd - 6);
      std::string blockName = line.substr(posBlkEnd+1, line.length());
      m_parserUtil.trim(blockName);

      ParserCmdBlock cmdBlock(blockType, blockName);
      m_dictionary.isValidBlock(blockType);
      m_cmdBlocks.push_back(cmdBlock);
    }

    // Leaving a command block
    if (m_parserUtil.toLower(line).find("end") == 0) {
      --nestingCount;
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

// Get key-value pair from a line in the input deck, to be 
// assigned to a ParserCmdBlock object
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


// Get a ParserCmdBlock object by block type
// TODO: This will return the first block of the type it finds.. if
// there are multiple blocks of this type defined it will not work right.
// ==================================================================
ParserCmdBlock Parser::getCmdBlock(std::string type) 
{
  ParserCmdBlock cmdBlock;
  for (int i=0; i < m_cmdBlocks.size(); ++i) {
    cmdBlock = m_cmdBlocks[i];
    if (cmdBlock.getType().compare(type) == 0) { 
      break; 
    }
  }
  return cmdBlock;
}
