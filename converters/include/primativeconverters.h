#ifndef H_PRIMATIVECONVERTERS
#define H_PRIMATIVECONVERTERS
#include "converters.h"


namespace PrimativeConverters{
static bool is_digits(std::string &s)
  {
  if (strcmp(&s[0], "-") == 0 && s.size() > 1)
    {
    std::string t = s.substr(1, s.size());
    return(t.find_first_not_of("0123456789")== std::string::npos  );
    }
  else  return(s.find_first_not_of("0123456789")== std::string::npos  );
  }


static void IntConverter(std::vector<std::string>&s, void* data){
  std::vector<int>* returnInts = static_cast<std::vector<int>* >(data);
  for (size_t i=0; i<s.size(); i++)
  {
  if (!is_digits(s[i])) {std::cout << s[i] << " is not an int!" << std::endl; throw std::exception();}
  (*returnInts).push_back(strtol(s[i].c_str(), NULL ,10));
  }
  return;

}
static void BoolConverter(std::vector<std::string>&s, void* data){
  std::vector<bool>* returnBools = static_cast<std::vector<bool>* >(data);//recast back to correct type
   //must have on or off, does not support --flag  
  for (size_t i=0; i<s.size(); i++)
  {
  if (strcmp(s[i].c_str(), "on") == 0 || strcmp(s[i].c_str(), "On") == 0 || strcmp(s[i].c_str(), "ON") == 0) returnBools->push_back(true);
  else if (strcmp(s[i].c_str(), "off") == 0 || strcmp(s[i].c_str(), "Off") == 0 || strcmp(s[i].c_str(), "OFF") == 0) returnBools->push_back(false);
  //#define DRACOBOOL
  #ifdef DRACOBOOL
  else returnBools->push_back(false);
  #endif
  else {std::cout << "Expecting on or off for boolean variable. Check PrimativeConverters.cpp for details." << std::endl; throw std::exception();}
   }
  return;
}

static void DoubleConverter(std::vector<std::string>&s, void* data){
  std::vector<double>* returnDoubles = static_cast<std::vector<double>* >(data);
  for (size_t i=0; i<s.size(); i++)
  {
//    if (!is_digits(s[i])) {std::cout << s[i] << " is not a double!" << std::endl; throw std::exception();}
    returnDoubles->push_back(strtod(s[i].c_str(), NULL));
  }
  return;
  
}
static void StringConverter(std::vector<std::string>&s, void* data){
std::vector<std::string>* returnString = static_cast<std::vector<std::string>* >(data);
for (size_t i=0; i<s.size(); i++)
  {
  returnString->push_back(s[i]);//a little pointless to make this class, but why not
  }
return;
}


} 
REGISTER_CLASS(iint, "INT",&PrimativeConverters::IntConverter);
REGISTER_CLASS(ibool, "BOOL",&PrimativeConverters::BoolConverter);
REGISTER_CLASS(ireal, "REAL",&PrimativeConverters::DoubleConverter);
REGISTER_CLASS(istring, "STRING",&PrimativeConverters::StringConverter);
#endif
