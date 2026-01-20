#ifndef H_PARSEINPUT
#define H_PARSEINPUT

#include <vector> //obvious
#include <string> //obvious
#include <iostream>//for errors
#include <cstdlib> //for strtod and strtol, maybe will be phased out in the future.
#include "converters/include/converters.h" //vector of type conversion functions
#include "common/include/types.h"

//global definitions
typedef std::pair<std::string, std::vector<std::string> > ArgumentPair;
//operators for comparisions and checking 


//helper functions for argument pair that don't live in a class
//but I don't want to give global scope
//contains:
//printArgumentPair -> prints all info
//checkDoubleHyphen -> checks if -- term
//parseArgvToPair -> parses user input into usable ArgumentPair data 
namespace ParsingHelperNS {
void printArgumentPair(ArgumentPair argp);
bool checkDoubleHyphen(std::string const & argument);
//runs through all the arguments and parses them out into pairs
//returns a vectors of pairs, with each pair being a --term and all its constituent arguments.
std::vector<ArgumentPair> parseArgvToPair(int argcount, const char* arguments_[], bool cmd=true);
}//end namespace

//class
class ExpectedInput {
private:
public:
std::string parameter;
std::string errorMessage;
//by default true;
bool mandatory;
bool strict; //if false, you can include UP to expectedNumInputs, if true, you need exactly that many
size_t expectedNumInputs;
//lets us do some light checking for type, and convert
std::string type;
void PrintExpected();

//constructors
ExpectedInput(){};
ExpectedInput(std::string parameter_, int arguments_, bool mandatory_, std::string type_, std::string errorMessage_, bool strict_=true);
//const ExpectedInput arr[] = {ExpectedInput("file", 3, true, ParsingHelperNS::ISTRING, "Expecting 3 File Names")};

//I do not support different data types for a single --term, just make a new term in that case
//std::vector<int> InputToInt();
//std::vector<Real> InputToReal();
//std::vector<bool> InputToBool(std::string trueString);
//everything is already a string :)
//anything more complicated should probably be done outside of the class anyway
//but these will get you going
//future plan would be a pair containing string | point to conversion function
//with a macro to register your conversion function with some container class


};
//class

class ParsedInput {
private:
//checking functions
bool CheckDuplicate(std::vector<ArgumentPair> *a);
bool CheckDuplicate(std::vector<ExpectedInput> *a);
int FindNameIndexArg(const std::string &a);
int FindNameIndexDef(const std::string &a);
int FindNameIndexExp(const std::string &a);
bool CheckExtraArgument();
void CheckNumArgument();
public:
template <typename T>
void WriteOut(const std::string &name, std::vector<T>& myParsedData,const std::string& TYPE);

std::vector<ArgumentPair> *argp;
std::vector<ArgumentPair>* dargp; 
std::vector<ExpectedInput>* eargp;

ParsedInput(){};
ParsedInput(std::vector<ArgumentPair> *argp_, std::vector<ArgumentPair>* dargp_ ,std::vector<ExpectedInput>* eargp_); 
};
template <typename T>//don't define texmplate arguments in cpp file
void ParsedInput::WriteOut(const std::string& name, std::vector<T>& myParsedData, const std::string & TYPE)

{
int expInd = FindNameIndexExp(name);
//std::cout << (*eargp)[expInd].type << std::endl;
if (expInd < 0){std::cout << "Argument " << name << " not expected in WriteOut." << std::endl; throw std::exception();}
if ((*eargp)[expInd].type != TYPE){std::cout << "Type " << TYPE << " not expected in WriteOut, does not match expected arguments." << std::endl; throw std::exception();}
if (!Converters::IsRegistered(TYPE)){std::cout << "Type " << TYPE << " not expected in Converted, make sure macro matches expected type declaration?" << std::endl; throw std::exception();}
//done most I can, belt and suspenders
int argInd = FindNameIndexArg(name);
int defInd = FindNameIndexDef(name); 
if (argInd < 0 && defInd < 0){std::cout << "Argument" << name << " not known by default or cmd arguments. This indicates a bug in the parser, as it should have been caught in previous check. Report to Jake." << std::endl; throw std::exception();}
void * myVoidParsedData = &myParsedData;
fp func = Converters::ReturnConverter(TYPE);
if (argInd >= 0) {(*func)((*argp)[argInd].second , myVoidParsedData); return;}
if (defInd >= 0) {(*func)((*dargp)[defInd].second , myVoidParsedData); return;}

return;
}
//struct ExpectedInputHelpMessage
//{
//std::vector< ExpectedInput> eargpt;
//std::string *prologue;
//std::string *epilogue;
//ExpectedInputHelpMessage(const ExpectedInput *eargpt_[]);
//void regPro(std::string *prologue_);
//void regEpi(std::string *epilogue_);
//void PrintHelp();
//};

#endif
