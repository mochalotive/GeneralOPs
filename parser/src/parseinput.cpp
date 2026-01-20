//#include "parser/include/parseinput.h"
#include "parser/include/parseinput.h"


//helper functions for argument pair that don't live in a class
//but I don't want to give global scope
//contains:
//printArgumentPair -> prints all info
//checkDoubleHyphen -> checks if -- term
//parseArgvToPair -> parses user input into usable ArgumentPair data 
bool operator==(ArgumentPair& a, ArgumentPair& b)
{
if (a.first != b.first) return (false);
if (a.second.size() != b.second.size()) return(false);
for (size_t i=0; i<a.second.size(); i++)
  {
    if (a.second[i] != b.second[i]) return(false);
  }
  return(true);
};
//ExpectedInputHelpMessage::ExpectedInputHelpMessage(const ExpectedInput *eargpt_[]){  std::vector<ExpectedInput> eargpt((*eargpt_), (*eargpt_)+sizeof((*eargpt_))/sizeof((*eargpt_[0])));}
//void ExpectedInputHelpMessage::regPro(std::string *prologue_){prologue=prologue_;}
//void ExpectedInputHelpMessage::regEpi(std::string *epilogue_){epilogue=epilogue_;}
//void ExpectedInputHelpMessage::PrintHelp()
//{std::cout << (*prologue) << std::endl;  
//  for (size_t i=0; i<eargpt.size(); i++)
//  {eargpt[i].PrintExpected();}
// std::cout << (*epilogue) << std::endl;}



namespace ParsingHelperNS {
void printArgumentPair(ArgumentPair argp) {
std::cout << "Argument: "<<argp.first << std::endl;
std::cout << "supplied values: " << std::endl;
for (size_t i =0; i < argp.second.size(); i++) {
  std::cout << argp.second[i] << std::endl;
  }
  return;
}
//determines if we have a --term, and does some light error checking
bool checkDoubleHyphen(std::string const & argument)
{
  int i = argument.size();
  switch (i){

    case 0:  {std::cout << "Error, size zero argument in parseinput.h" << std::endl; throw std::exception();}
    case 1:  {return(false);}
    case 2:  {if (argument != "--") return(false);
            else {std::cout << "Error, blank argument -- in parseunput.h" << std::endl; throw std::exception();}
             }
    default: {if (argument[0] == '-' && argument[1] == '-') return(true);
              else return(false);
             }
          }
}
//runs through all the arguments and parses them out into pairs
//returns a vectors of pairs, with each pair being a --term and all its constituent arguments.
std::vector<ArgumentPair> parseArgvToPair(int argcount, const char* arguments_[], bool cmd)
{
  int move = 1;
  if (!cmd)move = 0;//for default arguments
  std::vector<std::string> arguments(arguments_+move, arguments_ + argcount );//convert evil char array to vector of strings
  std::vector<ArgumentPair> returnPairs;//define ArgumentPair we will be returning

  bool hyphenTerm = false; //indicates we need to start a new pair
  //immediately, if first argument isnt --term, no dice. 
  if (!checkDoubleHyphen(arguments[0])){std::cout << "Error, first argument needs to be a --term! from parseinput.h" << std::endl; throw std::exception();}
  if (!cmd) argcount++;
  for (int i=0; i< argcount-1; i++) //filthy integer comparison :vomit:
  {
  //if (hyphenTerm) {returnPairs.push_back(std::pair<std::string, std::vector<std::string> >(std::string(arguments[i].begin()+2), arguments[i].end()), std::vector<std::string>); hyphenTerm=false;}
  hyphenTerm = checkDoubleHyphen(arguments[i]);
  if (hyphenTerm) {returnPairs.push_back(std::pair<std::string, std::vector<std::string> >(std::string(arguments[i].substr(2)), std::vector<std::string>())); hyphenTerm=false;}
  else returnPairs[returnPairs.size()-1].second.push_back(arguments[i]);
  

  }
return returnPairs;
}
}//end namespace

//ExpectedInput
ExpectedInput::ExpectedInput(std::string parameter_, int arguments_, bool mandatory_, std::string type_, std::string errorMessage_, bool strict_)
{
parameter=parameter_;
expectedNumInputs=size_t(arguments_);
if (expectedNumInputs <1) {std::cout << "Must pass at least one additional argument for input --" << parameter_ << ", if you want a boolean flag, do --foo on or --foo off explicitly." << std::endl; throw std::exception();}
mandatory=mandatory_;
type=type_;
errorMessage=errorMessage_;
strict=strict_;
}
void ExpectedInput::PrintExpected()
{
std::cout << "Argument: " << parameter;
if (mandatory && strict) std::cout << " is a mandatory argument. It expects "  << expectedNumInputs << " arguments." << std::endl;
else if (mandatory && !strict) std::cout << " is a mandatory argument. It has a variable number of inputs. It expects a maximum of "  << expectedNumInputs << " arguments." << std::endl;
else if (!mandatory && strict) std::cout << " is an optioanl argument. It expects "  << expectedNumInputs << " arguments." << std::endl;
else std::cout << " is an optional argument. It has a variable number of inputs. It expects a maxiumum of "  << expectedNumInputs << " arguments." << std::endl;

return;

}

//ParsedInput
ParsedInput::ParsedInput(std::vector<ArgumentPair> *argp_, std::vector<ArgumentPair>* dargp_ ,std::vector<ExpectedInput>* eargp_){
argp=argp_;
dargp=dargp_;
eargp=eargp_;
if (CheckDuplicate(argp)) {std::cout << "Duplicate Value in cmd arguments!" << std::endl; throw std::exception();};
if (CheckDuplicate(dargp)) {std::cout << "Duplicate Value in default arguments!" << std::endl; throw std::exception();};
if (CheckDuplicate(eargp)) {std::cout << "Duplicate Value in expected arguments!" << std::endl; throw std::exception();};
;
//now for the checking :)
CheckExtraArgument();//makes sure no mandatory arguments are un-accounted for, either mandatory or default
//now we know everyone is accounted for: first check and fill in defaults
//then apply 
CheckNumArgument();


};

bool ParsedInput::CheckDuplicate(std::vector<ArgumentPair> *a)
{
  for (size_t i=0; i<a->size(); i++){
  for (size_t j=0; j<a->size(); j++){
    if ((*a)[i].first == (*a)[j].first && i!=j)
      return(true);
    }
    }
return(false);
}
bool ParsedInput::CheckDuplicate(std::vector<ExpectedInput> *a)
{
  for (size_t i=0; i<a->size(); i++){
  for (size_t j=0; j<a->size(); j++){

    if ((*a)[i].parameter == (*a)[j].parameter && i!=j)
      return(true);
    }
    }
return(false);
}
int ParsedInput::FindNameIndexArg(const std::string &a)
{
for (size_t i=0; i<argp->size(); i++){if ((*argp)[i].first==a) return(i);}
return(-1);
}
int ParsedInput::FindNameIndexDef(const std::string &a)
{
for (size_t i=0; i<dargp->size(); i++){if ((*dargp)[i].first==a) return(i);}
return(-1);
}
int ParsedInput::FindNameIndexExp(const std::string &a)
{
for (size_t i=0; i<eargp->size(); i++){if ((*eargp)[i].parameter==a) return(i);}
return(-1);
}
bool ParsedInput::CheckExtraArgument()
{
std::vector<std::string> argNames;
std::vector<std::string> defNames;
std::vector<std::string> expNames;
std::vector<bool> manArgs;


//read names into list
for (size_t i=0; i<argp->size(); i++){argNames.push_back((*argp)[i].first);}
for (size_t i=0; i<dargp->size(); i++){defNames.push_back((*dargp)[i].first);}
for (size_t i=0; i<eargp->size(); i++){expNames.push_back((*eargp)[i].parameter);
if ((*eargp)[i].mandatory) manArgs.push_back(true);
else manArgs.push_back(false);
}
//now quickly check for extraneous values:
bool found=false;
for (size_t i=0; i<argNames.size(); i++)
{
  found=false;
  for (size_t j=0; j<expNames.size(); j++)
  {
  if (argNames[i]==expNames[j]) found=true; 
  }
if (!found){std::cout << "Error " << argNames[i] << " not found in expected values list. " <<std::endl; throw std::exception();}
}
for (size_t i=0; i<defNames.size(); i++)
{
  found=false;
  for (size_t j=0; j<expNames.size(); j++)
  {
  if (defNames[i]==expNames[j]) found=true; 
  }
if (!found){std::cout << "Error " << defNames[i] << " not found in default values list. " <<std::endl; throw std::exception();}

}
//so everything is in the expected values, now check if all defaults are covered
//and if all mandatory arguments are present
bool foundArg = true;
bool foundDef = true;
for (size_t k=0; k<expNames.size(); k++) 
{
if (manArgs[k]) 
{
  foundArg = false;
  for (size_t i=0; i<argNames.size(); i++) {if (expNames[k]==argNames[i] ) {foundArg=true;}}
if (!foundArg) {std::cout << "Mandatory Argument " << expNames[k] << " not found in user input" << std::endl; throw std::exception();}
}
else {
  foundDef = false;
  for (size_t i=0; i<defNames.size(); i++) { if (expNames[k]==defNames[i]) {foundDef=true;}}
if (!foundDef) {std::cout << "Default Argument value " << expNames[k] << " not found in default input list" << std::endl; throw std::exception();}
}
}

return(false); 
}

void ParsedInput::CheckNumArgument()
{
std::vector<std::string> argNames;
std::vector<std::string> defNames;
std::vector<std::string> expNames;
for (size_t i=0; i<argp->size(); i++){argNames.push_back((*argp)[i].first);}
for (size_t i=0; i<dargp->size(); i++){defNames.push_back((*dargp)[i].first);}
for (size_t i=0; i<eargp->size(); i++){expNames.push_back((*eargp)[i].parameter);}

//default args first:
for (size_t i=0; i<dargp->size(); i++)
  {  bool isStrict = (*eargp)[FindNameIndexExp((*dargp)[i].first)].strict;
  if (!isStrict &&  (*dargp)[i].second.size() == 0) {std::cout << "Need at least one argument given to " << (*dargp)[i].first << std::endl; throw std::exception();}
  if (!isStrict &&  (*dargp)[i].second.size() > (*eargp)[FindNameIndexExp((*dargp)[i].first)].expectedNumInputs ){std::cout << "Too many arguments given to " << (*dargp)[i].first << std::endl; throw std::exception();}

  //check if sizes are correct
  if (isStrict && ((*dargp)[i].second.size() != (*eargp)[FindNameIndexExp((*dargp)[i].first)].expectedNumInputs)){std::cout << 
   "Error: number of supplied inputs in defaut argument " << (*dargp)[i].first << " does not match expected number." << std::endl; throw std::exception();}
}
//  //now switch types (pain in the ass)
//  switch((*eargp)[FindNameIndexExp((*dargp)[i].first)].type){
//  case ParsingHelperNS::ISTRING:
//  {break;}//nothing actually to do in this case, so just break.
//  case ParsingHelperNS::IINT:
//  {
//  for (int j=0; j<(*dargp)[i].second.size(); j++){if ( std::strtol((*dargp)[i].second[j].c_str(), NULL, 10)==0 && ((*dargp)[i].second[j] != std::string("0")) ) {std::cout << "Error: " << (*dargp)[i].second[j] << " is not an int" << std::endl; throw std::exception();}
//  }
//  break;
//  }
//  case ParsingHelperNS::IREAL:
//  {
//for (int j=0; j<(*dargp)[i].second.size(); j++){if (std::strtod((*dargp)[i].second[j].c_str(), NULL) == 0 && (*dargp)[i].second[j] != std::string("0")) {std::cout << "Error: " << (*dargp)[i].second[j] << " is not a real!" << std::endl; throw std::exception();}
//}
//break;
//}
//  case ParsingHelperNS::IBOOL:
//  {
//  break;//also nothing really to do here
//  }
//  default: {std::cout << "Error: No Type set for " << (*dargp)[i].first << std::endl; throw std::exception(); break;}
//}

   

//now main input arguments

for (size_t i=0; i<argp->size(); i++)
  {bool isStrict = (*eargp)[FindNameIndexExp((*argp)[i].first)].strict;
  if (!isStrict && ((*argp)[i].second.size()==0) ){std::cout << "Need at least one argument given to " << (*argp)[i].first << std::endl; throw std::exception();}
  if (!isStrict && ((*argp)[i].second.size() > (*eargp)[FindNameIndexExp((*argp)[i].first)].expectedNumInputs) ){std::cout << "Too many arguments given to " << (*argp)[i].first << std::endl; throw std::exception();}

  //check if sizes are correct
  if (isStrict && ((*argp)[i].second.size() != (*eargp)[FindNameIndexExp((*argp)[i].first)].expectedNumInputs)){std::cout << 
   "Error: number of supplied inputs in defaut argument " << (*argp)[i].first << " does not match epxected number." << std::endl; throw std::exception();}
}
//  //now switch types (pain in the ass)
//  switch((*eargp)[FindNameIndexExp((*argp)[i].first)].type){
//  case ParsingHelperNS::ISTRING:
//  {break;}//nothing actually to do in this case, so just break.
//  case ParsingHelperNS::IINT:
//  {
//  for (int j=0; j<(*argp)[i].second.size(); j++)
//    { 
//      if ( std::strtol((*argp)[i].second[j].c_str(), NULL, 10)==0 && ((*argp)[i].second[j].c_str()!="0")) {std::cout << "Error: " << (*argp)[i].second[j] << " is not an int" << std::endl; throw std::exception();}
//  }
//  break;
//  }
//  case ParsingHelperNS::IREAL:
//  {
//for (int j=0; j<(*argp)[i].second.size(); j++){if ( std::strtod((*argp)[i].second[j].c_str(), NULL)==0 && ((*argp)[i].second[j].c_str() !="0")) {std::cout << "Error: " << (*argp)[i].second[j] << " is not a real!" << std::endl; throw std::exception();}
//}
//break;
//}
//  case ParsingHelperNS::IBOOL:
//  {
//  break;//also nothing really to do here
//  }
//  default: {std::cout << "Error: No Type set for " << (*argp)[i].first << std::endl; throw std::exception(); break;}
//}

   
//  }
return;
};


//{ 
//switch(t){
//
//  case ParsingHelperNS::ISTRING:
//{
//  if (FindNameIndexArg(name) > -1)
//  {
//  int i = FindNameIndexArg(name);
//  for (int j=0; j < (*argp)[i].second.size(); j++)
//    {
//    a.push_back((*argp)[i].second[j]);
//    } 
//    return;
//
//  }
//  else {
//    int i = FindNameIndexDef(name);
//    for (int j=0; j < (*dargp)[i].second.size(); j++)
//      {
//        a.push_back((*dargp)[i].second[j]);
//      } 
//      return;
//     }
//}//easy case
//  case ParsingHelperNS::IINT:{}
//  case ParsingHelperNS::IREAL:{}
//
////bool
//  case ParsingHelperNS::IBOOL:
//{
//  if (FindNameIndexArg(name) > -1){
//  int i = FindNameIndexArg(name);
//  if ((*argp)[i].second.size() ==0) {a.push_back(true); return;}
//  for (int j=0; j < (*argp)[i].second.size(); j++)
//    {
//        if ((*argp)[i].second[j] == (*eargp)[FindNameIndexExp(name)].boolString) {a.push_back(true);}
//        else a.push_back(false);
//
//    }
//  return;
//  }
//  else 
//  {
//  int i = FindNameIndexDef(name);
//  if ((*dargp)[i].second.size() ==0) {a.push_back(true); return;}
//  for (int j=0; j < (*dargp)[i].second.size(); j++) 
//    {   
//        if ((*dargp)[i].second[j] == (*eargp)[FindNameIndexExp(name)].boolString) {a.push_back(true);}
//        else a.push_back(false);
//
//    } 
//  return;
//  }
//  }//easy case
////bool
//  default: break;
//}
//return;
//}

//printArgumentPair
//checkDoubleHyphen
//parseArgvToPair
//InputType
//for testing, compile by itself with g++


