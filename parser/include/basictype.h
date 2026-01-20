#ifndef H_BASICTYPE
#define H_BASICTYPE


#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

#include "converters/include/convertermacro.h"

class BasicType {
public:
int a;
std::string b;
double c;
BasicType(){};
BasicType(int a_, std::string b_, double c);
void Print();
};

BasicType::BasicType(int a_, std::string b_, double c_){
a = a_;
b = b_;
c = c_;
}

void BasicType::Print()
{
std::cout << a << " " << b << " " << c << std::endl;
}

namespace BasicConverter {
//static void converter(std::vector<std::string> &s, void *data);

void converter(std::vector<std::string> &s, void *data) {
std::vector<BasicType>* datap = static_cast<std::vector<BasicType>* >(data);
  //in this case, I know to expect data in groups of 3
int a =  std::strtol(s[0].c_str(), NULL, 10);
std::string b = s[1];
double c = std::strtod(s[2].c_str(), NULL);
BasicType T(a,b,c);
(*datap).push_back(T);
return;
  }
}
REGISTER_CLASS(BasicType, "BASICTYPE",&BasicConverter::converter)




#endif 
