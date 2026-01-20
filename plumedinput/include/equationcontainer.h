#ifndef H_EQUATIONCONTAINER
#define H_EQUATIONCONTAINER

#include <vector>
#include <string>
#include <iostream>

class EquationContainer {
private:
std::vector<std::string> equations;
public:
EquationContainer();
EquationContainer(std::vector<std::string> &s);
EquationContainer(std::string &s);
size_t getSize();
std::string getString(size_t i);
void addString(std::string s);
void addString(std::vector<std::string> t);
//friend std::ostream &operator<<(std::ostream &str, EquationContainer &container);

};
inline EquationContainer::EquationContainer(){}
inline EquationContainer::EquationContainer(std::vector<std::string> &s)
{
for (size_t i=0; i< s.size(); ++i) equations.push_back(s[i]);
}
inline EquationContainer::EquationContainer(std::string &s)
{
equations.push_back(s);
}
//inline std::ostream operator<<( std::ostream &str, EquationContainer &container)
//{
//for (size_t i=0; i < container.getSize(); ++i)
//  str << container.getString(i) << std::endl;
//}

inline void EquationContainer::addString(std::string s) {equations.push_back(s);}
inline void EquationContainer::addString(std::vector<std::string> t) 
  {
  for (size_t i=0; i< t.size(); i++) addString(t[i]);
  }

inline std::string EquationContainer::getString(size_t i) {return(equations[i]);}
inline size_t EquationContainer::getSize() {return(equations.size());}



#endif
