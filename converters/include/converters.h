#ifndef H_CONVERTERS
#define H_CONVERTERS

#include <string>//needed in arguments
#include <vector>//needed for members
#include <cstdlib>//for strtod and strtol
#include <cstring>//for strcmp
#include <iostream>//needed for error i/o
//#include "convertermacro.h"
//data, container, data name



typedef void (*fp)(std::vector<std::string> &, void *);//pointers to functions which return void and take a string and a void pointer and another string.
namespace Converters 
{
  class MySingleton {
  public:
    static MySingleton& getInstance()
    {
      static MySingleton instance;
      return instance;
    }
    std::vector<std::pair<std::string ,fp> >& getConverters()
    {
    return converters;
    }

    private:
    MySingleton() : converters() {} // Constructor initializes myVariable
    
    std::vector<std::pair<std::string ,fp> >converters;//stores converter functions
  };
  inline void PrintAll()//prints all members (which should be initialized before main!)
  {
   MySingleton& singletonInstance = MySingleton::getInstance();
    // Access the vector instance
   std::vector<std::pair<std::string, fp> >& converters = singletonInstance.getConverters(); 
    for (size_t i=0; i<converters.size(); i++)
    {
      std::cout << converters[i].first << std::endl;
    }
    return;
  } 
  inline void AddConverter(std::string name, fp c)
  {
   MySingleton& singletonInstance = MySingleton::getInstance();
    // Access the vector instance
   std::vector<std::pair<std::string, fp> >& converters = singletonInstance.getConverters(); 
    if (converters.size() > 0 ) {
    for (size_t i=0; i<converters.size(); i++)
    { if (strcmp(converters[i].first.c_str(), name.c_str()) == 0) {Converters::PrintAll(); std::cout << "Converter " << name << " already represented in converters.h!" << std::endl; throw std::exception();}  }
    }
    std::pair<std::string, fp> temp;
    temp.first = name;
    temp.second = c;
    converters.push_back(temp); return;
  }//no duplicate checking at the moment, handle that later

  inline fp ReturnConverter(const std::string &name)//returns pointer to converter
  { 
   MySingleton& singletonInstance = MySingleton::getInstance();
    // Access the vector instance
   std::vector<std::pair<std::string, fp> >& converters = singletonInstance.getConverters(); 

    for (size_t i=0; i<converters.size(); i++) 
    {
      if (converters[i].first == name) return(converters[i].second);
    }
    throw std::exception();
    return(fp());//shouldn't ever get here
  }
 
  inline bool IsRegistered(const std::string& name)
  {
   MySingleton& singletonInstance = MySingleton::getInstance();
    // Access the vector instance
   std::vector<std::pair<std::string, fp> >& converters = singletonInstance.getConverters(); 

    for (size_t i=0; i<converters.size(); i++) 
    {
      if (name == converters[i].first) return(true);
        }
    return(false);
  }

}
    //namespace {
    //class initconv { 
    //public: 
    //    initconv() { 
    //        Converters::MySingleton& singletonInstance = Converters::MySingleton::getInstance();
    //        std::vector<std::pair<std::string, fp> >& converters = singletonInstance.getConverters();
    //    } 
    //}; 
    //static initconv initconvInstance;
    //}  //quickly register
    


// Macro to simplify registration. Defies a class called classNameRegistar (className really doesn't matter too much, expecially because its an anonymous namespace)
// then defines a static instance of this class, and within the constructor we register w/ Converters. 
//define int real string and double converters 
#define REGISTER_CLASS(className,classString, fps) \
  namespace Converters {\
  void AddConverter(std::string, fp);}\
  namespace{ \
    class className##Registrar { \
    public: \
        className##Registrar() { \
            Converters::AddConverter(classString , fps); \
        } \
    }; \
    static className##Registrar className##RegistrarInstance;\
    }  

#endif












