#ifndef VARIABLE_H_
#define VARIABLE_H_

#include <string>

/*----------------------------------------------------------------------------*/

class Variable {
private:
   std::string valsorig;
   std::string vals;
   float valf;
   int vali;
   bool valb;

protected:
public:
   Variable(const std::string &value);

   void set(const std::string &value);
   float getFloat() const;
   int getInt() const;
   bool getBool() const;
   const std::string &str() const;
   const char* cStr() const;
   void reset();

};

#endif /* VARIABLE_H_ */
