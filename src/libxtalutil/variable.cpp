#include <cstdlib>
#include "variable.h"

using namespace std;

/*----------------------------------------------------------------------------*/

Variable::Variable(const string &value) {
   this->valsorig = value;
   set(value);
}

/*----------------------------------------------------------------------------*/

void Variable::set(const string &value) {
   vals = value;
   valf = atof(value.c_str());
   vali = atoi(value.c_str());
   valb = (value == "true") || (value == "yes") || (value == "1");
}

/*----------------------------------------------------------------------------*/

const string& Variable::str() const {
   return vals;
}

/*----------------------------------------------------------------------------*/

const char* Variable::cStr() const {
   return vals.c_str();
}

/*----------------------------------------------------------------------------*/

float Variable::getFloat() const {
   return valf;
}

/*----------------------------------------------------------------------------*/

int Variable::getInt() const {
   return vali;
}

/*----------------------------------------------------------------------------*/

bool Variable::getBool() const {
   return valb;
}

/*----------------------------------------------------------------------------*/

void Variable::reset() {
   set(valsorig);
}
