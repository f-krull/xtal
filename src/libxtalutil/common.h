#ifndef COMMON_H_
#define COMMON_H_

#include <string>
#include <vector>
#include <fstream>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

namespace common {

static const float E = 2.71828182845f; /* eulers number */

static const float PI = 3.14159265358979323846264338327f; /* 180 degrees */
static const float PI_HALF = 1.570796326794896558f; /*  90 degrees */
static const float PI_THIRD = 1.047197551f;
static const float PI_QUATER = 0.785398163397448279f; /*  45 degrees */
static const float PI_DOUBLE = 6.283185307179586232f;
static const float PI_180TH = 0.01745329f; /* deg to rad: deg * PI_180TH */

/* 10 ms is smallest guaranteed wait */
void wait(unsigned int ms = 10);

bool exists(std::string filename);

std::string removePath(const std::string &filename);
std::string removeExtension(const std::string & filename);
std::string s_printf(const char *str, ...);
std::vector<std::string> readList(const char *filename, char comment = '#');

template<typename T>
bool readTableColumns(const std::string &filename, uint32_t nCols,
      std::vector<std::vector<T> > &table) {
   std::ifstream infile;

   if (nCols == 0) {
      return false;
   }
   table.clear();

   infile.open(filename.c_str(), std::ifstream::in);
   std::vector<T> row(nCols);
   while (infile.good()) {
      for (uint32_t i = 0; i < row.size(); ++i) {
         infile >> row[i];
      }
      if (infile.good() == false) {
         break;
      }
      table.push_back(row);
   }
   infile.close();
   return table[0].size() > 0 && table[0].size() == nCols;
}

}

//static void printSizeOf() {
//   std::cout << "size of char: " << sizeof(char) << std::endl;
//   std::cout << "size of void*: " << sizeof(void*) << std::endl;
//   std::cout << "size of unsigned int: " << sizeof(unsigned int) << std::endl;
//   std::cout << "size of float: " << sizeof(float) << std::endl;
//   std::cout << "size of double: " << sizeof(double) << std::endl;
//   std::cout << "size of long: " << sizeof(long) << std::endl;
//   std::cout << "size of unsigned long: " << sizeof(unsigned long) << std::endl;
//}}

#endif /*COMMON_H_*/
