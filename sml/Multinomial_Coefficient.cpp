#include <vector>
#include <map>
#include "SML.hpp"

long Multinomial_Coefficient(std::vector<int> &Array) {
   
   long array_size = Array.size();
   std::map<int, int> Number;
   
   for (long i = 0; i < array_size; i++) {
      Number[Array[i]] = 0;
   }
   
   for (long i = 0; i < array_size; i++) {
      Number[Array[i]] += 1;
   }
   
   long out  = 1;
   
   for (auto iter = Number.begin(); iter != Number.end(); iter++) {
      out  *= Binomial_Coefficient(array_size, iter->second);
      array_size -= iter->second;
   }
   
   return out;
   
}
