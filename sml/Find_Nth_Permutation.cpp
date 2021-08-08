#include <iostream>
#include <vector>
#include <map>
#include "SML.hpp"

void Find_Nth_Permutation(std::vector<int> &Array, long target_number) {
   
   if (Multinomial_Coefficient(Array) <= target_number) {
      std::cout << "Error in Find_Nth_Permutation" << std::endl;
      std::cout << "target_number=" << target_number << std::endl;
      std::exit(0);
   }
   
   long array_size = Array.size();
   std::map<int, int> Number;
   
   for (long i = 0; i < array_size; i++) {
      Number[Array[i]] = 0;
   }
   
   for (long i = 0; i < array_size; i++) {
      Number[Array[i]] += 1;
   }
      
   for (long i = 0; i < array_size; i++) {
      
      long temp1 = 0;
      for (auto iter1 = Number.begin(); iter1 != Number.end(); iter1++) {
         
         if (iter1->second > 0) {
            iter1->second -= 1;
            
            long temp2 = 1;
            long size  = array_size - (i + 1);
            for (auto iter2 = Number.begin(); iter2 != Number.end(); iter2++) {
               temp2 *= Binomial_Coefficient(size, iter2->second);
               size  -= iter2->second;
            }
            
            temp1 += temp2;
            if (temp1 > target_number) {
               Array[i] = iter1->first;
               temp1 -= temp2;
               break;
            }
            
            iter1->second += 1;
         }
         
      }
      target_number -= temp1;
   }
   
   if (target_number != 0) {
      std::cout << "Error in Find_Nth_Permutation" << std::endl;
      std::cout << "Can't find corresponding permutation, " << target_number <<  std::endl;
      std::exit(0);
   }
      
}
