#include "SML.hpp"

void Partition_Integer(int div_num, int max_num, std::vector<std::vector<int>> &Out_Array) {
   
   std::vector<int> Temp;
   
   if (Out_Array.size() > 0 && div_num >= max_num) {
      Temp = Out_Array[Out_Array.size() - 1];
   }
   
   if (div_num == 0 && Out_Array.size() == 0) {
      Out_Array.push_back(std::vector<int>());
   }
   else if (div_num == 1) {
      if (Out_Array.size() == 0) {
         Out_Array.push_back(std::vector<int>());
      }
      Out_Array[Out_Array.size() - 1].push_back(1);
   }
   else if (max_num == 1) {
      if (Out_Array.size() == 0) {
         Out_Array.push_back(std::vector<int>());
      }
      for (int i = 0; i < div_num; i++) {
         Out_Array[Out_Array.size() - 1].push_back(1);
      }
   }
   else {
      if (div_num >= max_num) {
         if (Out_Array.size() == 0) {
            Out_Array.push_back(std::vector<int>());
         }
         Out_Array[Out_Array.size() - 1].push_back(max_num);
         Partition_Integer(div_num - max_num, max_num, Out_Array);
         Out_Array.push_back(std::vector<int>());
         for (auto& v: Temp) {
            Out_Array[Out_Array.size() - 1].push_back(v);
         }
      }
      Partition_Integer(div_num, max_num - 1, Out_Array);
   }
   
}
