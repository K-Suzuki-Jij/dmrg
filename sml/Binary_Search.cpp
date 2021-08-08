#include <vector>

template <typename T> long Binary_Search(const std::vector<T> &Array, long min, long max, T target_val) {
   
   if (max < min) {
      return -1;
   }
   else {
      long mid = min + (max - min)/2;
      if (Array[mid] > target_val) {
         return Binary_Search(Array, min, mid - 1, target_val);
      }
      else if (Array[mid] < target_val) {
         return Binary_Search(Array, mid + 1, max, target_val);
      }
      else {
         return mid;
      }
   }
   
}

template long Binary_Search<char>  (const std::vector<char>  &Array, long min, long max, char  target_val);
template long Binary_Search<short> (const std::vector<short> &Array, long min, long max, short target_val);
template long Binary_Search<int>   (const std::vector<int>   &Array, long min, long max, int   target_val);
template long Binary_Search<long>  (const std::vector<long>  &Array, long min, long max, long  target_val);
template long Binary_Search<unsigned char>  (const std::vector<unsigned char>  &Array, long min, long max, unsigned char  target_val);
template long Binary_Search<unsigned short> (const std::vector<unsigned short> &Array, long min, long max, unsigned short target_val);
template long Binary_Search<unsigned int>   (const std::vector<unsigned int>   &Array, long min, long max, unsigned int   target_val);
template long Binary_Search<unsigned long>  (const std::vector<unsigned long>  &Array, long min, long max, unsigned long  target_val);
