template <typename T> int Delta_Function(T i, T j) {
   if (i == j) {
      return 1;
   }
   else {
      return 0;
   }
}

template int Delta_Function<char>  (char   i, char   j);
template int Delta_Function<short> (short  i, short  j);
template int Delta_Function<int  > (int    i, int    j);
template int Delta_Function<long > (long   i, long   j);

