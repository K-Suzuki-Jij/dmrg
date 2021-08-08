long Binomial_Coefficient(long n, long k) {
   
   long r = 1;
   
   for (long d = 1; d <= k; d++) {
       r *= n--;
       r /= d;
   }
   
   return r;
   
}

