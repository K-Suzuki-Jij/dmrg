#include "DMRG.hpp"

void DMRG_Expectation_Onsite(std::vector<double> &Out, const std::vector<CRS> &M_LL, const CRS &M_On, const std::vector<CRS> &M_RR, const std::vector<double> &Vec, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads) {
   
   std::vector<double> Temp_V;
   
   Out.resize(LL_site + RR_site + 4);
   
   //LL_site
   for (int site = 0; site <= LL_site; site++) {
      DMRG_Matrix_Vector_Product("LL", M_LL[site], Vec, Temp_V, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
      Out[site] = Inner_Product(Vec, Temp_V, p_threads);
   }
   
   //LR_site
   DMRG_Matrix_Vector_Product("LR", M_On, Vec, Temp_V, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   Out[LL_site + 1] = Inner_Product(Vec, Temp_V, p_threads);
   
   //RL_site
   DMRG_Matrix_Vector_Product("RL", M_On, Vec, Temp_V, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   Out[LL_site + 2] = Inner_Product(Vec, Temp_V, p_threads);
   
   //RR_site
   for (int site = RR_site; site >= 0; site--) {
      DMRG_Matrix_Vector_Product("RR", M_RR[site], Vec, Temp_V, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
      Out[RR_site + LL_site + 3 - site] = Inner_Product(Vec, Temp_V, p_threads);
   }
   
}
