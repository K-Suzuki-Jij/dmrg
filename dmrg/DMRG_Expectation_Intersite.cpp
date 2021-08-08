#include "DMRG.hpp"

void DMRG_Expectation_Intersite(std::vector<double> &Out, int site_ref, std::vector<CRS> &M_CF_LL, std::vector<CRS> &M_LL, CRS &M_On, std::vector<CRS> &M_RR, std::vector<double> &Vec, DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads) {
   
   std::vector<double> Vec_Ref, Vec_R;
   
   if (site_ref > LL_site + 1) {
      Out.resize(0);
      return;
   }
   
   Out.resize(LL_site + RR_site + 4 - site_ref + 1);
   
   if (site_ref == LL_site + 1) {
      DMRG_Matrix_Vector_Product("LR", M_On, Vec, Vec_Ref, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   }
   else {
      DMRG_Matrix_Vector_Product("LL", M_LL[site_ref], Vec, Vec_Ref, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   }
   
   //O*O to O*LL_site
   for (int site = site_ref; site <= LL_site; site++) {
      int r = site - site_ref;
      DMRG_Matrix_Vector_Product("LL", M_CF_LL[r], Vec, Vec_R, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
      Out[r] = Inner_Product(Vec, Vec_R, p_threads);
   }
   
   //O*LR_site
   DMRG_Matrix_Vector_Product("LR", M_On, Vec, Vec_R, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   Out[LL_site + 1 - site_ref] = Inner_Product(Vec_Ref, Vec_R, p_threads);
   
   //O*RL_site
   DMRG_Matrix_Vector_Product("RL", M_On, Vec, Vec_R, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
   Out[LL_site + 2 - site_ref] = Inner_Product(Vec_Ref, Vec_R, p_threads);
   
   //O*RR_site
   for (int site = RR_site; site >= 0; site--) {
      int r = RR_site + LL_site + 3 - site - site_ref;
      DMRG_Matrix_Vector_Product("RR", M_RR[site], Vec, Vec_R, Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
      Out[r] = Inner_Product(Vec_Ref, Vec_R, p_threads);
   }
   
}

void DMRG_Expectation_Intersite(std::vector<double> &Out, int site_ref, std::vector<CRS> &M_CF_LL, std::vector<CRS> &M_LL, CRS &M_On, std::vector<CRS> &M_RR, std::vector<double> &Vec, DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, std::vector<DMRG_Basis_LLLRRRRL> &W_Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads) {
   
   int sector_size = (int)W_Basis_LLLRRRRL.size();
   
   std::vector<std::vector<double>> Vec_Ref(sector_size), Vec_R(sector_size);
   
   if (site_ref > LL_site + 1) {
      Out.resize(0);
      return;
   }
   
   Out.resize(LL_site + RR_site + 4 - site_ref + 1);
   
   if (site_ref == LL_site + 1) {
      for (int sector = 0; sector < sector_size; sector++) {
         DMRG_Matrix_Vector_Product("LR", M_On, Vec, Vec_Ref[sector], Basis_LLLRRRRL.Inv, W_Basis_LLLRRRRL[sector], p_threads);
      }
   }
   else {
      for (int sector = 0; sector < sector_size; sector++) {
         DMRG_Matrix_Vector_Product("LL", M_LL[site_ref], Vec, Vec_Ref[sector], Basis_LLLRRRRL.Inv, W_Basis_LLLRRRRL[sector], p_threads);
      }
   }
   
   //O*O to O*LL_site
   for (int site = site_ref; site <= LL_site; site++) {
      int r = site - site_ref;
      DMRG_Matrix_Vector_Product("LL", M_CF_LL[r], Vec, Vec_R[0], Basis_LLLRRRRL.Inv, Basis_LLLRRRRL, p_threads);
      Out[r] = Inner_Product(Vec, Vec_R[0], p_threads);
   }
   
   //O*LR_site
   Out[LL_site + 1 - site_ref] = 0.0;
   for (int sector = 0; sector < sector_size; sector++) {
      DMRG_Matrix_Vector_Product("LR", M_On, Vec, Vec_R[sector], Basis_LLLRRRRL.Inv, W_Basis_LLLRRRRL[sector], p_threads);
      Out[LL_site + 1 - site_ref] += Inner_Product(Vec_Ref[sector], Vec_R[sector], p_threads);
   }
   
   //O*RL_site
   Out[LL_site + 2 - site_ref] = 0.0;
   for (int sector = 0; sector < sector_size; sector++) {
      DMRG_Matrix_Vector_Product("RL", M_On, Vec, Vec_R[sector], Basis_LLLRRRRL.Inv, W_Basis_LLLRRRRL[sector], p_threads);
      Out[LL_site + 2 - site_ref] += Inner_Product(Vec_Ref[sector], Vec_R[sector], p_threads);
   }
   
   //O*RR_site
   for (int site = RR_site; site >= 0; site--) {
      int r = RR_site + LL_site + 3 - site - site_ref;
      Out[r] = 0.0;
      for (int sector = 0; sector < sector_size; sector++) {
         DMRG_Matrix_Vector_Product("RR", M_RR[site], Vec, Vec_R[sector], Basis_LLLRRRRL.Inv, W_Basis_LLLRRRRL[sector], p_threads);
         Out[r] += Inner_Product(Vec_Ref[sector], Vec_R[sector], p_threads);
      }
   }
   
}
