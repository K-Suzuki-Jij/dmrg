#include <omp.h>
#include <iostream>
#include "DMRG.hpp"

void DMRG_Transform_Matrix_Two(const CRS &M_On, std::vector<CRS> &M_Out, const DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, std::vector<int> &Ele_LL, std::string Sign_Flag, int LL_site_now, int cf_origin, int p_threads) {
   
   if (cf_origin > LL_site_now) {
      return;
   }
   
   M_Out.resize(LL_site_now + 1);
   
   std::vector<int> Dummy;
   CRS M_Origin;
   CRS M_On_T;
   
   Matrix_Transpose(M_On, M_On_T);
   
   if (cf_origin == 0) {
      M_Origin = M_On_T;
      Matrix_Matrix_Product(M_Origin, M_On, M_Out[0]);
   }
   else {
      DMRG_Basis_LLLR Basis_LLLR;
      Basis_LLLR.LL  = Basis_System.LL_LLLR_Stored[cf_origin - 1];
      Basis_LLLR.LR  = Basis_System.LR_LLLR_Stored[cf_origin - 1];
      Basis_LLLR.Inv = Basis_System.Inv_LLLR_Stored[cf_origin - 1];
      
      DMRG_Renormalize_Matrix("LR", M_On_T, Basis_System.Trans_Mat_Stored[cf_origin - 1], M_Origin, Block, "No", Basis_LLLR, std::vector<int>(), p_threads);
      
      CRS Temp_M;
      Matrix_Matrix_Product(M_On_T, M_On, Temp_M);
      DMRG_Renormalize_Matrix("LR", Temp_M, Basis_System.Trans_Mat_Stored[cf_origin - 1], M_Out[0], Block, "No", Basis_LLLR, std::vector<int>(), p_threads);
   }
      
   for (int LL_site = cf_origin; LL_site < LL_site_now; LL_site++) {
      
      std::vector<CCS> Trans_Mat(p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (int thread_num = 0; thread_num < p_threads; thread_num++) {
         Trans_Mat[thread_num] = Basis_System.Trans_Mat_Stored[LL_site];
      }
      
      std::vector<DMRG_Basis_LLLR> Basis_LLLR(p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (int thread_num = 0; thread_num < p_threads; thread_num++) {
         Basis_LLLR[thread_num].LL  = Basis_System.LL_LLLR_Stored[LL_site];
         Basis_LLLR[thread_num].LR  = Basis_System.LR_LLLR_Stored[LL_site];
         Basis_LLLR[thread_num].Inv = Basis_System.Inv_LLLR_Stored[LL_site];
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int r = 0; r <= LL_site - cf_origin; r++) {
         int thread_num = omp_get_thread_num();
         DMRG_Renormalize_Matrix("LL", M_Out[r], Trans_Mat[thread_num], M_Out[r], Block, "No", Basis_LLLR[thread_num], std::vector<int>(), p_threads);
      }
      
      DMRG_Renormalize_Matrix_LLLR(M_Origin, M_On, Trans_Mat[0], M_Out[LL_site - cf_origin + 1], Ele_LL, Sign_Flag, Basis_LLLR[0]);
      
      DMRG_Renormalize_Matrix("LL", M_Origin, Trans_Mat[0], M_Origin, Block, "No", Basis_LLLR[0], std::vector<int>(), p_threads);

   }
   
   
}
