#include <omp.h>
#include "DMRG.hpp"

void DMRG_Transform_Matrix_One(const CRS &M_On, std::vector<CRS> &M_Out, const DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, int LL_site_now, int p_threads) {
   
   M_Out.resize(LL_site_now + 1);
   
   M_Out[0] = M_On;
   
   for (int LL_site = 0; LL_site <= LL_site_now - 1; LL_site++) {
      
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
      for (int site = 0; site <= LL_site; site++) {
         int thread_num = omp_get_thread_num();
         DMRG_Renormalize_Matrix("LL", M_Out[site], Trans_Mat[thread_num], M_Out[site], Block, "No", Basis_LLLR[thread_num], std::vector<int>(), 1);
      }
      
      DMRG_Renormalize_Matrix("LR", M_On, Trans_Mat[0], M_Out[LL_site + 1], Block, "No", Basis_LLLR[0], std::vector<int>(), p_threads);
      
   }
   
}
