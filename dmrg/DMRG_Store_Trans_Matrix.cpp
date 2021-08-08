#include "DMRG.hpp"

void DMRG_Store_Trans_Matrix(CCS &Trans_Mat, DMRG_Basis_Stored &Basis_System, DMRG_Basis_LLLR &Basis_LLLR, DMRG_Param &Dmrg_Param, int LL_site, int system_size) {
   
   
   int c1 = (Dmrg_Param.Enviro_Copy == "Yes") && (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c2 = (Dmrg_Param.Enviro_Copy == "No" ) && (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 3);
   int c3 = (LL_site >= system_size/2 - 2);
   
   if ((c1 || c2 || c3) && Dmrg_Param.Initial_Guess != "Yes") {
      return;
   }
      
   Basis_System.Trans_Mat_Stored[LL_site] = Trans_Mat;
   Basis_System.LL_LLLR_Stored[LL_site]   = Basis_LLLR.LL;
   Basis_System.LR_LLLR_Stored[LL_site]   = Basis_LLLR.LR;
   Basis_System.Inv_LLLR_Stored[LL_site]  = Basis_LLLR.Inv;
   
}
