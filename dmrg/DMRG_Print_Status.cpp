#include <iostream>
#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "DMRG.hpp"

void DMRG_Print_Status(int num_of_qn, std::string BC, const DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block,const  DMRG_Param &Dmrg_Param, const DMRG_Time &Time) {
   
   int LL_site    = Block.LL_site;
   int RR_site    = Block.RR_site;
   int dim_LL     = Block.dim_LL;
   int dim_RR     = Block.dim_RR;
   int dim_onsite = Block.dim_onsite;

   if (BC == "OBC") {
      
      std::cout << "["        << BC << "]";
      std::cout << "E/N="     << std::fixed      << std::setprecision(15) << GS.val/(LL_site + RR_site + 4);
      std::cout << "("        << std::scientific << std::setprecision(1)  << GS.error << ",";
      std::cout << std::scientific << std::setprecision(1)    << GS.tr_error << "),";
      std::cout << "N="       << std::left    << std::setw(3) << LL_site + RR_site + 4;
      std::cout << "("        << std::left    << std::setw(3) << LL_site + 1;
      std::cout << ",1,1,"    << std::right   << std::setw(3) << RR_site + 1 << "),";
      std::cout << "dim="     << std::left    << std::setw(8) << (int)Basis_LLLRRRRL.LL.size();
      std::cout << "("        << std::left    << std::setw(4) << dim_LL << ",";
      std::cout << std::left  << std::setw(2) << dim_onsite << ",";
      std::cout << std::right << std::setw(2) << dim_onsite << ",";
      std::cout << std::right << std::setw(4) << dim_RR << "),";
      
      int c_ele = 0, c_p = 0, c_sz = 0;
      for (int i = 0; i < num_of_qn; i++) {
         if (Basis_System.qn_LL_LL_stored_ele_start <= i && i < Basis_System.qn_LL_LL_stored_ele_end) {
            std::cout << "Nc" << c_ele << "=" << std::left << std::setw(3) << Basis_LLLRRRRL.QN[i] << ",";c_ele++;
         }
         else if (Basis_System.qn_LL_LL_stored_parity_start <= i && i < Basis_System.qn_LL_LL_stored_parity_end) {
            std::cout << "P" << c_p << "=" << Basis_LLLRRRRL.QN[i] << ",";c_p++;
         }
         else {
            std::cout << "Sz" << c_sz << "=" << std::left << std::setw(5) << std::fixed << std::setprecision(1) << Basis_LLLRRRRL.QN[i]/2.0 << ",";c_sz++;
         }
      }
      
      std::cout << "T:"       << std::left    << std::setw(5) << std::setprecision(1) << Time.total;
      std::cout << "(D:"      << std::left    << std::setw(5) << std::setprecision(1) << Time.diag << ",";
      std::cout << "II:"      << std::left    << std::setw(5) << std::setprecision(1) << Time.inv_iter << ",";
      std::cout << "Ham:"     << std::left    << std::setw(5) << std::setprecision(1) << Time.make_ham << ",";
      std::cout << "Oth:"     << std::left    << std::setw(4) << std::setprecision(1) << Time.total - (Time.diag + Time.inv_iter + Time.make_ham) << ")";
      std::cout << ","        << std::right   << std::setw(3) << Dmrg_Param.renorm_now_iter << "/" << Dmrg_Param.renorm_tot_iter;
      std::cout << ","        << std::right   << std::setw(3) << Dmrg_Param.param_now_iter + 1 << "/" << Dmrg_Param.param_tot_iter;
      std::cout << std::endl;
      
      mkdir("./result", 0777);
      std::string Out_Name = "./result/log.txt";
      std::ofstream file(Out_Name, std::ios::app);
      
      file << "["        << BC << "]";
      file << "E/N="     << std::fixed      << std::setprecision(15) << GS.val/(LL_site + RR_site + 4);
      file << "("        << std::scientific << std::setprecision(1)  << GS.error << ",";
      file << std::scientific << std::setprecision(1)    << GS.tr_error << "),";
      file << "N="       << std::left    << std::setw(3) << LL_site + RR_site + 4;
      file << "("        << std::left    << std::setw(3) << LL_site + 1;
      file << ",1,1,"    << std::right   << std::setw(3) << RR_site + 1 << "),";
      file << "dim="     << std::left    << std::setw(8) << (int)Basis_LLLRRRRL.LL.size();
      file << "("        << std::left    << std::setw(4) << dim_LL << ",";
      file << std::left  << std::setw(2) << dim_onsite << ",";
      file << std::right << std::setw(2) << dim_onsite << ",";
      file << std::right << std::setw(4) << dim_RR << "),";
    
      c_ele = 0; c_p = 0; c_sz = 0;
      for (int i = 0; i < num_of_qn; i++) {
         if (Basis_System.qn_LL_LL_stored_ele_start <= i && i < Basis_System.qn_LL_LL_stored_ele_end) {
            file << "Nc" << c_ele << "=" << std::left << std::setw(3) << Basis_LLLRRRRL.QN[i] << ",";c_ele++;
         }
         else if (Basis_System.qn_LL_LL_stored_parity_start <= i && i < Basis_System.qn_LL_LL_stored_parity_end) {
            file << "P=" << c_p << Basis_LLLRRRRL.QN[i] << ",";c_p++;
         }
         else {
            file << "Sz=" << c_sz << std::left  << std::setw(5) << Basis_LLLRRRRL.QN[i]/2.0 << ",";c_sz++;
         }
      }
      
      file << "T:"       << std::left    << std::setw(5) << std::setprecision(1) << Time.total;
      file << "(D:"      << std::left    << std::setw(5) << std::setprecision(1) << Time.diag << ",";
      file << "II:"      << std::left    << std::setw(5) << std::setprecision(1) << Time.inv_iter << ",";
      file << "Ham:"     << std::left    << std::setw(5) << std::setprecision(1) << Time.make_ham << ",";
      file << "Oth:"     << std::left    << std::setw(4) << std::setprecision(1) << Time.total - (Time.diag + Time.inv_iter + Time.make_ham) << ")";
      file << ","        << std::right   << std::setw(3) << Dmrg_Param.renorm_now_iter << "/" << Dmrg_Param.renorm_tot_iter;
      file << ","        << std::right   << std::setw(3) << Dmrg_Param.param_now_iter + 1 << "/" << Dmrg_Param.param_tot_iter;
      file << std::endl;
      
      file.close();
      
   }
   
}
