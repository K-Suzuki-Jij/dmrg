#include "DMRG.hpp"

void DMRG_Make_Elem_Ham_LLLRRRRL(const DMRG_Onsite_Basis &Basis_Onsite, DMRG_A_Basis_Set &A_Basis, const DMRG_Basis &Basis, const DMRG_Block_Hamiltonian &Ham, std::string Sign) {
   
   DMRG_Make_Elem_Onsite_LLLRRRRL("LLLR", Basis_Onsite.LLLR, Basis_Onsite, A_Basis, Basis.LLLR.LL, Basis.LLLR.LR, Basis.LLLRRRRL.Inv, Ham.LLLR, 1.0, Ham.Ele_RR, Ham.Ele_On, "No");
   DMRG_Make_Elem_Onsite_LLLRRRRL("LLRR", Basis_Onsite.LLRR, Basis_Onsite, A_Basis, Basis.LLRR.LL, Basis.LLRR.RR, Basis.LLLRRRRL.Inv, Ham.LLRR, 1.0, Ham.Ele_RR, Ham.Ele_On, Sign);
   DMRG_Make_Elem_Onsite_LLLRRRRL("LLRL", Basis_Onsite.LLRL, Basis_Onsite, A_Basis, Basis.LLRL.LL, Basis.LLRL.RL, Basis.LLLRRRRL.Inv, Ham.LLRL, 1.0, Ham.Ele_RR, Ham.Ele_On, Sign);
   DMRG_Make_Elem_Onsite_LLLRRRRL("LRRR", Basis_Onsite.LRRR, Basis_Onsite, A_Basis, Basis.LRRR.LR, Basis.LRRR.RR, Basis.LLLRRRRL.Inv, Ham.LRRR, 1.0, Ham.Ele_RR, Ham.Ele_On, "No");
   DMRG_Make_Elem_Onsite_LLLRRRRL("LRRL", Basis_Onsite.LRRL, Basis_Onsite, A_Basis, Basis.LRRL.LR, Basis.LRRL.RL, Basis.LLLRRRRL.Inv, Ham.LRRL, 1.0, Ham.Ele_RR, Ham.Ele_On, Sign);
   DMRG_Make_Elem_Onsite_LLLRRRRL("RRRL", Basis_Onsite.RRRL, Basis_Onsite, A_Basis, Basis.RRRL.RR, Basis.RRRL.RL, Basis.LLLRRRRL.Inv, Ham.RRRL, 1.0, Ham.Ele_RR, Ham.Ele_On, "No");
   
   DMRG_Make_Elem_Zero_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv);
   DMRG_Clear_Check_Basis(A_Basis, Basis.LLLRRRRL.Inv);
   
}
