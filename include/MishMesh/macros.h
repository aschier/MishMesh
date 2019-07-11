#pragma once

#define FOR_CVE(e_it, vh) for(auto e_it = mesh.cve_ccwbegin(vh); e_it != mesh.cve_ccwend(vh); e_it++)
#define FOR_CVF(f_it, vh) for(auto f_it = mesh.cvf_ccwbegin(vh); f_it != mesh.cvf_ccwend(vh); f_it++)
#define FOR_CVV(v_it, vh) for(auto v_it = mesh.cvv_ccwbegin(vh); v_it != mesh.cvv_ccwend(vh); v_it++)
#define FOR_CVOH(h_it, vh) for(auto h_it = mesh.cvoh_ccwbegin(vh); h_it != mesh.cvoh_ccwend(vh); h_it++)
#define FOR_CVIH(h_it, vh) for(auto h_it = mesh.cvih_ccwbegin(vh); h_it != mesh.cvih_ccwend(vh); h_it++)

#define FOR_VE(e_it, vh) for(auto e_it = mesh.ve_ccwbegin(vh); e_it != mesh.ve_ccwend(vh); e_it++)
#define FOR_VF(f_it, vh) for(auto f_it = mesh.vf_ccwbegin(vh); f_it != mesh.vf_ccwend(vh); f_it++)
#define FOR_VV(v_it, vh) for(auto v_it = mesh.vv_ccwbegin(vh); v_it != mesh.vv_ccwend(vh); v_it++)
#define FOR_VOH(h_it, vh) for(auto h_it = mesh.voh_ccwbegin(vh); h_it != mesh.voh_ccwend(vh); h_it++)
#define FOR_VIH(h_it, vh) for(auto h_it = mesh.vih_ccwbegin(vh); h_it != mesh.vih_ccwend(vh); h_it++)

#define FOR_CFV(v_it, fh) for(auto v_it = mesh.cfv_ccwbegin(fh); v_it != mesh.cfv_ccwend(fh); v_it++)
#define FOR_CFE(e_it, fh) for(auto e_it = mesh.cfe_ccwbegin(fh); e_it != mesh.cfe_ccwend(fh); e_it++)
#define FOR_CFH(h_it, fh) for(auto h_it = mesh.cfh_ccwbegin(fh); h_it != mesh.cfh_ccwend(fh); h_it++)
#define FOR_CFF(f_it, fh) for(auto f_it = mesh.cff_ccwbegin(fh); f_it != mesh.cff_ccwend(fh); f_it++)

#define FOR_FV(v_it, fh) for(auto v_it = mesh.fv_ccwbegin(fh); v_it != mesh.fv_ccwend(fh); v_it++)
#define FOR_FE(e_it, fh) for(auto e_it = mesh.fe_ccwbegin(fh); e_it != mesh.fe_ccwend(fh); e_it++)
#define FOR_FH(h_it, fh) for(auto h_it = mesh.fh_ccwbegin(fh); h_it != mesh.fh_ccwend(fh); h_it++)
#define FOR_FF(f_it, fh) for(auto f_it = mesh.ff_ccwbegin(fh); f_it != mesh.ff_ccwend(fh); f_it++)
