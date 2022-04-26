#include "measure/measure_methods.h"
#include "model/model_base.h"
#include "lattice/lattice_base.h"
#include "dqmc_walker.h"


namespace Measure {

    // some aliases
    using RealScalar = double; 
    using GreensFunc = Eigen::MatrixXd;


    // -----------------------------  Method routines for equal-time measurements  -------------------------------

    void Methods::measure_equaltime_config_sign( ScalarObs& equaltime_sign, 
                                                 const DqmcWalker& walker,
                                                 const ModelBase& model,
                                                 const LatticeBase& lattice )
    {
        equaltime_sign.tmp_value() += walker.vecConfigSign().sum();
        equaltime_sign.counts() += walker.TimeSize();
    }


    void Methods::measure_filling_number( ScalarObs& filling_number, 
                                          const DqmcWalker& walker,
                                          const ModelBase& model,
                                          const LatticeBase& lattice )
    {
        // loop over equivalent time slices
        for (auto t = 0; t < walker.TimeSize(); ++t) {
            filling_number.tmp_value() += walker.ConfigSign(t) * 
                ( 2 - (walker.GreenttUp(t).trace()+walker.GreenttDn(t).trace())/lattice.SpaceSize() );
            ++filling_number;
        }
    }


    void Methods::measure_double_occupancy( ScalarObs& double_occupancy,
                                            const DqmcWalker& walker,
                                            const ModelBase& model,
                                            const LatticeBase& lattice )
    {
        for (auto t = 0; t < walker.TimeSize(); ++t) {
            const GreensFunc& gu = walker.GreenttUp(t);
            const GreensFunc& gd = walker.GreenttDn(t);

            RealScalar tmp_double_occu = 0.0;
            for (int i = 0; i < lattice.SpaceSize(); ++i) {
                tmp_double_occu += walker.ConfigSign(t) * (1 - gu(i, i)) * (1 - gd(i, i));
            }
            double_occupancy.tmp_value() += tmp_double_occu/lattice.SpaceSize();
            ++double_occupancy;
        }
    }


    void Methods::measure_kinetic_energy( ScalarObs& kinetic_energy, 
                                 const DqmcWalker& walker,
                                 const ModelBase& model,
                                 const LatticeBase& lattice )
    {
        // todo
    }




    // ------------------------------  Method routines for dynamic measurements  ---------------------------------
    
    void Methods::measure_dynamic_config_sign( ScalarObs& dynamic_sign, 
                                               const DqmcWalker& walker,
                                               const ModelBase& model,
                                               const LatticeBase& lattice )
    {
        dynamic_sign.tmp_value() += walker.ConfigSign();
        ++dynamic_sign;
    }

} // namespace Measure




//     void Methods::measure_kinetic_energy(Observable<double> &kinetic_energy, Measure &measure, const Model::Hubbard &hubbard) {
//         const int ll = hubbard.ll;
//         for (int t = 0; t < hubbard.lt; ++t) {
//             const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
//             const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

//             double tmp_kinetic_energy = 0.0;
//             for (int x = 0; x < ll; ++x) {
//                 for (int y = 0; y < ll; ++y) {
//                     tmp_kinetic_energy += 2 * hubbard.t * (*hubbard.vec_config_sign)[t] *
//                                         ( gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll))
//                                         + gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)) );
//                 }
//             }
//             kinetic_energy.tmp_value() += tmp_kinetic_energy / hubbard.ls;
//             ++kinetic_energy;
//         }
//     }

//     void Methods::measure_local_spin_corr(Observable<double> &local_spin_corr, Measure &measure, const Model::Hubbard &hubbard) {
//         for (int t = 0; t < hubbard.lt; ++t) {
//             const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
//             const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

//             double tmp_local_spin_corr = 0.0;
//             for (int i = 0; i < hubbard.ls; ++i) {
//                 tmp_local_spin_corr += (*hubbard.vec_config_sign)[t] * ( gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i) );
//             }
//             local_spin_corr.tmp_value() += tmp_local_spin_corr / hubbard.ls;
//             ++local_spin_corr;
//         }
//     }

//     void Methods::measure_momentum_distribution(Observable<double> &momentum_dist, Measure &measure, const Model::Hubbard &hubbard) {
//         const int ll = hubbard.ll;
//         for (int t = 0; t < hubbard.lt; ++t) {
//             const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
//             const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

//             double tmp_momentum_dist = 0.0;
//             // base point
//             for (int xi = 0; xi < ll; ++xi) {
//                 for (int yi = 0; yi < ll; ++yi) {
//                     const int i = xi + ll * yi;
//                     // displacement
//                     for (int dx = 0; dx < ll; ++dx) {
//                         for (int dy = 0; dy < ll; ++dy) {
//                             const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
//                             const Eigen::Vector2d r(dx, dy);
//                             tmp_momentum_dist += cos(-r.dot(measure.q)) * (gu(j, i) + gd(j, i));
//                         }
//                     }
//                 }
//             }
//             momentum_dist.tmp_value() += (*hubbard.vec_config_sign)[t] * (1 - 0.5 * tmp_momentum_dist / hubbard.ls);
//             ++momentum_dist;
//         }
//     }

//     void Methods::measure_spin_density_structure_factor(Observable<double> &sdw_factor, Measure &measure, const Model::Hubbard &hubbard) {
//         const int ll = hubbard.ll;
//         const int ls = hubbard.ls;
//         for (int t = 0; t < hubbard.lt; ++t) {
//             const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
//             const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

//             // g(i,j)  = < c_i * c^+_j >
//             // gc(i,j) = < c^+_i * c_j >
//             const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - gu.transpose();
//             const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - gd.transpose();

//             // loop for site i and average
//             double tmp_sdw = 0.0;
//             for (int xi = 0; xi < ll; ++xi) {
//                 for (int yi = 0; yi < ll; ++yi) {
//                     const int i = xi + ll * yi;
//                     // displacement
//                     for (int dx = 0; dx < ll; ++dx) {
//                         for (int dy = 0; dy < ll; ++dy) {
//                             const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
//                             const Eigen::Vector2d r(dx, dy);
//                             const double factor = (*hubbard.vec_config_sign)[t] * cos(-r.dot(measure.q));
//                             // factor 1/4 comes from spin 1/2
//                             tmp_sdw += 0.25 * factor * ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
//                                                          + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
//                                                          - gdc(i, i) * guc(j, j) - guc(i, i) * gdc(j, j) );
//                         }
//                     }
//                 }
//             }
//             sdw_factor.tmp_value() += tmp_sdw / hubbard.ls / hubbard.ls;
//             ++sdw_factor;
//         }
//     }

//     void Methods::measure_charge_density_structure_factor(Observable<double> &cdw_factor, Measure &measure, const Model::Hubbard &hubbard) {
//         const int ll = hubbard.ll;
//         const int ls = hubbard.ls;
//         for (int t = 0; t < hubbard.lt; ++t) {
//             const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
//             const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

//             // g(i,j)  = < c_i * c^+_j >
//             // gc(i,j) = < c^+_i * c_j >
//             const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - gu.transpose();
//             const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - gd.transpose();

//             // loop for site i and average
//             double tmp_cdw = 0.0;
//             for (int xi = 0; xi < ll; ++xi) {
//                 for (int yi = 0; yi < ll; ++yi) {
//                     const int i = xi + ll * yi;
//                     // displacement
//                     for (int dx = 0; dx < ll; ++dx) {
//                         for (int dy = 0; dy < ll; ++dy) {
//                             const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
//                             const Eigen::Vector2d r(dx, dy);
//                             const double factor = (*hubbard.vec_config_sign)[t] * cos(-r.dot(measure.q));
//                             tmp_cdw += factor * ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
//                                                   + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
//                                                   + gdc(i, i) * guc(j, j) + guc(i, i) * gdc(j, j) );
//                         }
//                     }
//                 }
//             }
//             cdw_factor.tmp_value() += tmp_cdw / hubbard.ls / hubbard.ls;
//             ++cdw_factor;
//         }
//     }

//     void Methods::measure_s_wave_pairing_corr(Observable<double> &s_wave_pairing, Measure &measure, const Model::Hubbard &hubbard) {
//         const int ll = hubbard.ll;
//         const int ls = hubbard.ls;
//         for (int t = 0; t < hubbard.lt; ++t) {
//             // The s-wave pairing order parameter is defined as \Delta = 1/\sqrt(N) \sum_i c_up_i * c_dn_i
//             // Accordingly, the s-wave pairing correlation function P_s is defined as  
//             //   P_s = 1/2 ( \Delta^+ * \Delta + h.c. )
//             //       = 1/N \sum_ij ( (\delta_ij - G_up(t,t)_ji) * (\delta_ij - G_dn(t,t)_ji) )
//             // which is a extensive quantity. The 1/2 prefactor in definition of P_s cancels the duplicated counting of ij.
            
//             // g(i,j)  = < c_i * c^+_j >
//             // gc(i,j) = < c^+_i * c_j >
//             const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - (*hubbard.vec_green_tt_up)[t].transpose();
//             const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - (*hubbard.vec_green_tt_dn)[t].transpose();

//             // loop for site i and average
//             double tmp_s_wave_pairing = 0.0;
//             for (int xi = 0; xi < ll; ++xi) {
//                 for (int yi = 0; yi < ll; ++yi) {
//                     const int i = xi + ll * yi;
//                     // displacement
//                     for (int dx = 0; dx < ll; ++dx) {
//                         for (int dy = 0; dy < ll; ++dy) {
//                             const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
//                             const double factor = (*hubbard.vec_config_sign)[t];
//                             tmp_s_wave_pairing += (*hubbard.vec_config_sign)[t] * ( guc(i, j)* gdc(i, j) );
//                         }
//                     }
//                 }
//             }
//             // entensive quantity
//             s_wave_pairing.tmp_value() += tmp_s_wave_pairing / hubbard.ls;
//             ++s_wave_pairing;
//         }
//     }



//     /** Time-displaced ( Dynamical ) Measurements */


//     void Methods::measure_greens_functions(Observable<Eigen::MatrixXd> &greens_functions, Measure &measure, const Model::Hubbard &hubbard) {
//         for (int t = 0; t < hubbard.lt; ++t) {
//             // factor 1/2 comes from two degenerate spin states
//             const Eigen::MatrixXd gt0 = ( t == 0 )?
//                     0.5 * ((*hubbard.vec_green_tt_up)[hubbard.lt-1] + (*hubbard.vec_green_tt_dn)[hubbard.lt-1])
//                   : 0.5 * ((*hubbard.vec_green_t0_up)[t-1] + (*hubbard.vec_green_t0_dn)[t-1]);

//             // base point i
//             for (int xi = 0; xi < hubbard.ll; ++xi) {
//                 for (int yi = 0; yi < hubbard.ll; ++yi) {
//                     const int i = xi + hubbard.ll * yi;
//                     // displacement
//                     for (int dx = 0; dx < hubbard.ll; ++dx) {
//                         for (int dy = 0; dy < hubbard.ll; ++dy) {
//                             const int j = (xi + dx) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
//                             const Eigen::Vector2d r(dx, dy);
//                             // loop for momentum in qlist
//                             for (int idq = 0; idq < measure.q_list.size(); ++idq) {
//                                 greens_functions.tmp_value()(idq, t) += hubbard.config_sign * cos(-r.dot(measure.q_list[idq])) * gt0(j, i) / hubbard.ls;
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//         ++greens_functions;
//     }

//     void Methods::measure_density_of_states(Observable<Eigen::VectorXd> &density_of_states, Measure &measure, const Model::Hubbard &hubbard) {
//         for (int t = 0; t < hubbard.lt; ++t) {
//             // spin degenerate model
//             const Eigen::MatrixXd gt0 = ( t == 0 )?
//                 0.5 * ((*hubbard.vec_green_tt_up)[hubbard.lt-1] + (*hubbard.vec_green_tt_dn)[hubbard.lt-1])
//               : 0.5 * ((*hubbard.vec_green_t0_up)[t-1] + (*hubbard.vec_green_t0_dn)[t-1]);
//             density_of_states.tmp_value()(t) += hubbard.config_sign * gt0.trace() / hubbard.ls;
//         }
//         ++density_of_states;
//     }

//     void Methods::measure_superfluid_stiffness(Observable<double> &superfluid_stiffness, Measure &measure, const Model::Hubbard &hubbard) {
//         // momentum qx and qy
//         const Eigen::VectorXd qx = ( Eigen::VectorXd(2) << 2 * M_PI / hubbard.ll, 0.0 ).finished();
//         const Eigen::VectorXd qy = ( Eigen::VectorXd(2) << 0.0, 2 * M_PI / hubbard.ll ).finished();
        
//         // Superfluid stiffness \rho_s = 1/4 * ( Gamma^L - Gamma^T )
//         double tmp_rho_s = 0.0;
//         Eigen::MatrixXd gt0_up, g0t_up, gtt_up;
//         Eigen::MatrixXd gt0_dn, g0t_dn, gtt_dn;
//         const Eigen::MatrixXd g00_up = (*hubbard.vec_green_tt_up)[hubbard.lt-1]; 
//         const Eigen::MatrixXd g00_dn = (*hubbard.vec_green_tt_dn)[hubbard.lt-1]; 

//         for (int l = 0; l < hubbard.lt; ++l) {
//             const int tau = (l == 0)? hubbard.lt-1 : l-1;
//             gt0_up = (*hubbard.vec_green_t0_up)[tau];
//             g0t_up = (*hubbard.vec_green_0t_up)[tau];
//             gtt_up = (*hubbard.vec_green_tt_up)[tau];
//             gt0_dn = (*hubbard.vec_green_t0_dn)[tau];
//             g0t_dn = (*hubbard.vec_green_0t_dn)[tau];
//             gtt_dn = (*hubbard.vec_green_tt_dn)[tau];

//             // space point i is chosen as our base point, which is going to be averaged
//             for (int xi = 0; xi < hubbard.ll; ++xi) {
//                 for (int yi = 0; yi < hubbard.ll; ++yi) {
//                     const int i = xi + hubbard.ll * yi;
//                     const int ipx = (xi + 1) % hubbard.ll + hubbard.ll * yi;

//                     // displacement
//                     for (int dx = 0; dx < hubbard.ll; ++dx) {
//                         for (int dy = 0; dy < hubbard.ll; ++dy) {
//                             /* for a given site l and time-slice tau
//                             * the current-current correlation Jx-Jx: \Gamma_xx (l, \tau) = < jx(l, \tau) * jx(0, 0) > */
//                             const int j = (xi + dx) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
//                             const int jpx = (xi + dx + 1) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
//                             const Eigen::Vector2d r(dx, dy);
//                             const double factor = hubbard.config_sign * (cos(r.dot(qx)) - cos(r.dot(qy)));

//                             tmp_rho_s += hubbard.t * hubbard.t * factor * (
//                                     // uncorrelated part
//                                     - ( gtt_up(j, jpx) - gtt_up(jpx, j) + gtt_dn(j, jpx) - gtt_dn(jpx, j) ) *
//                                       ( g00_up(i, ipx) - g00_up(ipx, i) + g00_dn(i, ipx) - g00_dn(ipx, i) )
                                    
//                                     // correlated part
//                                     - g0t_up(ipx, jpx) * gt0_up(j, i) - g0t_dn(ipx, jpx) * gt0_dn(j, i)
//                                     + g0t_up(i, jpx) * gt0_up(j, ipx) + g0t_dn(i, jpx) * gt0_dn(j, ipx)
//                                     + g0t_up(ipx, j) * gt0_up(jpx, i) + g0t_dn(ipx, j) * gt0_dn(jpx, i)
//                                     - g0t_up(i, j) * gt0_up(jpx, ipx) - g0t_dn(i, j) * gt0_dn(jpx, ipx) );
//                         }
//                     }
//                 }
//             }
//         }
//         // average over base point i
//         // the 1/4 prefactor is due to Cooper pairs with charge 2
//         // see https://arxiv.org/pdf/1912.08848.pdf
//         superfluid_stiffness.tmp_value() += 0.25 * tmp_rho_s / hubbard.ls / hubbard.ls;
//         ++superfluid_stiffness;
//     }

// } // namespace Measure