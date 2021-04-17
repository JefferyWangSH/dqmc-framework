#include "hubbard.h"
#include "StableGreens.hpp"
#include "ProgressBar.hpp"

/** TODO:
 *   1. cyclic sweeping (done)
 *   2. chemical potential mu (done)
 *   3. stable multiplication of ill-conditioned matrices (done)
 *   4. check greens (done, unnecessary)
 *   6. detQMC class (done)
 *   7. check-board decomposition (missing)
 *   8. time-displaced dynamic measurements: Green_t_0_up/dn (done)
 *   9. attractive interaction U < 0 (done)
 *   10. ...
 */

Hubbard::Hubbard(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap)
{

    this->ll = ll;
    this->ls = ll * ll;
    this->lt = lt;

    this->beta = beta;
    this->dtau = beta / lt;
    this->Uint = Uint;
    this->alpha = acosh(exp(0.5 * dtau * abs(Uint)));
    this->u_is_attractive = (Uint <= 0);

    // unnecessary, only used to generate expK
    this->t = t;
    this->mu = mu;

    this->nwrap = nwrap;
    this->current_tau = 0;

    // resize matrices and svdStacks
    s.resize(ls, lt);
    exppdtK.resize(ls, ls);
    expmdtK.resize(ls, ls);
    GreenU.resize(ls, ls);
    GreenD.resize(ls, ls);
    Green_t0_up.resize(ls, ls);
    Green_t0_dn.resize(ls, ls);

    stackLeftU.resize(ls, lt);
    stackLeftD.resize(ls, lt);
    stackRightU.resize(ls, lt);
    stackRightD.resize(ls, lt);

    vecGreenU.reserve(lt);
    vecGreenD.reserve(lt);
    vecGreen_t0_up.reserve(lt);
    vecGreen_t0_dn.reserve(lt);

    for (int l = 0; l < lt; ++l) {
        vecGreenU.emplace_back(ls, ls);
        vecGreenD.emplace_back(ls, ls);
        vecGreen_t0_up.emplace_back(ls, ls);
        vecGreen_t0_dn.emplace_back(ls, ls);
    }

    // set field config to random
    initRandom();

    // compute exp of Kinetic matrix K
    make_expdtK();

    // initialize udv stacks for sweep use, stabilize every nwrap slices
    initStacks(nwrap);
}

void Hubbard::initRandom() {
    // set field configuration to random
    assert(s.rows() == ls && s.cols() == lt);

    std::bernoulli_distribution dist(0.5);
    for(int i = 0; i < ls; ++i) {
        for(int l = 0; l < lt; ++l) {
            s(i,l) = dist(gen)? +1.0:-1.0;
        }
    }
}

void Hubbard::make_expdtK() {
    // kinetic matrix K: depends on geometry and hopping
    assert(expmdtK.cols() == ls && expmdtK.rows() == ls);
    assert(exppdtK.cols() == ls && exppdtK.rows() == ls);

    matXd K = matXd::Zero(ls, ls);
    for(int x = 0; x < ll; ++x) {
        for(int y = 0; y < ll; ++y) {
            K(x + ll*y, ((x+1)%ll) + ll*y) = -t;
            K(((x+1)%ll) + ll*y, x + ll*y) = -t;
            K(x + ll*y, x + ll*((y+1)%ll)) = -t;
            K(x + ll*((y+1)%ll), x + ll*y) = -t;
            if (x==y) { K(x, y) = - mu;}
        }
    }
    expmdtK = (-dtau * K).exp();
    exppdtK = (+dtau * K).exp();
}

matXd Hubbard::make_Bl(int l, int sigma) {
    /*
     *  Compute B matrix at time slice l for given spin-1/2 state.
     *  definition of B(l, sigma):
     *  B_l = exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K)
     */
    assert( l >= 0 && l <= lt);
    assert( sigma == 1 || sigma == -1);

    matXd r = expmdtK;
    const int tau = (l==0)? lt-1 : l-1;
    const int Sigma = (u_is_attractive)? 1 : sigma;
    for (int i=0; i < ls; ++i) {
        r.row(i) *= exp(+ Sigma * alpha * s(i,tau));
    }
    return r;
}

void Hubbard::multB_fromL(matXd& A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the left by B_l
     *  A ->  B_l * A = exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K) * A
     *  Matrix A is changed in place.
     */
    assert(A.rows() == ls && A.cols() == ls);
    assert(l >= 0 && l <= lt);
    assert( sigma == 1 || sigma == -1);

    const int tau = (l==0)? lt-1 : l-1;
    const int Sigma = (u_is_attractive)? 1 : sigma;
    A = expmdtK * A;
    for (int i=0; i < ls; ++i) {
        A.row(i) *= exp(+ Sigma * alpha * s(i,tau));
    }
}

void Hubbard::multB_fromR(matXd& A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the right by B_l
     *  A ->  A * B_l = A * exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K)
     *  Matrix A is changed in place.
     */
    assert(A.rows() == ls && A.cols() == ls);
    assert(l >= 0 && l <= lt);
    assert( sigma == 1 || sigma == -1);

    const int tau = (l==0)? lt-1 : l-1;
    const int Sigma = (u_is_attractive)? 1 : sigma;
    for (int i=0; i < ls; ++i) {
        A.col(i) *= exp(+ Sigma * alpha * s(i,tau));
    }
    A = A * expmdtK;
}

void Hubbard::multinvB_fromL(matXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the left by B_l^{-1}
     *  A -> B_{l}^{-1} * A = exp(+\\Delta \\tau K) * exp(+\\Delta \\tau V^{\\sigma}(l))  * A
     *  Matrix A is changed in place.
     */
    assert(A.rows() == ls && A.cols() == ls);
    assert(l >= 0 && l <= lt);
    assert( sigma == 1 || sigma == -1);

    const int tau = (l==0)? lt-1 : l-1;
    const int Sigma = (u_is_attractive)? 1 : sigma;
    for (int i=0; i < ls; ++i) {
        A.row(i) *= exp(- Sigma * alpha * s(i,tau));
    }
    A = exppdtK * A;
}

void Hubbard::multinvB_fromR(matXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the right by B_l^{-1}
     *  A -> A * B_{l}^{-1} = A * exp(+\\Delta \\tau K) * exp(+\\Delta \\tau V^{\\sigma}(l))
     *  Matrix A is changed in place.
     */
    assert(A.rows() == ls && A.cols() == ls);
    assert(l >= 0 && l <= lt);
    assert( sigma == 1 || sigma == -1);

    const int tau = (l==0)? lt-1 : l-1;
    const int Sigma = (u_is_attractive)? 1 : sigma;
    A = A * exppdtK;
    for (int i=0; i < ls; ++i) {
        A.col(i) *= exp(- Sigma * alpha * s(i,tau));
    }
}

void Hubbard::Metropolis_update(int l) {
    /*
     * Update a HS field at space-time position (i,l) for all i with Metropolis
     * probability, and - if the update is accepted - perform a
     * in-place update of the Green's function.
     * Record the updated green's function at the end of function.
     */
    assert(current_tau == l);
    assert(l >= 0 && l <= lt);

    const int tau = (l==0)? lt-1 : l-1;
    for (int i = 0; i < ls; ++i) {
        // radio of flipping aux field s(i,l)
        double p;
        if (!u_is_attractive) {
            p = (1 + (1 - GreenU(i,i)) * (exp(-2 * alpha * s(i,tau)) - 1))
                       * (1 + (1 - GreenD(i,i)) * (exp(+2 * alpha * s(i,tau)) - 1));
        }
        else {
            p = exp(2 * alpha * s(i, tau))  * (1 + (1 - GreenU(i,i)) * (exp(-2 * alpha * s(i,tau)) - 1))
                       * (1 + (1 - GreenD(i,i)) * (exp(-2 * alpha * s(i,tau)) - 1));
        }

        if(std::bernoulli_distribution(std::min(1.0, p))(gen)) {

            /** reference:
             *  Quantum Monte Carlo Methods (Algorithms for Lattice Models) Determinant method
             *  Here we use the sparseness of matrix \delta */
            // update greens function (which is wrapped such that the update is at timeslice 0 of g)
            // with a number of arithmetic operations proportional to N^2
            double factorU = (exp(-2*alpha*s(i,tau))-1)/(1 + (1-GreenU(i,i))*(exp(-2*alpha*s(i,tau))-1));
            GreenU -= factorU * GreenU.col(i) * (vecXd::Unit(ls,i).transpose() - GreenU.row(i));

            double factorD = (!u_is_attractive)? (exp(+2*alpha*s(i,tau))-1)/(1 + (1-GreenD(i,i))*(exp(+2*alpha*s(i,tau))-1)) : factorU;
            GreenD -= factorD * GreenD.col(i) * (vecXd::Unit(ls,i).transpose() - GreenD.row(i));

            // flip aux field
            s(i,tau) = -s(i,tau);
        }
    }

    // record equal-time green function of current time
    vecGreenU[tau] = GreenU;
    vecGreenD[tau] = GreenD;
}

void Hubbard::wrap_north(int l) {
    /*
     * Propagate the Green's function from the current time slice l
     * upward to the time slice l+1:
     * G(l+1) = B_{l+1} G(l) B_{l+1}^{-1}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert(l >= 0 && l <= lt);

    multB_fromL(GreenU, (l+1)%lt, +1);
    multinvB_fromR(GreenU, (l+1)%lt, +1);
    multB_fromL(GreenD, (l+1)%lt, -1);
    multinvB_fromR(GreenD, (l+1)%lt, -1);
}

void Hubbard::wrap_south(int l) {
    /*
     * Propagate the Green's function from the current time slice l
     * downward to the time slice l-1:
     * G(l-1) = B_{l}^{-1} G(l) B_{l}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert(l >= 0 && l <= lt);

    multB_fromR(GreenU, l, +1);
    multinvB_fromL(GreenU, l, +1);
    multB_fromR(GreenD, l, -1);
    multinvB_fromL(GreenD, l, -1);
}

void Hubbard::initStacks(int istab) {
    /*
     *  initialize udv stacks for sweep use
     *  sweep process will start from 0 to beta, so we initialize stackRight here.
     *  stabilize the process every istab steps
     */
    assert(stackLeftU.empty() && stackLeftD.empty());
    assert(stackRightU.empty() && stackRightD.empty());

    matXd tmpU = matXd::Identity(ls, ls);
    matXd tmpD = matXd::Identity(ls, ls);

    // initial udv stacks for sweep use
    for (int l = lt; l >= 1; --l) {
        tmpU = make_Bl(l, +1).transpose() * tmpU;
        tmpD = make_Bl(l, -1).transpose() * tmpD;
        // stabilize every istab steps with svd decomposition
        if ((l-1) % istab == 0) {
            stackRightU.push(tmpU);
            stackRightD.push(tmpD);
            tmpU = matXd::Identity(ls, ls);
            tmpD = matXd::Identity(ls, ls);
        }
    }

    // initialize green function at l = 0
    compute_Green(stackLeftU, stackRightU, GreenU, Green_t0_up, false);
    compute_Green(stackLeftD, stackRightD, GreenD, Green_t0_dn, false);
}

void Hubbard::sweep_0_to_beta(int istab) {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l = 1,2...,lt  flip fields and propagate green's functions
     *  Stabilize every istab time slices
     */
    current_tau++;

    int nlen = (lt % istab == 0)? lt/istab : lt/istab+1;
    assert(current_tau == 1);
    assert(stackLeftU.empty() && stackLeftD.empty());
    assert(stackRightU.len==nlen && stackRightD.len==nlen);

    // temporary matrices
    matXd tmpU = matXd::Identity(ls, ls);
    matXd tmpD = matXd::Identity(ls, ls);

    // sweep up from 0 to beta
    for (int l = 1; l <= lt; ++l) {
        // wrap green function to current time slice l
        wrap_north(l-1);

        // update aux field and record new greens
        Metropolis_update(l);

        tmpU = make_Bl(l, +1) * tmpU;
        tmpD = make_Bl(l, -1) * tmpD;

        if (l % istab == 0 || l == lt) {
            // wrap greens function
            stackRightU.pop();
            stackRightD.pop();
            stackLeftU.push(tmpU);
            stackLeftD.push(tmpD);

            // compute fresh greens every istab steps: g = (1 + stackLeft * stackRight^T)^-1
            // stackLeft = B(l-1) *...* B(0)
            // stackRight = B(l)^T *...* B(L-1)^T
            compute_Green(stackLeftU, stackRightU, GreenU, Green_t0_up, false);
            compute_Green(stackLeftD, stackRightD, GreenD, Green_t0_dn, false);

            tmpU = matXd::Identity(ls, ls);
            tmpD = matXd::Identity(ls, ls);
        }

        // finally stop at l = lt + 1
        current_tau++;
    }

    // end with fresh green functions
    vecGreenU[lt-1] = GreenU;
    vecGreenD[lt-1] = GreenD;
}

void Hubbard::sweep_beta_to_0(int istab) {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l=lt,lt-1,...,1  flip fields and propagate green's functions
     *  Stabilize every istab time slices
     */
    current_tau--;

    int nlen = (lt % istab == 0)? lt/istab : lt/istab+1;
    assert(current_tau == lt);
    assert(stackRightU.empty() && stackRightD.empty());
    assert(stackLeftU.len==nlen && stackLeftD.len==nlen);

    // temporary matrices
    matXd tmpU = matXd::Identity(ls, ls);
    matXd tmpD = matXd::Identity(ls, ls);

    // sweep down from beta to 0
    for (int l = lt; l >= 1; --l) {

        // update aux field and record new greens
        Metropolis_update(l);

        tmpU = make_Bl(l, +1).transpose() * tmpU;
        tmpD = make_Bl(l, -1).transpose() * tmpD;

        if ((l-1) % istab == 0) {
            // update udv stacks
            stackLeftU.pop();
            stackLeftD.pop();
            stackRightU.push(tmpU);
            stackRightD.push(tmpD);

            compute_Green(stackLeftU, stackRightU, GreenU, Green_t0_up, false);
            compute_Green(stackLeftD, stackRightD, GreenD, Green_t0_dn, false);

            tmpU = matXd::Identity(ls, ls);
            tmpD = matXd::Identity(ls, ls);
        }
        else
            wrap_south(l);

        current_tau--;
    }

    // end with fresh green functions
    vecGreenU[lt-1] = GreenU;
    vecGreenD[lt-1] = GreenD;
}

void Hubbard::sweep_0_to_beta_displaced(int istab) {
    /*
     *  Calculate time-displaced green function, while the aux field remains unchanged.
     *  For l = 1,2...,lt, recompute SvdStacks every istab time slices.
     *  Data is to be stored in vecGreen_t0_up and vecGreen_t0_dn.
     */
    current_tau++;

    int nlen = (lt % istab == 0)? lt/istab : lt/istab+1;
    assert(current_tau == 1);
    assert(stackLeftU.empty() && stackLeftD.empty());
    assert(stackRightU.len==nlen && stackRightD.len==nlen);

    // initialize: at l = 0, gt0 = g00
    Green_t0_up = GreenU;
    Green_t0_dn = GreenD;

    // temporary matrices
    matXd tmpU = matXd::Identity(ls, ls);
    matXd tmpD = matXd::Identity(ls, ls);

    // sweep up from 0 to beta
    for (int l = 1; l <= lt; ++l) {
        
        // calculate and record time-displaced green functions at different time slices
        Green_t0_up = make_Bl(l, +1) * Green_t0_up;
        Green_t0_dn = make_Bl(l, -1) * Green_t0_dn;
        vecGreen_t0_up[l-1] = Green_t0_up;
        vecGreen_t0_dn[l-1] = Green_t0_dn;

        tmpU = make_Bl(l, +1) * tmpU;
        tmpD = make_Bl(l, -1) * tmpD;

        if (l % istab == 0 || l == lt) {
            // wrap greens function
            stackRightU.pop();
            stackRightD.pop();
            stackLeftU.push(tmpU);
            stackLeftD.push(tmpD);

            // compute fresh greens every istab steps
            // stackLeft = B(l-1) *...* B(0)
            // stackRight = B(l)^T *...* B(L-1)^T
            compute_Green(stackLeftU, stackRightU, GreenU, Green_t0_up, true);
            compute_Green(stackLeftD, stackRightD, GreenD, Green_t0_dn, true);

            vecGreen_t0_up[l-1] = Green_t0_up;
            vecGreen_t0_dn[l-1] = Green_t0_dn;

            tmpU = matXd::Identity(ls, ls);
            tmpD = matXd::Identity(ls, ls);
        }

        // finally stop at l = lt + 1
        current_tau++;
    }
}
