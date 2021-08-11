#include "CheckerBoard.h"
#include "hubbard.h"


bool CheckerBoard::is_checker_board() const {
    return this->is_checkerboard;
}

void CheckerBoard::check_checker_board() {
    if (is_checkerboard && ll % 2 != 0) {
        is_checkerboard = false;
        std::cerr << "Checkerboard break-up only supported for even lattice sizes." << std::endl;
        std::cerr << "  Simulating with direct multiplication algorithm ..." << std::endl;
    }
}

void CheckerBoard::init_from_model(const Hubbard &hubbard) {
    this->ll = hubbard.ll;
    this->ls = hubbard.ls;
    this->t = hubbard.t;
    this->mu = hubbard.mu;
    this->dtau = hubbard.dtau;
    this->is_checkerboard = hubbard.is_checkerboard;

    exp_dtK.resize(ls, ls);
    inv_exp_dtK.resize(ls, ls);

    check_checker_board();
    if (is_checkerboard) { make_expK_checkerboard();}
    else { make_expK_direct();}
}

void CheckerBoard::make_expK_direct() {
    // hopping matrix K: depends on geometry and hopping
    assert(exp_dtK.cols() == ls && exp_dtK.rows() == ls);
    assert(inv_exp_dtK.cols() == ls && inv_exp_dtK.rows() == ls);

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ls, ls);
    for(int x = 0; x < ll; ++x) {
        for(int y = 0; y < ll; ++y) {
            // chemical potential 'mu' on the diagonal
            K(x + ll * y, x + ll * y) -= mu;
            K(x + ll * y, ((x + 1) % ll) + ll * y) -= t;
            K(((x + 1) % ll) + ll * y, x + ll * y) -= t;
            K(x + ll * y, x + ll * ((y + 1) % ll)) -= t;
            K(x + ll * ((y + 1) % ll), x + ll * y) -= t;
        }
    }
    exp_dtK = (-dtau * K).exp();
    inv_exp_dtK = (+dtau * K).exp();
}

void CheckerBoard::make_expK_checkerboard() {
    // construct exponent of hopping in a reduced sub Hilbert-space
    double ch, sh;
    Eigen::Matrix4d reduced_K;

    ch = cosh(dtau * t);
    sh = sinh(dtau * t);
    reduced_K << ch * ch, ch * sh, ch * sh, sh * sh,
                 ch * sh, ch * ch, sh * sh, ch * sh,
                 ch * sh, sh * sh, ch * ch, ch * sh,
                 sh * sh, ch * sh, ch * sh, ch * ch;
    // factor 0.5 comes from the double counting of sites
    exp_dtK_reduced =  exp(0.5 * dtau * mu) * reduced_K;

    ch = cosh(-dtau * t);
    sh = sinh(-dtau * t);
    reduced_K << ch * ch, ch * sh, ch * sh, sh * sh,
                 ch * sh, ch * ch, sh * sh, ch * sh,
                 ch * sh, sh * sh, ch * ch, ch * sh,
                 sh * sh, ch * sh, ch * sh, ch * ch;
    inv_exp_dtK_reduced = exp(-0.5 * dtau * mu) * reduced_K;
}

int xy_to_label(int x, int y, int ll) {
    return (x % ll) + ll * (y % ll);
}

// ref: Max H. Gerlach, 2017
// (x,y) counts the upper left conner of a plaquette.
// plaquette: labeled by site i (x, y), the upper-left conner
//   i ---- j
//   |      |
//   |      |
//   k ---- l
// and the corresponding effective hopping matrix (4*4) reads
//   0.0, 1.0, 1.0, 0.0,
//   1.0, 0.0, 0.0, 1.0,
//   1.0, 0.0, 0.0, 1.0,
//   0.0, 1.0, 1.0, 0.0.

void CheckerBoard::mult_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const {
    assert(A.rows() == ls && A.cols() == ls);

    const int label_xy = xy_to_label(x, y, ll);
    const int label_xy_right = xy_to_label(x+1, y, ll);
    const int label_xy_down = xy_to_label(x, y+1, ll);
    const int label_xy_diagonal = xy_to_label(x+1, y+1, ll);

    Eigen::MatrixXd tmp_mat(4, A.cols());
    tmp_mat.row(0) = A.row(label_xy);
    tmp_mat.row(1) = A.row(label_xy_right);
    tmp_mat.row(2) = A.row(label_xy_down);
    tmp_mat.row(3) = A.row(label_xy_diagonal);

    tmp_mat = exp_dtK_reduced * tmp_mat;

    A.row(label_xy) = tmp_mat.row(0);
    A.row(label_xy_right) = tmp_mat.row(1);
    A.row(label_xy_down) = tmp_mat.row(2);
    A.row(label_xy_diagonal) = tmp_mat.row(3);
}

void CheckerBoard::mult_inv_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const {
    assert(A.rows() == ls && A.cols() == ls);

    const int label_xy = xy_to_label(x, y, ll);
    const int label_xy_right = xy_to_label(x+1, y, ll);
    const int label_xy_down = xy_to_label(x, y+1, ll);
    const int label_xy_diagonal = xy_to_label(x+1, y+1, ll);

    Eigen::MatrixXd tmp_mat(4, A.cols());
    tmp_mat.row(0) = A.row(label_xy);
    tmp_mat.row(1) = A.row(label_xy_right);
    tmp_mat.row(2) = A.row(label_xy_down);
    tmp_mat.row(3) = A.row(label_xy_diagonal);

    tmp_mat = inv_exp_dtK_reduced * tmp_mat;

    A.row(label_xy) = tmp_mat.row(0);
    A.row(label_xy_right) = tmp_mat.row(1);
    A.row(label_xy_down) = tmp_mat.row(2);
    A.row(label_xy_diagonal) = tmp_mat.row(3);
}

void CheckerBoard::mult_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const {
    assert(A.rows() == ls && A.cols() == ls);

    const int label_xy = xy_to_label(x, y, ll);
    const int label_xy_right = xy_to_label(x+1, y, ll);
    const int label_xy_down = xy_to_label(x, y+1, ll);
    const int label_xy_diagonal = xy_to_label(x+1, y+1, ll);

    Eigen::MatrixXd tmp_mat(A.rows(), 4);
    tmp_mat.col(0) = A.col(label_xy);
    tmp_mat.col(1) = A.col(label_xy_right);
    tmp_mat.col(2) = A.col(label_xy_down);
    tmp_mat.col(3) = A.col(label_xy_diagonal);

    tmp_mat = tmp_mat * exp_dtK_reduced;

    A.col(label_xy) = tmp_mat.col(0);
    A.col(label_xy_right) = tmp_mat.col(1);
    A.col(label_xy_down) = tmp_mat.col(2);
    A.col(label_xy_diagonal) = tmp_mat.col(3);
}

void CheckerBoard::mult_inv_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const {
    assert(A.rows() == ls && A.cols() == ls);

    const int label_xy = xy_to_label(x, y, ll);
    const int label_xy_right = xy_to_label(x+1, y, ll);
    const int label_xy_down = xy_to_label(x, y+1, ll);
    const int label_xy_diagonal = xy_to_label(x+1, y+1, ll);

    Eigen::MatrixXd tmp_mat(A.rows(), 4);
    tmp_mat.col(0) = A.col(label_xy);
    tmp_mat.col(1) = A.col(label_xy_right);
    tmp_mat.col(2) = A.col(label_xy_down);
    tmp_mat.col(3) = A.col(label_xy_diagonal);

    tmp_mat = tmp_mat * inv_exp_dtK_reduced;

    A.col(label_xy) = tmp_mat.col(0);
    A.col(label_xy_right) = tmp_mat.col(1);
    A.col(label_xy_down) = tmp_mat.col(2);
    A.col(label_xy_diagonal) = tmp_mat.col(3);
}

void CheckerBoard::mult_expK_from_left(Eigen::MatrixXd &A) const {
    if (is_checkerboard) {
        // checkerboard break-up only supported for even lattice sizes
        assert(ll % 2 == 0);

        // board B
        for (int x = 1; x < ll; x+=2) {
            for (int y = 1; y < ll; y+=2) {
                mult_expK_plaquette_from_left(A, x, y);
            }
        }
        // board A
        for (int x = 0; x < ll; x+=2) {
            for (int y = 0; y < ll; y+=2) {
                mult_expK_plaquette_from_left(A, x, y);
            }
        }
    }
    else { A = exp_dtK * A;}
}

void CheckerBoard::mult_inv_expK_from_left(Eigen::MatrixXd &A) const {
    if (is_checkerboard) {
        // checkerboard break-up only supported for even lattice sizes
        assert(ll % 2 == 0);

        // board A
        for (int x = 0; x < ll; x+=2) {
            for (int y = 0; y < ll; y+=2) {
                mult_inv_expK_plaquette_from_left(A, x, y);
            }
        }
        // board B
        for (int x = 1; x < ll; x+=2) {
            for (int y = 1; y < ll; y+=2) {
                mult_inv_expK_plaquette_from_left(A, x, y);
            }
        }
    }
    else { A = inv_exp_dtK * A;}
}

void CheckerBoard::mult_expK_from_right(Eigen::MatrixXd &A) const {
    if (is_checkerboard) {
        // checkerboard break-up only supported for even lattice sizes
        assert(ll % 2 == 0);

        // board A
        for (int x = 0; x < ll; x+=2) {
            for (int y = 0; y < ll; y+=2) {
                mult_expK_plaquette_from_right(A, x, y);
            }
        }
        // board B
        for (int x = 1; x < ll; x+=2) {
            for (int y = 1; y < ll; y+=2) {
                mult_expK_plaquette_from_right(A, x, y);
            }
        }
    }
    else { A = A * exp_dtK;}
}

void CheckerBoard::mult_inv_expK_from_right(Eigen::MatrixXd &A) const {
    if (is_checkerboard) {
        // checkerboard break-up only supported for even lattice sizes
        assert(ll % 2 == 0);

        // board B
        for (int x = 1; x < ll; x+=2) {
            for (int y = 1; y < ll; y+=2) {
                mult_inv_expK_plaquette_from_right(A, x, y);
            }
        }
        // board A
        for (int x = 0; x < ll; x+=2) {
            for (int y = 0; y < ll; y+=2) {
                mult_inv_expK_plaquette_from_right(A, x, y);
            }
        }
    }
    else { A = A * inv_exp_dtK;}
}

void CheckerBoard::mult_trans_expK_from_left(Eigen::MatrixXd &A) const {
    if (is_checkerboard) {
        // FIXME: rewrite in a direct way
        A.transposeInPlace();
        mult_expK_from_right(A);
        A.transposeInPlace();
    }
    else { A = exp_dtK.transpose() * A;}
}
