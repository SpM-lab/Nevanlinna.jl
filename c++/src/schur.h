#include "nevanlinna.h"


template <class T>
class Schur : precision_<T> {
private:
    using typename precision_<T>::nev_complex;
    using typename precision_<T>::nev_complex_vector;
    using typename precision_<T>::nev_complex_matrix;
    using typename precision_<T>::nev_complex_matrix_vector;
public:
    //check Nevanlinna/contractive interpolant existence condition
    Schur (std::string ifile, int imag_num, std::string ofile);
    //evaluation with 0 parametric function 
    void evaluation ();
private:
    int M; //number of Matsubara points
    imag_domain_data <T> imag; //theta values at Matsubara points (G -> NG -> theta)
    real_domain_data <T> real; //real frequency NG storage, at omega + i*eta
    nev_complex_vector phis; //phi_1 to phi_M
    nev_complex_matrix_vector abcds; //intermediate {a, b, c, d}s used to calculate phis
    //memoize intermediate abcds and calculate phis by iteration
    void core ();
};


template <class T>
Schur<T>::Schur (std::string ifile, int imag_num, std::string ofile) : imag(ifile, imag_num), real(ofile)  {
    M = imag_num;
    //fill the Pick matrix
    nev_complex_matrix Pick (M, M);
    nev_complex I {0., 1.};
    nev_complex One {1., 0.};
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            nev_complex freq_i = (imag.freq()[i] - I) / (imag.freq()[i] + I);
            nev_complex freq_j = (imag.freq()[j] - I) / (imag.freq()[j] + I);
            nev_complex one {1., 0.};
            nev_complex nom = one - imag.val()[i] * std::conj(imag.val()[j]);
            nev_complex den = one - freq_i * std::conj(freq_j);
            Pick(i, j) = nom / den;
        }
    }
    //check the positive semi-definiteness of the Pick matrix using Cholesky decomposition
    Eigen::LLT<nev_complex_matrix> lltOfPick(Pick + nev_complex_matrix::Identity(M, M) * 1e-250);
    if(lltOfPick.info() == Eigen::NumericalIssue) 
        std::cerr << "Pick matrix is non positive semi-definite matrix in Schur method." << std::endl;
    else std::cerr << "Pick matrix is positive semi-definite." << std::endl;
}


template <class T>
void Schur<T>::core() {
    phis.resize(M);
    abcds.resize(M);
    phis[0] = imag.val()[0];
    for (int k = 0; k < M; k++) abcds[k] = nev_complex_matrix::Identity(2, 2);
    for (int j = 0; j < M - 1; j++) {
        for (int k = j; k < M; k++) {
            nev_complex_matrix prod(2, 2);
            prod(0, 0) = (imag.freq()[k] - imag.freq()[j]) / (imag.freq()[k] - std::conj(imag.freq()[j]));
            prod(0, 1) = phis[j];
            prod(1, 0) = std::conj(phis[j])*
                        ((imag.freq()[k] - imag.freq()[j]) / (imag.freq()[k] - std::conj(imag.freq()[j])));
            prod(1, 1) = nev_complex{1., 0.};
            abcds[k] *= prod;
        }
        phis[j + 1] = (- abcds[j + 1](1, 1) * imag.val()[j + 1] + abcds[j + 1](0, 1)) /
                        (abcds[j + 1](1, 0) * imag.val()[j + 1] - abcds[j + 1](0, 0));   
    }

    std::stringstream filename;
    filename << "../result/phis.dat";
    std::string result = filename.str();
    std::ofstream f(result.c_str());
    for(int k=0; k<M; k++){
        f << std::fixed;
        f << std::setprecision(100) << std::real(phis[k]) << '\t' << std::imag(phis[k]) << std::endl;
    }
    f.close();
}


template <class T>
void Schur<T>::evaluation () {
    core();
    nev_complex I {0., 1.};
    nev_complex One {1., 0.};

    std::stringstream filename;
    filename << "../result/abcd.dat";
    std::string result = filename.str();
    std::ofstream f(result.c_str());

    for (int i = 0; i < real.N_real(); i++) {
        nev_complex_matrix result = nev_complex_matrix::Identity(2, 2);
        nev_complex z = real.freq()[i];
        for (int j = 0; j < M; j++) {
            nev_complex_matrix prod(2, 2);
            prod(0, 0) = (z - imag.freq()[j]) / (z - std::conj(imag.freq()[j]));
            prod(0, 1) = phis[j];
            prod(1, 0) = std::conj(phis[j])*
                        ((z - imag.freq()[j]) / (z - std::conj(imag.freq()[j])));
            prod(1, 1) = nev_complex{1., 0.};
            result *= prod;
        }
        nev_complex param {0., 0.}; //theta_{M+1}, choose to be constant function 0 here
        nev_complex theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
        //can output "real.freq(), a.real(), a.imag(), ..., d.imag()\n" into a file for optimization convenience
        real.val()[i] = I * (One + theta) / (One - theta); //inverse Mobius transform from theta to NG

        f << std::fixed;
        f << std::setprecision(100) << std::real(result(0,0)) << '\t' << std::imag(result(0,0)) << '\t';
        f << std::setprecision(100) << std::real(result(0,1)) << '\t' << std::imag(result(0,1)) << '\t';
        f << std::setprecision(100) << std::real(result(1,0)) << '\t' << std::imag(result(1,0)) << '\t';
        f << std::setprecision(100) << std::real(result(1,1)) << '\t' << std::imag(result(1,1)) << std::endl;
    }
    real.write();
    f.close();
}
