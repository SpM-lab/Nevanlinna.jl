#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <fstream>
#include <gmpxx.h>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>


//precision class is used to define typenames
//template T can be any precision type, e.g. double or mpf_class
template <class T> 
class precision_ {
protected:
    using nev_real = T;
    using nev_complex = std::complex<T>;
    using nev_complex_vector = std::vector<nev_complex>;
    using nev_complex_matrix = Eigen::Matrix <nev_complex, Eigen::Dynamic, Eigen::Dynamic>;
    using nev_complex_matrix_vector = std::vector<nev_complex_matrix>;
};


//Matsubara data storage (theta values)
template <class T>
class imag_domain_data : precision_<T> {
private:
    using typename precision_<T>::nev_real;
    using typename precision_<T>::nev_complex;
    using typename precision_<T>::nev_complex_vector;
public:
    //calculate theta (G -> NG -> theta) and store Matsubara frequencies and theta
    imag_domain_data (std::string ifile, int imag_num) : N_imag_(imag_num) {
        std::ifstream ifs(ifile);
        val_.resize(N_imag_);
        freq_.resize(N_imag_);
        nev_real freq, re, im;
        nev_complex I {0., 1.};
        for (int i = 0; i < N_imag_; i++) {
            ifs >> freq >> re >> im;
            nev_complex val = nev_complex{-re, -im}; //minus signs to transform G to NG
            freq_[i] = nev_complex{0., freq};
            val_[i] = (val - I) / (val + I); //Mobius transform from NG to theta
        }
        //reverse input frequency order (decreasing then) and the corresponding thetas, 
        //which tests to be the most robust interpolation order with Schur algorithm
        std::reverse(freq_.begin(),freq_.end());
        std::reverse(val_.begin(), val_.end());
    }
    //number of Matsubara points
    int N_imag() const { return N_imag_; }
    //contractive interpolant theta values at Matsubara points
    const nev_complex_vector &val() const { return val_; }
    //Matsubra frequencies
    const nev_complex_vector &freq() const { return freq_; }
private:
    int N_imag_;
    nev_complex_vector val_;
    nev_complex_vector freq_;
};


//real frequency NG storage, at omega+i*eta
template <class T>
class real_domain_data : precision_<T> {
private:
    using typename precision_<T>::nev_real;
    using typename precision_<T>::nev_complex;
    using typename precision_<T>::nev_complex_vector;
public:
    //calculate and store real frequencies (at omega+i*eta), uniform grid
    //***change N_real_, omega_min, omega_max and eta as needed***
    real_domain_data (std::string ofile) : ofs(ofile), N_real_(6000), omega_min(-10), omega_max(10), eta(0.001) {
        val_.resize(N_real_);
        freq_.resize(N_real_);
        nev_real inter = (omega_max - omega_min) / (N_real_ - 1);
        nev_real temp = omega_min;
        freq_[0] = nev_complex{omega_min, eta};
        for (int i = 1; i < N_real_; i++) {
            temp += inter;
            freq_[i] = nev_complex{temp, eta};
        }
    }
    //number of real frequencies
    int N_real() const { return N_real_; }
    //NG values at real frequencies
    nev_complex_vector &val() { return val_; }
    //real frequencies
    const nev_complex_vector &freq() const { return freq_; }
    //write real frequencies and spectral function A(omega) values to the output file
    void write () {
        for(int i = 0;i < N_real_; i++){
            ofs << std::fixed << std::setprecision(15);
            ofs << freq_[i].real() << " " << 1 / M_PI * val_[i].imag() <<std::endl;
        }
    }
private:
    std::ofstream ofs;
    int N_real_;
    nev_real omega_min;
    nev_real omega_max;
    nev_real eta;
    nev_complex_vector val_;
    nev_complex_vector freq_;
};
