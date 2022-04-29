#include "schur.h"


int main (int argc, char * argv[]) {
    std::string ifile, ofile;
    int imag_num;
    //prompt user for input parameters
    std::cin >> ifile >> imag_num >> ofile;
    //set calculation precision
    mpf_set_default_prec(128);
    //begin evaluation
    Schur<mpf_class> NG(ifile, imag_num, ofile);
    NG.evaluation();
    return 0;
}