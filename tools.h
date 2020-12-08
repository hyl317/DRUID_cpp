#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "Eigen/Dense"
#include "zlib.h"

// wrappers for file I/O
template<typename IO_TYPE>
class FileOrGZ {
  public:
    bool open(const char *filename, const char *mode);
    int getline();
    int printf(const char *format, ...);
    int close();

    static const int INIT_SIZE = 1024 * 50;

    // IO_TYPE is either FILE* or gzFile;
    IO_TYPE fp;

    // for I/O:
    char *buf;
    size_t buf_size;
    size_t buf_len;

  private:
    void alloc_buf();
};

template<typename IO_TYPE>
void FileOrGZ<IO_TYPE>::alloc_buf() {
  buf = (char *) malloc(INIT_SIZE);
  if (buf == NULL) {
    fprintf(stderr, "ERROR: out of memory\n");
    exit(1);
  }
  buf_size = INIT_SIZE;
  buf_len = 0;
}

template<> bool FileOrGZ<FILE *>::open(const char *filename, const char *mode);
template<> bool FileOrGZ<gzFile>::open(const char *filename, const char *mode);
template<> int FileOrGZ<FILE *>::getline();
template<> int FileOrGZ<gzFile>::getline();
template<> int FileOrGZ<FILE *>::printf(const char *format, ...);
template<> int FileOrGZ<gzFile>::printf(const char *format, ...);
template<> int FileOrGZ<FILE *>::close();
template<> int FileOrGZ<gzFile>::close();

using ibdMapType = std::map<std::string, std::vector<std::pair<double, double>>*>;
struct Pair{
  double ibd1_tot;
  double ibd2_tot;
  ibdMapType *ibd1_map = new ibdMapType();
  ibdMapType *ibd2_map = new ibdMapType();
};


void print_help();

void parse_command_line(int argc, char **argv, std::string &ibdFile, std::string &bimFile, std::string &NeFile,
        std::string &prefix, int &maxDeg, double &minIBD);

Eigen::VectorXd readBimFile(const std::string &bimFile, 
  std::map<std::string, std::map<int, double>*> &snpmap, FileOrGZ<FILE *> &logFile);

int readIBDFile(const std::string &ibdFile, 
  std::map<std::pair<std::string, std::string>, Pair*> &allsegs, FileOrGZ<FILE *> &logFile);

double calc_bkg_sharing(const std::string &NeFile, const Eigen::VectorXd &chrLens, const double &minIBD);
void cumsum_eigen_colvector(const Eigen::VectorXd &source, Eigen::VectorXd &dest);
double log_expectedIBD_beyond_maxGen_given_Ne(const Eigen::VectorXd &Ne, 
    const Eigen::VectorXd &chrLens, int n_p, double minIBD);
inline double logaddexp(double d1, double d2){
    if (d1 > d2){return d1 + log1p(exp(d2-d1));}
    else{return d2 + log1p(exp(d1-d2));}
}


#endif