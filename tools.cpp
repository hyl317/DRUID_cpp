#include <string.h>
#include <set>
#include "tools.h"
#include "assert.h"
#include "stdarg.h"

void print_help(){
  exit(0);
}

void parse_command_line(int argc, char **argv, std::string &ibdFile, std::string &bimFile, 
    std::string &NeFile, std::string &prefix, int &maxDeg, double &minIBD){
    int i = 1;
    for(; i < argc; i++){
        char *token = argv[i];
        if (strcmp(token, "-h") == 0 || strcmp(token, "--help") == 0){
            print_help();
        }else if (strcmp(token, "-i") == 0){
            ibdFile = argv[i+1];
            i++;
        }else if (strcmp(token, "--bim") == 0){
            bimFile = argv[i+1];
            i++;
        }else if (strcmp(token, "--Ne") == 0){
            NeFile = argv[i+1];
            i++;
        }else if (strcmp(token, "-o") == 0){
            prefix = argv[i+1];
            i++;
        }else if (strcmp(token, "--minIBD") == 0){
            minIBD = std::stod(argv[i+1]);
            i++;
        }else if (strcmp(token, "--max") == 0){
            maxDeg = std::stoi(argv[i+1]);
            i++; 
        }else{
            fprintf(stderr, "unrecognized token %s\n", token);
            print_help();
        }
    }
}


Eigen::VectorXd readBimFile(const std::string &bimFile, 
    std::map<std::string, std::map<int, double>*> &snpmap, FileOrGZ<FILE *> &logFile){
  FileOrGZ<FILE *> in;
  bool ret = in.open(bimFile.c_str(), "r");
  if (!ret){
    fprintf(stderr, "cannot open %s\n", bimFile.c_str());
    exit(1);
  }

  while(in.getline() >= 0){
    char chr[50];
    double genetic_pos;
    int bp_pos;
    sscanf(in.buf, "%s %*s %lf %d %*s %*s", chr, &genetic_pos, &bp_pos);
    std::string chr_str = chr;
    if (snpmap.find(chr_str) != snpmap.end()){
      snpmap[chr_str]->insert(std::pair<int, double>(bp_pos, genetic_pos));
    }else{
      snpmap.insert(std::pair<std::string, std::map<int, double>*>(chr_str, new std::map<int, double>()));
      snpmap[chr_str]->insert(std::pair<int, double>(bp_pos, genetic_pos));
    }
  }

  // log some basic info about chromosomes
  logFile.printf("\tNumber of chromosomes: %d\n", snpmap.size());
  int numChr = snpmap.size();
  Eigen::VectorXd chrLens(numChr);
  int counter = 0;
  for (auto it = snpmap.begin(); it != snpmap.end(); it++){
    // C++ container is sorted, so the following is fine
    double snp_start = it->second->begin()->second;
    double snp_end = (--it->second->end())->second;
    chrLens[counter++] = snp_end - snp_start;
  }

  logFile.printf("\ttotal genome length: %lfcM\n", chrLens.sum());
  return chrLens;
}

int readIBDFile(const std::string &ibdFile, 
  std::map<std::pair<std::string, std::string>, Pair*> &allsegs, FileOrGZ<FILE *> &logFile){

  FileOrGZ<gzFile> in;
  bool ret = in.open(ibdFile.c_str(), "r");
  if (!ret){
    fprintf(stderr, "cannot open %s\n", ibdFile.c_str());
    exit(1);
  }

  std::set<std::string> inds;
  while(in.getline() >= 0){
    char id1_[50];
    char id2_[50];
    char chr_[50];
    char ibd12_[5];
    double start, end;
    sscanf(in.buf, "%s %s %s %s %lf %lf", id1_, id2_, chr_, ibd12_, &start, &end);

    std::string id1 = id1_;
    std::string id2 = id2_;
    std::string chr = chr_;
    std::string ibd12 = ibd12_;
    double segLen = end - start;
    std::pair<std::string, std::string> pair;
    if (id1 < id2){pair = std::make_pair(id1, id2);}
    else{pair = std::make_pair(id2, id1);}
    inds.insert(id1);
    inds.insert(id2);
    if (allsegs.find(pair) == allsegs.end()){
      allsegs.insert(std::make_pair(pair, new Pair()));
    }

    Pair *p = allsegs[pair];
    if (ibd12 == "IBD1"){
      p->ibd1_tot += segLen;
      if (p->ibd1_map->find(chr) == p->ibd1_map->end()){
        (*(p->ibd1_map)).insert(std::make_pair(chr, new std::vector<std::pair<double, double>>()));
      }
      (*(p->ibd1_map))[chr]->push_back(std::make_pair(start, end));
    }
    else{
      p->ibd2_tot += segLen;
      if (p->ibd2_map->find(chr) == p->ibd2_map->end()){
        (*(p->ibd2_map)).insert(std::make_pair(chr, new std::vector<std::pair<double, double>>()));
      }
      (*(p->ibd2_map))[chr]->push_back(std::make_pair(start, end));
    }
  }

  logFile.printf("\tFinished reading segments for %d samples\n", inds.size());
  return inds.size();
}

double calc_bkg_sharing(const std::string &NeFile, const Eigen::VectorXd &chrLens, const double &minIBD){
  // read the Ne File first; Don't know G yet, so need to store in a vector first
  std::vector<double> Ne_;
  // open Ne file:
  FILE *in = fopen(NeFile.c_str(), "r");
  if (!in) {
    printf("ERROR: could not open def file %s!\n", NeFile.c_str());
    exit(1);
  }

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  const char *delim = " \t\n";

  int line = 0;
  int prev_gen = -1;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    line++;
    char *token, *saveptr;
    token = strtok_r(buffer, delim, &saveptr);
    if (token == NULL || token[0] == '#') {
      // blank line or comment -- skip
      continue;
    }else{
      // parse generation number and Ne at that generation
      int gen = std::stoi(token);
      if (prev_gen != -1 && gen != prev_gen +1){
        fprintf(stderr, "ERROR: We expect NeFile's generation to increase one by one on each line, starting from the present. Line %d\n", line);
        exit(1);
      }
      prev_gen = gen;
      char *ne = strtok_r(NULL, delim, &saveptr);
      if (ne == NULL){
        fprintf(stderr, "ERROR: the NeFile expects two fields per line. Line %d\n", line);
        exit(1);
      }
      double ne_d = std::stod(ne);
      Ne_.push_back(ne_d);
    }
  }
  fclose(in);

  int maxGen = Ne_.size();
  Eigen::VectorXd Ne(maxGen);
  for(int i = 0; i < maxGen; i++){Ne[i] = Ne_[i];}
  Eigen::RowVectorXd gen = Eigen::VectorXd::LinSpaced(maxGen, 1, maxGen);
  Eigen::MatrixXd tmp1 = chrLens*gen;
  tmp1 = minIBD*tmp1/50.0;
  tmp1.colwise() += chrLens;
  tmp1.rowwise() -= minIBD*minIBD*gen/50.0;
  Eigen::RowVectorXd log_term3 = tmp1.colwise().sum().array().log();

  // Calculate expected background sharing
  Eigen::VectorXd aux(maxGen+1);
  aux(0) = 0;
  aux.tail(maxGen) = log((1 - 1.0/(2.0*Ne.array())));
  Eigen::VectorXd sum_log_prob_not_coalesce(maxGen+1);
  cumsum_eigen_colvector(aux, sum_log_prob_not_coalesce);
  Eigen::VectorXd log_expectation(maxGen+1);
  Eigen::VectorXd gen_c = Eigen::VectorXd::LinSpaced(maxGen, 1, maxGen);
  log_expectation.head(maxGen) = log(4) + sum_log_prob_not_coalesce.head(maxGen).array() - (2.0*Ne).array().log()
        - minIBD*gen_c.array()/50.0 + log_term3.transpose().array();
  double log_IBD_beyond_maxGen = log_expectedIBD_beyond_maxGen_given_Ne(Ne, chrLens, 4, minIBD);
  log_expectation(maxGen) = log_IBD_beyond_maxGen;
  return log_expectation.array().exp().sum();
}

void cumsum_eigen_colvector(const Eigen::VectorXd &source, Eigen::VectorXd &dest){
    double cum = 0;
    int n = source.rows();
    for(int i = 0; i < n; i++){
        cum += source(i);
        dest(i) = cum;
    }
}

double log_expectedIBD_beyond_maxGen_given_Ne(const Eigen::VectorXd &Ne, 
  const Eigen::VectorXd &chrLens, int n_p, double minIBD){
    int G = Ne.rows();
    double totg = chrLens.sum();
    int num_chr = chrLens.rows();
    double N_past = Ne(G-1);
    double alpha = log(1.0 - 1.0/(2.0*N_past)) - minIBD/50.0;
    double log_part_A = log(totg) + (G+1.0)*log((2.0*N_past)/(2.0*N_past-1)) + 
            alpha*(G+1) - log(1.0-exp(alpha));
    double D = (1.0*minIBD/50.0)*totg - 1.0*minIBD*minIBD*num_chr/50.0;
    double log_part_B = log(D) + (G+1.0)*log((2.0*N_past)/(2.0*N_past-1.0)) + alpha*(G+1.0) +
            log(1.0+G*(1.0-exp(alpha))) - 2.0*log(1.0-exp(alpha));
    double sum_not_coalesce = (1 - 1.0/(2.0*Ne.array())).log().sum();
    return log(n_p) + sum_not_coalesce - log(2.0*N_past) + logaddexp(log_part_A, log_part_B);

}


// specialization for FILE I/O wrapper
// open <filename> using standard FILE *
template<>
bool FileOrGZ<FILE *>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = fopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

// open <filename> as a gzipped file
template<>
bool FileOrGZ<gzFile>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = gzopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

template<>
int FileOrGZ<FILE *>::getline() {
  return ::getline(&buf, &buf_size, fp);
}

template<>
int FileOrGZ<gzFile>::getline() {
  int n_read = 0;
  int c;

  while ((c = gzgetc(fp)) != EOF) {
    // About to have read one more, so n_read + 1 needs to be less than *n.
    // Note that we use >= not > since we need one more space for '\0'
    if (n_read + 1 >= (int) buf_size) {
      const size_t GROW = 1024;
      char *tmp_buf = (char *) realloc(buf, buf_size + GROW);
      if (tmp_buf == NULL) {
	fprintf(stderr, "ERROR: out of memory!\n");
	exit(1);
      }
      buf_size += GROW;
      buf = tmp_buf;
    }
    buf[n_read] = (char) c;
    n_read++;
    if (c == '\n')
      break;
  }

  if (c == EOF && n_read == 0)
    return -1;

  buf[n_read] = '\0';

  return n_read;
}

template<>
int FileOrGZ<FILE *>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  ret = vfprintf(fp, format, args);

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<gzFile>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  // NOTE: one can get automatic parallelization (in a second thread) for
  // gzipped output by opening a pipe to gzip (or bgzip). For example:
  //FILE *pipe = popen("gzip > output.vcf.gz", "w");
  // Can then fprintf(pipe, ...) as if it were a normal file.

  // gzvprintf() is slower than the code below that buffers the output.
  // Saw 13.4% speedup for processing a truncated VCF with ~50k lines and
  // 8955 samples.
  //  ret = gzvprintf(fp, format, args);
  ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  if (ret < 0) {
    printf("ERROR: could not print\n");
    perror("printf");
    exit(10);
  }

  if (buf_len + ret > buf_size - 1) {
    // didn't fit the text in buf
    // first print what was in buf before the vsnprintf() call:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
    // now ensure that redoing vsnprintf() will fit in buf:
    if ((size_t) ret > buf_size - 1) {
      do { // find the buffer size that fits the last vsnprintf() call
	buf_size += INIT_SIZE;
      } while ((size_t) ret > buf_size - 1);
      free(buf);
      buf = (char *) malloc(buf_size);
      if (buf == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
    }
    // redo:
    ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  }

  buf_len += ret;
  if (buf_len >= buf_size - 1024) { // within a tolerance of MAX_BUF?
    // flush:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
  }

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<FILE *>::close() {
  assert(buf_len == 0);
  // should free buf, but I know the program is about to end, so won't
  return fclose(fp);
}

template<>
int FileOrGZ<gzFile>::close() {
  if (buf_len > 0)
    gzwrite(fp, buf, buf_len);
  // should free buf, but I know the program is about to end, so won't
  return gzclose(fp);
}