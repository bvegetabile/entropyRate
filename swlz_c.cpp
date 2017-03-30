#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame swlz_er(std::string x, int start_pos = 1) {
  // int start_pos = 1;
  
  int str_len = x.length();
  int trials = str_len - start_pos;
  // Rprintf("Input String Length: %d\n", str_len);
  IntegerVector i_vec(trials);
  IntegerVector MatchLength(trials);
  NumericVector InstaEntropy(trials);
  StringVector NewString(trials);
  LogicalVector Include(trials);
  
  int left_pointer = start_pos;
  int right_pointer = start_pos;
  int max_pointer = str_len - 1;
  int trial = 0;
  
  std::string left_word;
  std::string right_word;
  
  while(right_pointer < max_pointer){
    left_word = x.substr(0, left_pointer);
    while(right_pointer <= max_pointer){
      int nc = right_pointer - left_pointer + 1;
      right_word = x.substr(left_pointer, nc);
      if(left_word.find(right_word) == -1){
        i_vec(trial) = left_pointer;
        MatchLength(trial) = nc;
        NewString(trial) = right_word;
        Include(trial) = TRUE;
        trial += 1;
        left_pointer += 1;
        right_pointer = left_pointer;
        break;
      } else {
        right_pointer += 1;
      }
    }
  }
  return DataFrame::create(_["i"] = i_vec, 
                           _["MatchLength"] = MatchLength, 
                           _["NewString"] = NewString, 
                           _["Include"] = Include);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

# /*** R
# swlz_er("13131213232331313332")
# */
