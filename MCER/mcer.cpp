// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
# include <Rcpp.h>

using namespace Rcpp;

int find_pos(CharacterVector v, std::string s){
  CharacterVector::iterator it;
  int pos = 0;
  for(it = v.begin(); it < v.end(); it++){
    if(*it == s){
      return pos;
    } else {
      pos += 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
double mcer(std::string mc_seq,
               int mc_order,
               CharacterVector uniq_states,
               CharacterVector state_space) {
  arma::mat tc(state_space.length(),
               uniq_states.length(),
               arma::fill::zeros);
  arma::mat tm(state_space.length(),
               uniq_states.length(),
               arma::fill::zeros);
  arma::vec sm(state_space.length(),
               arma::fill::zeros);

  int mc_len = mc_seq.length();
  std::string x_minus;
  std::string x_plus;
  int r_pos;
  int c_pos;

  for(int i = 0; i < mc_len - mc_order; i++){
    x_minus = mc_seq.substr(i, mc_order);
    x_plus = mc_seq.substr(i+mc_order, 1);
    r_pos = find_pos(state_space, x_minus);
    c_pos = find_pos(uniq_states, x_plus);
    tc(r_pos, c_pos) += 1;
    sm(r_pos) += 1;
  }

  sm = sm / arma::accu(sm);

  arma::mat rs(state_space.length(),
               uniq_states.length(),
               arma::fill::zeros);
  arma::mat sm_mat(state_space.length(),
                   uniq_states.length(),
                   arma::fill::zeros);

  rs.each_col() = arma::sum(tc, 1);
  sm_mat.each_col() = sm;

  tm = tc / rs;

  arma::mat res = - tm % log2(tm) % sm_mat;

  res.replace(arma::datum::nan, 0);

  return arma::accu(res);
}


// [[Rcpp::export]]
double lz77entropy(std::string x, int window_size = 10){
  int str_len = x.length();
  int start_pos = window_size;
  int right_pointer = start_pos;
  int pos = start_pos;
  int max_pointer = str_len - 1;
  int n_char;
  int counter = 0;
  double L = 0;
  double ent = 0;

  std::string window;
  std::string right_word;

  while(right_pointer < max_pointer){
    window = x.substr(pos-window_size, window_size);
    while(right_pointer <= max_pointer){
      n_char = right_pointer - pos + 1;
      right_word = x.substr(pos, n_char);
      if(window.find(right_word) == -1){
        counter += 1;
        L += n_char;
        ent += std::log2(window_size) / n_char;
        pos += 1;
        right_pointer = pos;
        // std::cout << right_word << '\n';
        break;
      } else {
        right_pointer += 1;
      }
    }
  }




  return ent / counter;
}


