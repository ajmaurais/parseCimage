
#include <Rcpp.h>
#include <vector>

template <typename _Tp>
double median(std::vector<_Tp>& v){
  size_t size = v.size();
  
  if (size == 0){
    return 0;  // Undefined, really.
  }
  else if (size == 1){
    return v[0];
  }
  else
  {
    std::sort(v.begin(), v.end());
    if (size % 2 == 0){
      return (v[size / 2 - 1] + v[size / 2]) / 2;
    }
    else 
    {
      return v[size / 2];
    }
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector normRat(Rcpp::NumericVector vec, double mult = 1) {
  std::vector<double> tempVec;
  size_t len = vec.size();
  for(int i = 0; i < len; i++){
    if(vec[i] != 0)
      tempVec.push_back(vec[i]);
  }
  double med = median(tempVec);
  
  for(int i = 0; i < len; i++){
    if(vec[i] != 20)
      vec[i] = (vec[i] / med) * mult;
  }
  return vec;
}


