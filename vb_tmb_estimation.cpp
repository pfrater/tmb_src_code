#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(len);
  DATA_VECTOR(age);
  PARAMETER(log_linf);
  PARAMETER(log_k);
  PARAMETER(t0);
  PARAMETER(log_sigma);
  int n = length.size();
  
  // backtransform parameters
  Type linf = exp(log_linf);
  Type k = exp(log_k);
  Type sigma = exp(log_sigma);

  vector<Type> pred_length(n);
  pred_length = linf * (Type(1) - exp(-k * (age - t0)));
  // for (int i = 0; i < n; i++) {
  //   pred_length(i) = linf * (1 - exp(-k * (age(i) - t0)));
  // }
  
  Type nll = -sum(dnorm(len, pred_length, sigma, true));
  return nll;
}
