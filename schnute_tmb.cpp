#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(len);
  DATA_VECTOR(age);
  PARAMETER(L1);
  PARAMETER(L2);
  PARAMETER(k);
  PARAMETER(log_sigma);
  
  // compute mean LAA
  Type t1 = min(age);
  Type t2 = max(age);
  Type denominator = Type(1) - exp(-k * (t2 - t1));
  int n = age.size();

  vector<Type> pred_length(n);
  vector<Type> numerator(n);

  numerator = Type(1) - exp(-k * (age - t1));
  pred_length = L1 + ((L2 - L1) * (numerator / denominator));
  
  // calculate negative log-likelihood
  Type nll = -sum(dnorm(len, pred_length, exp(log_sigma), true));
  return nll;
}
