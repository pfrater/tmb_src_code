// code for this file not complete

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(length);
  DATA_VECTOR(age);
  PARAMETER(L1);
  PARAMETER(L2);
  PARAMETER(L3);
  PARAMETER(t1);
  PARAMETER(t2);
  PARAMETER(t3)
  PARAMETER(log_sigma);
  int n = length.size();
  int r = (L3 - L2) / (L2 - L1);

  vector<Type> pred_length(n);
  vector<Type> numerator(n);
  vector<Type> denominator(n);
  for (int i = 0; i < n; i++) {
    numerator(i) = (1 - (r^2)^((age(i) - t1) / (t3 - t1)))
    denominator(i) = 1 - (r^2)
    pred_length(i) = L1 + ((L3 - L1) * (numerator / denominator));
  }
  
  Type nll = -sum(dnorm(length, pred_length, exp(log_sigma), true));
  return nll;
}
