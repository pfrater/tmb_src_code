#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(len);
  DATA_VECTOR(age);
  PARAMETER(L1);
  PARAMETER(L2);
  PARAMETER(k);
  PARAMETER(sig_L1);
  PARAMETER(sig_L2);
  PARAMETER(sig_k);
  
  // compute mean LAA
  Type t1 = min(age);
  Type t2 = max(age);
  Type denominator = Type(1) - exp(-k * (t2 - t1));
  int n = len.size();

  vector<Type> pred_len(n);
  vector<Type> numerator(n);

  numerator = Type(1) - exp(-k * (age - t1));
  pred_len = L1 + ((L2 - L1) * (numerator / denominator));
  
  // compute sigma at age
  Type sig_denom = Type(1) - exp(-sig_k * (t2 - t1));
  vector<Type> sig_at_age(n);
  vector<Type> sig_numerator(n);
  
  sig_numerator = Type(1) - exp(-sig_k * (age - t1));
  sig_at_age = sig_L1 + ((sig_L2 - sig_L1) * (sig_numerator / sig_denom));
  
  REPORT(len);
  REPORT(age);
  REPORT(pred_len);
  REPORT(sig_at_age);
  REPORT(L1);
  REPORT(L2);
  REPORT(k);
  REPORT(sig_L1);
  REPORT(sig_L2);
  REPORT(sig_k);
  REPORT(numerator);
  REPORT(denominator);
  REPORT(sig_numerator);
  REPORT(sig_denom);
  
  // compute and return likelihood
  Type nll;
  nll = -sum(dnorm(len, pred_len, sig_at_age, true));
  REPORT(nll);
  return nll;
}
