#define TMB_LIB_INIT R_init_VBLB
#include <TMB.hpp>

using namespace density;
using namespace Eigen;
using namespace tmbutils;


template<class Type>
Type objective_function<Type>::operator() () {
  // Individual data
  DATA_IVECTOR(A);		// Vector of measured ages (sans offset)
  DATA_VECTOR(L);		// Vector of measured lengths
  DATA_VECTOR(S);		// Vector of estimated selectivity
  DATA_IVECTOR(bin);            // Vector of bin membership
  DATA_VECTOR(Aoff);            // Ageing offset

  // Length bin description
  DATA_VECTOR(Lbrk);            // Length bin boundaries
  DATA_VECTOR(Sbin);            // Approximate selectivity for each bin
  DATA_VECTOR(Fbin);            // Sampling fraction for each bin

  // Ageing
  DATA_MATRIX(Aerr);            // Ageing error matrix (probabilities)
  DATA_INTEGER(Amin);           // Minimum age

  DATA_INTEGER(model);          // Model indicator
  PARAMETER_VECTOR(logp1);      // Log constrained parameters
  PARAMETER_VECTOR(p2);         // Unconstrained parameters
  PARAMETER(logCV)		// Log CV

  // Backtransform constrained parameters
  vector<Type> p1 = exp(logp1);
  Type CV = exp(logCV);

  // Negative log likelihood
  Type nll = Type(0.0);

  // Loop over individuals
  for(int i=0; i < A.size(); i++) {
    Type lik = Type(0.0);

    // Sum over the true ages
    for(int a=Amin; a<Amin+Aerr.rows(); a++) {
      Type num = Type(0.0);
      Type den = Type(0.0);

      // Probability of the measured age for this true age
      Type prAge = Aerr(a-Amin,A(i)-Amin);

      if(prAge > 0.0) {
	Type age = a + Aoff(i);

	// Mean length for this true age
	Type mu = 0.0;
	switch(model) {
	case 1: // von Bertalanffy growth
	  mu = p1(0)*(Type(1.0)-exp(-p1(1)*(age-p2(0))));
	  break;
	case 2: // Gompertz
	  mu = p1(0)*exp(-exp(-p1(1)*(age-p2(0)))/p1(1));
	  break;
	case 3: // Schnute-Richards
	  mu = p1(0)*(Type(1.0)+p1(3)*exp(-p1(2)*(age^p1(3))))^(Type(-1.0)/p1(4))
	  break;
	case 4: // logistic
	  mu = p1(0)*(Type(1.0)+exp(-p1(1)*(age-p2(0))))^Type(-1.0)
	  break;
	}
	Type sigma = CV*mu;

	// Probability of the measured length adjusted for selectivity
	// and sampling fraction for this true age
	// Pr(Length) * Pr(Sel) * Pr(aged)
	num = dnorm(L(i),mu,sigma)*S(i)*Fbin(bin(i)-1);

	// Compute the total probability of being selected and sampled for this
	// true age by summing over all length bins
	den += (pnorm(Lbrk(0),mu,sigma)-pnorm(Type(0.0),mu,sigma))*Sbin(0)*Fbin(0);
	for(int b=1; b<Lbrk.size(); b++)
	  den += (pnorm(Lbrk(b),mu,sigma)-pnorm(Lbrk(b-1),mu,sigma))*Sbin(b)*Fbin(b);
	den += (Type(1.0)-pnorm(Lbrk(Lbrk.size()-1),mu,sigma))*Sbin(Lbrk.size())*Fbin(Lbrk.size());


	// Contribution to marginal likelihood
	lik += prAge*num/den;;
      }
    }
    // Contribution to neg log likelihood
    nll -= log(lik);
  }

  // Report transformed parameters
  ADREPORT(p1);
  ADREPORT(p2);
  ADREPORT(CV);

  //Return negative log likelihood
  return nll;
}
