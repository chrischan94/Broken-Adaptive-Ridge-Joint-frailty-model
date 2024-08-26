// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

arma::rowvec smhrCpp(NumericVector y, double yi) {

	/*(Constant piecewise) this function finds the quintile which y is in*/
	NumericVector Yr(y);
	NumericVector Pr = { 0,0.2,0.4,0.6,0.8,1 }; /* don't manually create the numeric vector. Put this as an argument */

	arma::rowvec Y(Yr.begin(), Yr.size(), false);
	arma::rowvec P(Pr.begin(), Pr.size(), false);
	arma::rowvec Ys = arma::quantile(Y, P);

	Ys[Ys.size() - 1] = Ys[Ys.size() - 1] + 1e-04; /* Add small number to prevent calculation error */
	int l = Ys.size() - 1;
	arma::rowvec h(l);

	for (int j = 0; j < l; j++)
		if (yi >= Ys[j] && yi < Ys[j + 1])
			h[j] = 1;

	return h;
}

//[[Rcpp::export]]
arma::rowvec smhrCpp2(NumericVector y, double yi) {
  
  /*(Constant piecewise) this function finds the quintile which y is in*/
  NumericVector Yr(y);
  NumericVector Pr = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 }; /* don't manually create the numeric vector. Put this as an argument */
  
  arma::rowvec Y(Yr.begin(), Yr.size(), false);
  arma::rowvec P(Pr.begin(), Pr.size(), false);
  arma::rowvec Ys = arma::quantile(Y, P);
  
  Ys[Ys.size() - 1] = Ys[Ys.size() - 1] + 1e-04; /* Add small number to prevent calculation error */
  int l = Ys.size() - 1;
  arma::rowvec h(l);
  
  for (int j = 0; j < l; j++) {
	  if (yi >= Ys[j] && yi < Ys[j + 1])
		  h[j] = 1;
  }

    return h;
}

arma::rowvec caphrCpp(NumericVector y, double yi) { /* Change name of function to avoid confusion*/
	/*H_0(t) or R_0(t)*/

	NumericVector Yr(y);
	NumericVector Pr = { 0,0.2,0.4,0.6,0.8,1 };

	arma::rowvec Y(Yr.begin(), Yr.size(), false);
	arma::rowvec P(Pr.begin(), Pr.size(), false);
	arma::rowvec Ys = arma::quantile(Y, P);

	Ys[0] = Ys[0] - 1e-06;
	int l = Ys.size() - 1;
	arma::rowvec H(l);

	for (int j = 0; j < l; j++) {
		double A = Ys[j + 1] - Ys[j];
		double B = yi - Ys[j];
		double M = std::min((float)A, (float)B);
		H[j] = std::max((float)M, (float)0);
	}
	// change double to float to avoid numerical instability. More investigation is needed. 
	return H;
}

//[[Rcpp::export]]
arma::rowvec caphrCpp2(NumericVector y, double yi) { /* Change name of function to avoid confusion*/
	/*H_0(t) or R_0(t)*/
	
  NumericVector Yr(y);
  NumericVector Pr = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 };
  
  arma::rowvec Y(Yr.begin(), Yr.size(), false);
  arma::rowvec P(Pr.begin(), Pr.size(), false);
  arma::rowvec Ys = arma::quantile(Y, P);
  
  Ys[0] = Ys[0] - 1e-06;
  int l = Ys.size() - 1;
  arma::rowvec H(l);
  
  for (int j = 0; j < l; j++) {
    double A = Ys[j + 1] - Ys[j];
    double B = yi - Ys[j];
    double M = std::min((float)A, (float)B);
    H[j] = std::max((float)M, (float)0);
  }
  // change double to float to avoid numerical instability. More investigation is needed. 
  return H;
}

NumericVector findCpp2(NumericMatrix Z, int i) {
	/*This function returns the sequence of recurrent event times for id i.*/

	arma::rowvec id = Z(_, 1);
	arma::rowvec Tij = Z(_, 0);

	arma::uvec pos = arma::find(id == i);
	arma::vec Tij_pos = Tij.elem(pos);

	return wrap(Tij_pos);
}

double dnormCpp(double x, double m, double v) {

	return (R::dnorm(x, m, v, FALSE));

}

NumericVector Concvecelem(arma::rowvec a, arma::rowvec b, arma::rowvec c, arma::rowvec d, double e, double f) {

	int l1 = a.size();
	int l2 = b.size();
	int l3 = c.size();
	int l4 = d.size();

	arma::rowvec ans(l1 + l2 + l3 + l4 + 2, arma::fill::zeros);

	ans(arma::span(0, l1 + l2 + l3 + l4 - 1)) = arma::join_rows(a, b, c, d);
	ans[l1 + l2 + l3 + l4] = e;
	ans[l1 + l2 + l3 + l4 + 1] = f;

	return wrap(ans);

}

//[[Rcpp::export]]
double gCpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl,
	NumericMatrix z, NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericMatrix Zr(z);
	NumericVector theta(the);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}

	return ans;
}

//[[Rcpp::export]]
double g2Cpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl,
	NumericMatrix z, NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericMatrix Zr(z);
	NumericVector theta(the);
	NumericVector delta(dl);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * dnormCpp(u, mean, phi);
	}

	return ans;
}

//[[Rcpp::export]]
double g3Cpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl, NumericMatrix z,
	NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericMatrix Zr(z);
	NumericVector theta(the);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}

	return ans;
}

//[[Rcpp::export]]
double g4Cpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl,
	NumericMatrix z, NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericMatrix Zr(z);
	NumericVector theta(the);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (pow(ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u), 2) - (sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u))) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (pow(ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u), 2) - (sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u))) * dnormCpp(u, mean, phi);
	}

	return ans;
}

//[[Rcpp::export]]
double g5Cpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl, NumericMatrix z,
	NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector ni(r);
	NumericVector Tij = R(_, 0);
	NumericVector delta(dl);
	NumericVector theta(the);
	NumericMatrix Zr(z);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (pow(delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma*u), 2) - (sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma*u))) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (pow(delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma*u), 2) - (sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma*u))) * dnormCpp(u, mean, phi);
	}

	return ans;
}

//[[Rcpp::export]]
double g6Cpp(double u, NumericVector the, NumericVector y, NumericVector r, NumericVector dl, NumericMatrix z,
	NumericMatrix R, int j, int l1, int l2, int l3, int l4) {

	NumericVector ni(r);
	NumericVector Tij = R(_, 0);
	NumericVector delta(dl);
	NumericVector theta(the);
	NumericMatrix Zr(z);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	arma::rowvec Z = Zr(j, _);
	double mean = 0;
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * 
			(delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}
	else {
		ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) * pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[j])), delta[j]) *
			exp(-sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * (ni[j] - sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u)) *
			(delta[j] - sum(hs % caphrCpp(y, y[j])) * exp(sum(beta2 % Z) + gamma * u)) * dnormCpp(u, mean, phi);
	}
	return ans;
}

//[[Rcpp::export]]
arma::mat sigmaHCpp(NumericVector the, NumericVector y, NumericMatrix z, NumericVector r, NumericMatrix R,
	NumericVector dl, arma::mat muHat, int l1, int l2, int l3, int l4) {
	/*finds the standard deviations */

	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();
	arma::mat Ans(n, 6, arma::fill::zeros);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		Ans(i, 0) = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i,0))) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i,0))) - pow(phi, -2);
		Ans(i, 1) = Ans(i, 0) -
			(ni[i] * sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))) /
			pow((ni[i] - sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))), 2);
		Ans(i, 2) = Ans(i, 0) -
			(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))) /
			pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))), 2);
		Ans(i, 3) = Ans(i, 0) +
			(4 * pow(ni[i], 2) * pow((sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) -
				(pow((sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 3) * (2 * ni[i] + 1)) -
				(pow(ni[i], 2) * (sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))) * (1 + 2 * ni[i]))) /
			(pow(pow((ni[i] - sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) - 
				sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3)), 2));
		Ans(i, 4) = Ans(i, 0) + 
			(4 * pow(gamma, 2) * delta[i] * pow(sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2) -
			pow(gamma, 2) * pow(sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 3) * (1 + 2 * delta[i]) -
			3 * pow(gamma, 2) * delta[i] * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))) /
			(pow(pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))), 2) -
				sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2));
		Ans(i, 5) = Ans(i, 0) -
			(ni[i] * sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 5))) /
			pow((ni[i] - sum(rs % caphrCpp(y, y[i]) * exp(sum(beta1 % Z) + muHat(i, 5)))), 2) -
			(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))) /
			pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))), 2);
	}

	return Ans;
}

//[[Rcpp::export]]
arma::mat sigmaHCpp2(NumericVector the, NumericVector y, NumericMatrix z, NumericVector r, NumericMatrix R, NumericVector dl, arma::mat muHat, int l1, int l2, int l3, int l4) {
/*finds the standard deviations */

  NumericVector theta(the);
  NumericMatrix Zr(z);
  NumericVector ni(r);
  NumericVector delta(dl);

  arma::rowvec beta1(theta.begin(), l1);
  arma::rowvec beta2(theta.begin() + l1, l2);
  arma::rowvec hs(theta.begin() + l1 + l2, l3);
  arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
  double gamma = theta[l1 + l2 + l3 + l4];
  double phi = theta[l1 + l2 + l3 + l4 + 1];

  int n = delta.size();
  arma::mat Ans(n, 6, arma::fill::zeros);

  for (int i = 0; i < n; i++) {

	arma::rowvec Z = Zr(i, _);

	Ans(i, 0) = -pow(gamma, 2) * (sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 0))) -
		(sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 0))) - pow(phi, -2);
	Ans(i, 1) = Ans(i, 0) -
		(ni[i] * sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))) /
		pow((ni[i] - sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))), 2);
	Ans(i, 2) = Ans(i, 0) -
		(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))) /
		pow((delta[i] - sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))), 2);
	Ans(i, 3) = Ans(i, 0) +
		(4 * pow(ni[i], 2) * pow((sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) -
			(pow((sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 3) * (2 * ni[i] + 1)) -
			(pow(ni[i], 2) * (sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))) * (1 + 2 * ni[i]))) /
		(pow(pow((ni[i] - sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) -
			sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3)), 2));
	Ans(i, 4) = Ans(i, 0) +
		(4 * pow(gamma, 2) * delta[i] * pow(sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2) -
			pow(gamma, 2) * pow(sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 3) * (1 + 2 * delta[i]) -
			3 * pow(gamma, 2) * delta[i] * sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))) /
		(pow(pow((delta[i] - sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))), 2) -
			sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2));
	Ans(i, 5) = Ans(i, 0) -
		(ni[i] * sum(rs % caphrCpp2(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 5))) /
		pow((ni[i] - sum(rs % caphrCpp2(y, y[i]) * exp(sum(beta1 % Z) + muHat(i, 5)))), 2) -
		(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))) /
		pow((delta[i] - sum(hs % caphrCpp2(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))), 2);
   }

   return Ans;
}

//[[Rcpp::export]]
NumericVector BARCppwAGQ(NumericVector the, NumericVector y, NumericVector r, NumericVector dl,
	NumericMatrix z, NumericMatrix R, NumericMatrix muHat, arma::mat sdH, arma::rowvec x, arma::rowvec w, double lam, int l1, int l2, int l3, int l4) {

	NumericVector Tij = R(_, 0);
	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();
	double perb = 1e-06;
	arma::mat beta1prime(l1, n, arma::fill::zeros); /* first derivative of beta1 or beta_r */
	arma::mat beta2prime(l2, n, arma::fill::zeros); /* first derivative of beta2 or beta_h */
	arma::cube beta1doubprime(l1, l1, n, arma::fill::zeros); /* */
	arma::cube beta2doubprime(l2, l2, n, arma::fill::zeros); /* */
	arma::cube betacross(l1, l2, n, arma::fill::zeros); /*cross product of derivatives*/

	arma::rowvec likeAGQ(n);
	arma::rowvec g4(n);
	arma::rowvec g5(n);
	arma::rowvec g6(n);
	arma::rowvec g7(n);
	arma::rowvec g8(n);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);
		arma::colvec Zt = Z.t();

		arma::rowvec s(10);
		arma::rowvec s2(10);
		arma::rowvec s3(10);
		arma::rowvec s4(10);
		arma::rowvec s5(10);
		arma::rowvec s6(10);

		for (int j = 0; j < 10; j++) {
			s[j] = w[j] * exp(pow(x[j], 2)) * gCpp(muHat(i, 0) + pow(2, 0.5) * pow(-sdH(i, 0), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
			s2[j] = w[j] * exp(pow(x[j], 2)) * g2Cpp(muHat(i, 1) + pow(2, 0.5) * pow(-sdH(i, 1), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
			s3[j] = w[j] * exp(pow(x[j], 2)) * g3Cpp(muHat(i, 2) + pow(2, 0.5) * pow(-sdH(i, 2), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
			s4[j] = w[j] * exp(pow(x[j], 2)) * g4Cpp(muHat(i, 3) + pow(2, 0.5) * pow(-sdH(i, 3), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
			s5[j] = w[j] * exp(pow(x[j], 2)) * g5Cpp(muHat(i, 4) + pow(2, 0.5) * pow(-sdH(i, 4), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
			s6[j] = w[j] * exp(pow(x[j], 2)) * g6Cpp(muHat(i, 5) + pow(2, 0.5) * pow(-sdH(i, 5), -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
		}

		likeAGQ[i] = pow(2, 0.5) * pow(-sdH(i, 0), -0.5) * sum(s);
		g4[i] = pow(2, 0.5) * pow(-sdH(i, 1), -0.5) * sum(s);
		g5[i] = pow(2, 0.5) * pow(-sdH(i, 2), -0.5) * sum(s);
		g6[i] = pow(2, 0.5) * pow(-sdH(i, 3), -0.5) * sum(s);
		g7[i] = pow(2, 0.5) * pow(-sdH(i, 4), -0.5) * sum(s);
		g8[i] = pow(2, 0.5) * pow(-sdH(i, 5), -0.5) * sum(s);

		beta1prime.col(i) = Zt * (g4[i] / likeAGQ[i]);
		beta2prime.col(i) = Zt * (g5[i] / likeAGQ[i]);
		beta1doubprime.slice(i) = Zt * Z * ((g6[i] * likeAGQ[i] - pow(g4[i], 2)) / pow(likeAGQ[i], 2));
		beta2doubprime.slice(i) = Zt * Z * ((g7[i] * likeAGQ[i] - pow(g5[i], 2)) / pow(likeAGQ[i], 2));
		betacross.slice(i) = Zt * Z * ((g8[i] * likeAGQ[i] - g4[i] * g5[i]) / pow(likeAGQ[i], 2));
	}
	
	arma::rowvec fd1(l1);
	arma::rowvec fd2(l2);

	for (int k = 0; k < l1; k++) {
		for (int m = 0; m < n; m++) {
			fd1[k] += beta1prime(k, m);
			fd2[k] += beta2prime(k, m);
		}
	}

	arma::rowvec fd = arma::join_rows(fd1, fd2); /* Gradient vector */
	arma::mat sd1 = arma::sum(beta1doubprime, 2);
	arma::mat sd2 = arma::sum(betacross, 2);
	arma::mat sd4 = arma::sum(beta2doubprime, 2);
	arma::mat sd = arma::join_rows(arma::join_cols(sd1, sd2), arma::join_cols(sd2, sd4)); /* Hessian matrix */
    arma::mat Omega_n = -sd; /*Omega_n*/
	arma::rowvec betas = arma::join_rows(beta1, beta2);
	arma::rowvec betas_sq(l1 + l2);
	for (int a = 0; a < (l1 + l2); a++) {
		betas_sq[a] = pow((pow(betas[a], 2) + perb), -1);
	}
	arma::mat Db = arma::diagmat(betas_sq);
	arma::colvec v_n = fd.t() - sd * betas.t(); /*v_n*/
	arma::colvec betaBAR = arma::inv(Omega_n + lam * Db) * v_n;
	arma::colvec betar_NEW = betaBAR(arma::span(0, l1 - 1));
	arma::colvec betah_NEW = betaBAR(arma::span(l1, 2 * l2 - 1));
	NumericVector answer = Concvecelem(betar_NEW.t(), betah_NEW.t(), hs, rs, gamma, phi);
	
	return wrap(Omega_n);
}

//[[Rcpp::export]]
arma::mat sigmaHwGammaCpp(NumericVector the, NumericVector y, NumericMatrix z, NumericVector r, NumericMatrix R,
	NumericVector dl, arma::mat muHat, int l1, int l2, int l3, int l4) {
	/*finds the standard deviations */

	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();
	arma::mat Ans(n, 6, arma::fill::zeros);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		Ans(i, 0) = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 0))) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 0))) - phi * exp(muHat(i, 0));
		Ans(i, 1) = Ans(i, 0) -
			(ni[i] * sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))) /
			pow((ni[i] - sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 1))), 2);
		Ans(i, 2) = Ans(i, 0) -
			(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))) /
			pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 2))), 2);
		Ans(i, 3) = Ans(i, 0) +
			(4 * pow(ni[i], 2) * pow((sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) -
				(pow((sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 3) * (2 * ni[i] + 1)) -
				(pow(ni[i], 2) * (sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))) * (1 + 2 * ni[i]))) /
			(pow(pow((ni[i] - sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3))), 2) -
				sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 3)), 2));
		Ans(i, 4) = Ans(i, 0) +
			(4 * pow(gamma, 2) * delta[i] * pow(sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2) -
				pow(gamma, 2) * pow(sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 3) * (1 + 2 * delta[i]) -
				3 * pow(gamma, 2) * delta[i] * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))) /
			(pow(pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4))), 2) -
				sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 4)), 2));
		Ans(i, 5) = Ans(i, 0) -
			(ni[i] * sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + muHat(i, 5))) /
			pow((ni[i] - sum(rs % caphrCpp(y, y[i]) * exp(sum(beta1 % Z) + muHat(i, 5)))), 2) -
			(delta[i] * pow(gamma, 2) * sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))) /
			pow((delta[i] - sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * muHat(i, 5))), 2);
	}

	return Ans;
}
