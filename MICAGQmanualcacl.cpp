// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec smhrCpp(NumericVector y, double yi) {

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

// [[Rcpp::export]]
arma::rowvec caphrCpp(NumericVector y, double yi) { /* Change name of function to avoid confusion*/
	NumericVector Yr(y);
	NumericVector Pr = { 0,0.2,0.4,0.6,0.8,1 };

	arma::rowvec Y(Yr.begin(), Yr.size(), false);
	arma::rowvec P(Pr.begin(), Pr.size(), false);
	arma::rowvec Ys = arma::quantile(Y, P);

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

// [[Rcpp::export]]
NumericVector findCpp2(NumericMatrix Z, int i) {

	arma::rowvec id = Z(_, 1);
	arma::rowvec Tij = Z(_, 0);

	arma::uvec pos = arma::find(id == i);
	arma::vec Tij_pos = Tij.elem(pos);

	return wrap(Tij_pos);
	//This function returns the sequence of recurrent event times for id i.
}

//[[Rcpp::export]]
double dnormCpp(double x, double m, double v) {

	return (R::dnorm(x, m, v, FALSE));

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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
double MICAGQ(NumericVector the, arma::rowvec w, arma::rowvec x, NumericVector uh, NumericVector y, NumericMatrix R,
	NumericVector r, NumericVector dl, NumericMatrix z, int l1, int l2, int l3, int l4) {

	NumericVector theta(the); /* Vector of parameters */
	NumericVector Tij = R(_, 0); /* observed recurrent event times */
	NumericVector delta(dl); /* terminal event indicator */
	NumericVector ni(r); /* # of recurrent events */
	NumericMatrix Zr(z); /* Covariate matrix */

	arma::rowvec alpha1(theta.begin(), l1);
	arma::rowvec alpha2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();
	double an = n / 4;
	arma::rowvec omega_alphar(l1);
	arma::rowvec omega_alphah(l2);

	for (int a = 0; a < l1; a++) {
		omega_alphar[a] = (exp(2 * an * pow(alpha1[a], 2)) - 1) / (exp(2 * an * pow(alpha1[a], 2)) + 1);
		omega_alphah[a] = (exp(2 * an * pow(alpha2[a], 2)) - 1) / (exp(2 * an * pow(alpha2[a], 2)) + 1);
	}

	arma::rowvec beta1 = omega_alphar % alpha1;
	arma::rowvec beta2 = omega_alphah % alpha2;
	NumericVector theta_new = Concvecelem(beta1, beta2, hs, rs, gamma, phi);

	arma::rowvec sdH(n); /*standard deviation*/
	arma::rowvec ll(n);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		sdH[i] = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * uh[i])) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + uh[i])) - pow(phi, -2);

		arma::rowvec s(10);

		for (int j = 0; j < 10; j++) {
			s[j] = w[j] * exp(pow(x[j], 2)) * gCpp(uh[i] + pow(2, 0.5) * pow(-sdH[i], -0.5) * x[j], theta_new, y, r, dl, z, R, i, l1, l2, l3, l4);
		}

		ll[i] = log(pow(2, 0.5) * pow(-sdH[i], -0.5) * sum(s));
	}

	int m = Tij.size();
	double pen = (sum(omega_alphah) + sum(omega_alphar)) * (log(m) / 2);
	double ans = -sum(ll) + pen;
	return ans;
}