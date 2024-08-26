// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector TermIndCpp(NumericMatrix Z) {

	arma::rowvec id = Z(_, Z.ncol() - 2);
	arma::rowvec delta = Z(_, Z.ncol() - 3);

	arma::uvec pos = arma::find(id == 1);
	arma::vec delta_pos = delta.elem(pos);

	return wrap(delta_pos);
	//returns the terminal event indicator
}

// [[Rcpp::export]]
NumericVector RecurIndCpp(NumericMatrix Z, int n) {

	arma::rowvec id = Z(_, Z.ncol() - 4);
	arma::rowvec count(n);
	arma::rowvec indi(n);

	for (int i = 0; i < n; i++) {
		arma::uvec pos = arma::find(id == i + 1);
		count[i] = pos.n_elem;
		if (count[i] > 1) {
			indi[i] = 1;
		}
		else {
			indi[i] = 0;
		}
	}

	return wrap(indi);
	//returns recurrent event indicator;
}

//[[Rcpp::export]]
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

//[[Rcpp::export]]
arma::rowvec caphrCpp(NumericVector y, double yi) { /* Change name of function to avoid confusion*/
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

// [[Rcpp::export]]
arma::mat RecurDataCpp(NumericMatrix Z) {

	arma::colvec recur_times = Z(_, Z.ncol() - 5);
	arma::colvec Ddelta = Z(_, Z.ncol() - 2);
	arma::colvec id = Z(_, Z.ncol() - 4);

	arma::uvec pos = arma::find(Ddelta == 0);
	arma::vec id_pos = id.elem(pos);
	arma::vec recur_times_pos = recur_times.elem(pos);

	arma::mat RecurData = arma::join_rows(recur_times_pos, id_pos);
	return RecurData;
	//returns the recurrent events data matrix
}

// [[Rcpp::export]]
NumericMatrix ZsCpp(NumericMatrix Z, int p, int n) {

	arma::colvec Ddelta = Z(_, Z.ncol() - 2);
	arma::uvec pos = arma::find(Ddelta == 1);

	arma::mat Zs(n, p, arma::fill::zeros);

	for (int i = 0; i < p; i++) {
		arma::vec Zi = Z(_, i);
		arma::vec Zi_pos = Zi.elem(pos);
		Zs.col(i) = Zi_pos;
	}

	return wrap(Zs);
	//returns the design matrix for the covariates Z
}

// [[Rcpp::export]]
NumericVector YCpp(NumericMatrix Z) {

	arma::colvec Ddelta = Z(_, Z.ncol() - 2);
	arma::colvec Y = Z(_, Z.ncol() - 6);

	arma::uvec pos = arma::find(Ddelta == 1);

	arma::vec Y_pos = Y.elem(pos);

	return (wrap(Y_pos));
	/*returns the follow-up time*/

}

// [[Rcpp::export]]
NumericVector niCpp(NumericVector x, NumericMatrix recur) {

	arma::rowvec recur_times = recur(_, 0);
	arma::rowvec id = recur(_, 1);
	int n = x.length();
	arma::rowvec ni(n);

	for (int i = 0; i < n; i++) {
		if (x[i] == 0) {
			ni[i] = 0;
		}
		else {
			arma::uvec pos = arma::find(id == i + 1);
			arma::vec recur_times_id = recur_times.elem(pos);
			ni[i] = recur_times_id.n_elem;
		}
	}

	return (wrap(ni));
	/*returns vector of the number of recurrent events */
}


NumericVector findCpp2(NumericMatrix Z, int i) {

	arma::rowvec id = Z(_, 1);
	arma::rowvec Tij = Z(_, 0);

	arma::uvec pos = arma::find(id == i);
	arma::vec Tij_pos = Tij.elem(pos);

	return wrap(Tij_pos);
	//This function returns the sequence of recurrent event times for id i.
}

double dnormCpp(double x, double m, double v) {

	return (R::dnorm(x, m, v, FALSE));

}

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
double LikeAGQ(NumericVector the, arma::rowvec w, arma::rowvec x, NumericVector uh, NumericVector y, NumericMatrix R,
	NumericVector r, NumericVector dl, NumericMatrix z, int l1, int l2, int l3, int l4) {

	NumericVector theta(the); /* Vector of parameters */
	NumericVector Tij = R(_, 0); /* observed recurrent event times */
	NumericVector delta(dl); /* terminal event indicator */
	NumericVector ni(r); /* # of recurrent events */
	NumericMatrix Zr(z); /* Covariate matrix */

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();

	arma::rowvec sdH(n); /*standard deviation*/
	arma::rowvec ll(n);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		sdH[i] = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * uh[i])) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + uh[i])) - pow(phi, -2);

		arma::rowvec s(10);

		for (int j = 0; j < 10; j++) {
			s[j] = w[j] * exp(pow(x[j], 2)) * gCpp(uh[i] + pow(2, 0.5) * pow(-sdH[i], -0.5) * x[j], theta, y, r, dl, z, R, i, l1, l2, l3, l4);
		}

		ll[i] = log(pow(2, 0.5) * pow(-sdH[i], -0.5) * sum(s));
	}

	return (-sum(ll));
}

//[[Rcpp::export]]
double prodRecurCpp(arma::rowvec rs, NumericMatrix R, int nid, int i) {

	NumericVector Tij = R(_, 0);
	NumericVector recur_times = findCpp2(R, i);
	double pr = 1;

	for (int j = 0; j < nid; j++) {
		pr *= sum(rs % smhrCpp(Tij, recur_times[j]));
	}

	return pr;
}

//[[Rcpp::export]]
double prodRecurCpp2(arma::rowvec rs, NumericMatrix R, int nid, int i) {

	NumericVector Tij = R(_, 0);
	NumericVector recur_times = findCpp2(R, i);
	double pr = 1;

	for (int j = 0; j < nid; j++) {
		pr *= sum(rs % smhrCpp2(Tij, recur_times[j]));
	}

	return pr;
}

//[[Rcpp::export]]
double elemmultsum(arma::rowvec a, arma::rowvec b) {

	return (sum(a % b));

}

//[[Rcpp::export]]
NumericVector sdHCpp(NumericVector the, NumericVector uh, NumericVector y, NumericMatrix R,
	NumericVector r, NumericVector dl, NumericMatrix z, int l1, int l2, int l3, int l4) {

	NumericVector theta(the); /* Vector of parameters */
	NumericVector Tij = R(_, 0); /* observed recurrent event times */
	NumericVector delta(dl); /* terminal event indicator */
	NumericVector ni(r); /* # of recurrent events */
	NumericMatrix Zr(z); /* Covariate matrix */

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	int n = delta.size();
	arma::rowvec sdH(n); /*standard deviation*/

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		sdH[i] = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * uh[i])) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(beta1 % Z) + uh[i])) - pow(phi, -2);
	}

	return wrap(sdH);
}

//[[Rcpp::export]]
double NuisUpd(NumericVector the, arma::rowvec w, arma::rowvec x, NumericVector uh, NumericVector y, NumericMatrix R,
	NumericVector r, NumericVector dl, NumericMatrix z, arma::rowvec br, arma::rowvec bh, int l3, int l4) {

	/*x refers to the absciccas, w refers to the weights*/
	NumericVector theta(the); /* Vector of parameters */
	NumericVector Tij = R(_, 0); /* observed recurrent event times */
	NumericVector delta(dl); /* terminal event indicator */
	NumericVector ni(r); /* # of recurrent events */
	NumericMatrix Zr(z); /* Covariate matrix */

	arma::rowvec hs(theta.begin(), l3);
	arma::rowvec rs(theta.begin() + l3, l4);
	double gamma = theta[l3 + l4];
	double phi = theta[l3 + l4 + 1];

	int n = delta.size();
	int l1 = br.size();
	int l2 = bh.size();
	arma::rowvec sdH(n); /*standard deviation*/
	arma::rowvec ll(n);

	NumericVector thetanew = Concvecelem(br, bh, hs, rs, gamma, phi);

	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);

		sdH[i] = -pow(gamma, 2) * (sum(hs % caphrCpp(y, y[i])) * exp(sum(bh % Z) + gamma * uh[i])) -
			(sum(rs % caphrCpp(y, y[i])) * exp(sum(br % Z) + uh[i])) - pow(phi, -2);

		arma::rowvec s(10);

		for (int j = 0; j < 10; j++) {
			s[j] = w[j] * exp(pow(x[j], 2)) * gCpp(uh[i] + pow(2, 0.5) * pow(-sdH[i], -0.5) * x[j], thetanew, y, r, dl, z, R, i, l1, l2, l3, l4);
		}

		ll[i] = log(pow(2, 0.5) * pow(-sdH[i], -0.5) * sum(s));
	}

	return (-sum(ll));
}

