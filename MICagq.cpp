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

	return wrap(delta_pos); /*returns the terminal event indicator*/
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

}

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
arma::rowvec termlikCpp(NumericVector u, NumericMatrix z, NumericVector y, NumericVector dl, 
	              NumericVector the, int i, int l1, int l2, int l3, int l4) {
	
	NumericVector delta(dl);
	NumericVector theta(the);
	NumericMatrix Zr(z);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double gamma = theta[l1 + l2 + l3 + l4];

	int n = u.size();
	arma::rowvec ans(n);
	arma::rowvec Z = Zr(i, _);

	for (int j = 0; j < u.size(); j++) {
        ans[j] = pow(exp(sum(beta2 % Z) + gamma * u[j]) * sum(hs % smhrCpp(y, y[i])), delta[i])*
			exp(-sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * u[j]));
	}

	return ans;
}

//[[Rcpp::export]]
arma::rowvec reclikeCpp(NumericVector u, NumericVector the, NumericVector y, NumericVector r,
	NumericMatrix z, NumericMatrix R, int a, int l1, int l2, int l3, int l4) {
  /* returns a row vector*/
	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericMatrix Zr(z);
	NumericVector theta(the);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	
	arma::rowvec Z = Zr(a, _);
	int m = u.size();
	arma::rowvec ans(m);

	for (int j = 0; j < m; j++) {
		if (ni[a] > 0) {
			NumericVector recur_times = findCpp2(R, a + 1);
			double prod = 1;
			for (int k = 0; k < ni[a]; k++) {
				prod *= exp(sum(beta1 % Z) + u[j]) * sum(rs % smhrCpp(Tij, recur_times[k]));
			}
			ans[j] = prod * exp(-sum(rs % caphrCpp(y, y[a])) * exp(sum(beta1 % Z) + u[j]));
		}
		else {
		    ans[j] = exp(-sum(rs % caphrCpp(y, y[a])) * exp(sum(beta1 % Z) + u[j]));
		}
	}
	return ans;
}

//[[Rcpp::export]]
arma::rowvec dnormCpp(NumericVector x, double m, double v) {

	int n = x.size();
	arma::rowvec ans(n);

	for (int i = 0; i < n; i++) {
		ans[i] = R::dnorm(x[i], m, v, FALSE);
	}
   return ans;
}

//[[Rcpp::export]]
arma::rowvec MllCpp(NumericVector u, NumericVector the, NumericVector y, NumericVector dl, NumericVector r,
	NumericMatrix z, NumericMatrix R, int i, int l1, int l2, int l3, int l4) {

	NumericVector useq(u);
	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericVector Tij = R(_, 0);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	//int m = useq.size();
	double mn = 0;
	//arma::rowvec mll(m);
	//arma::rowvec Z = Zr(i, _);

	arma::rowvec mll = termlikCpp(useq, z, y, dl, the, i, l1, l2, l3, l4) % reclikeCpp(useq, the, y, r, z, R, i, l1, l2, l3, l4) % dnormCpp(useq, mn, phi);

	
	//for (int j = 0; j < m; j++) {
		//	mll[j] = termlikCpp(useq, z, y, dl, the, i, l1, l2, l3, l4) * reclikeCpp(useq, the, y, r, z, R, i, l1, l2, l3, l4) *
			//	dnormCpp(useq, mn, phi);
	//}


	return mll;
}

//[[Rcpp::export]]
double UhatCpp(NumericVector u, NumericVector the, NumericVector y, NumericVector dl, NumericVector r,
	NumericMatrix z, NumericMatrix R, int i, int l1, int l2, int l3, int l4) {

	NumericVector useq(u);
	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericVector Tij = R(_, 0);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	double mn = 0;

	arma::rowvec mll = termlikCpp(useq, z, y, dl, the, i, l1, l2, l3, l4) % reclikeCpp(useq, the, y, r, z, R, i, l1, l2, l3, l4) % dnormCpp(useq, mn, phi);

	arma::uword mx = arma::index_max(mll);

	//double ans = useq[mx];
	return useq[mx];
}

// [[Rcpp::export]]
arma::rowvec UhatVecCpp(NumericVector u, NumericVector the, NumericVector y, NumericVector dl, NumericVector r,
	NumericMatrix z, NumericMatrix R, int l1, int l2, int l3, int l4) {

	NumericVector useq(u);
	NumericVector theta(the);
	NumericMatrix Zr(z);
	NumericVector ni(r);
	NumericVector delta(dl);
	NumericVector Tij = R(_, 0);

	arma::rowvec beta1(theta.begin(), l1);
	arma::rowvec beta2(theta.begin() + l1, l2);
	arma::rowvec hs(theta.begin() + l1 + l2, l3);
	arma::rowvec rs(theta.begin() + l1 + l2 + l3, l4);
	double phi = theta[l1 + l2 + l3 + l4 + 1];

	double mn = 0;
	int n = ni.size();
	arma::rowvec uhat(n);
	
	for (int i = 0; i < n; i++) {

		arma::rowvec Z = Zr(i, _);
		arma::rowvec mll = termlikCpp(useq, z, y, dl, the, i, l1, l2, l3, l4) % reclikeCpp(useq, the, y, r, z, R, i, l1, l2, l3, l4) % dnormCpp(useq, mn, phi);
		arma::uword mx = arma::index_max(mll);
		uhat[i] = useq[mx];
	}

	return uhat;
}

//[[Rcpp::export]]
double termlikscCpp(double u, NumericMatrix z, NumericVector y, NumericVector dl, arma::rowvec beta2, arma::rowvec hs, double gamma, int i) {
   
	/*returns just a number*/
	NumericVector delta(dl);
	//NumericVector theta(the);
	NumericMatrix Zr(z);

	arma::rowvec Z = Zr(i, _);

	double ans = pow(exp(sum(beta2 % Z) + gamma * u) * sum(hs % smhrCpp(y, y[i])), delta[i]) *
		exp(-sum(hs % caphrCpp(y, y[i])) * exp(sum(beta2 % Z) + gamma * u));

	return ans;
}

//[[Rcpp::export]]
double reclikescCpp(double u, NumericVector y, NumericVector r, NumericMatrix z, NumericMatrix R, arma::rowvec beta1, arma::rowvec rs, int j) {

	/*returns just a number */
	NumericVector Tij = R(_, 0);
	NumericVector ni(r);
	NumericMatrix Zr(z);
	//NumericVector theta(the);

	arma::rowvec Z = Zr(j, _);
	double ans = 0;

	if (ni[j] > 0) {
		NumericVector recur_times = findCpp2(R, j + 1);
		double prod = 1;
		for (int k = 0; k < ni[j]; k++) {
			prod *= exp(sum(beta1 % Z) + u) * sum(rs % smhrCpp(Tij, recur_times[k]));
		}
		 ans = prod * exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u));
	}
	else {
		 ans = exp(-sum(rs % caphrCpp(y, y[j])) * exp(sum(beta1 % Z) + u));
	}

return ans;
}

//[[Rcpp::export]]
double dnormscCpp(double x, double m, double v) {

	return (R::dnorm(x, m, v, FALSE));

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
double elemmultsum(arma::rowvec a, arma::rowvec b) {

	return (sum(a % b));

}