#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "updateRP_genomic_cpp.h"

// [[Rcpp::export]]
arma::mat matProdVec(const arma::mat x, const arma::vec y)
{
    // multiply (element-wise) a matrix to a expanded vector
    
    arma::mat mat_y = arma::zeros<arma::mat>(y.n_elem, x.n_cols);
    mat_y.each_col() = y;
    arma::mat spanMat = x % mat_y; // elementwise product
    return spanMat;
}

// [[Rcpp::export]]
arma::vec sumMatProdVec(const arma::mat x, const arma::vec y)
{
    // compute "arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );"
    
    arma::vec spanVec = arma::zeros(x.n_cols);
    for (unsigned int i = 0; i < x.n_cols; ++i)
        spanVec(i) = arma::dot(x.col(i), y);
    return spanVec;
}

// [[Rcpp::export]]
Rcpp::List updateRP_genomic_cpp(const unsigned int p,
                                const arma::mat x_,
                                const unsigned int J_,
                                arma::mat ind_r_,
                                arma::mat ind_d_,
                                arma::mat ind_r_d_,
                                arma::vec be_,
                                arma::vec ga_,
                                arma::vec h_,
                                const double tau,
                                const double cb)
{
    // update coefficients of genomic variables via a MH sampler
    
    arma::uvec updatej = arma::randperm(p);
    
    arma::vec xbeta_ = x_ * be_;
    arma::vec sd_be_ = arma::ones<arma::vec>(ga_.n_elem);
    sd_be_.elem(arma::find(ga_ == 1.)).fill(cb);
    sd_be_ = sd_be_ * tau;
    arma::uvec sampleRPg_accept_ = arma::zeros<arma::uvec>(p);
    
    unsigned int j = 0;
    for (unsigned int j_id = 0; j_id < p; ++j_id)
    {
        j = updatej(j_id);
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        exp_xbeta = arma::exp(xbeta_);
        x_exp_xbeta = x_.col(j) % exp_xbeta;
        // D1_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta ).t(), 1 );
        D1_1st = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta);

        h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
        h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
        exp_h_exp_xbeta_mat = arma::exp(h_exp_xbeta_mat);
        D1_2nd_den = 1. - exp_h_exp_xbeta_mat;
        D1_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_exp_xbeta);
        D1_2nd = h_ % arma::sum((D1_2nd_num / D1_2nd_den % ind_d_).t(), 1);
        D1 = arma::sum(D1_1st + D1_2nd) - 1. / sd_be_(j) / sd_be_(j) * be_(j);

        x_sq_exp_xbeta = x_.col(j) % x_.col(j) % exp_xbeta;
        // D2_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta ).t(), 1 );
        D2_1st = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta);
        D2_2nd_den = D1_2nd_den % D1_2nd_den;
        D2_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_sq_exp_xbeta) % (1. - exp_h_exp_xbeta_mat + h_exp_xbeta_mat);
        D2_2nd = h_ % arma::sum((D2_2nd_num / D2_2nd_den % ind_d_).t(), 1);
        D2 = arma::accu(D2_1st + D2_2nd) - 1. / sd_be_(j) / sd_be_(j);

        be_prop_me = be_(j) - D1 / D2;
        be_prop_sd = 2.4 / sqrt(-D2);
        be_prop = be_;

        // genomic version:
        // be_prop(j) = R::rnorm( be_prop_me, be_prop_sd );
        be_prop(j) = arma::randn(arma::distr_param(be_prop_me, be_prop_sd));
        xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop(j);
        xbeta_prop.elem(arma::find(xbeta_prop > 700)).fill(700.);
        exp_xbeta_prop = arma::exp(xbeta_prop);
        x_exp_xbeta_prop = x_.col(j) % exp_xbeta_prop;
        // D1_1st_prop = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta_prop ).t(), 1 );
        D1_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta_prop);

        h_exp_xbeta_prop_mat = -arma::kron(exp_xbeta_prop, h_.t());
        h_exp_xbeta_prop_mat.elem(arma::find(h_exp_xbeta_prop_mat > -1.0e-7)).fill(-1.0e-7);
        exp_h_exp_xbeta_prop_mat = arma::exp(h_exp_xbeta_prop_mat);
        D1_2nd_den_prop = 1. - exp_h_exp_xbeta_prop_mat;
        D1_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_exp_xbeta_prop);
        D1_2nd_prop = h_ % arma::sum((D1_2nd_num_prop / D1_2nd_den_prop % ind_d_).t(), 1);
        D1_prop = arma::accu(D1_1st_prop + D1_2nd_prop) - 1. / sd_be_(j) / sd_be_(j) * be_prop(j);

        x_sq_exp_xbeta_prop = x_.col(j) % x_.col(j) % exp_xbeta_prop;
        // D2_1st_prop = -h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta_prop ).t(), 1);
        D2_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta_prop);
        D2_2nd_den_prop = D1_2nd_den_prop % D1_2nd_den_prop;
        D2_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_sq_exp_xbeta_prop) % (1. - exp_h_exp_xbeta_prop_mat + h_exp_xbeta_prop_mat);
        D2_2nd_prop = h_ % arma::sum((D2_2nd_num_prop / D2_2nd_den_prop, ind_d_).t(), 1);
        D2_prop = arma::accu(D2_1st_prop + D2_2nd_prop) - 1. / sd_be_(j) / sd_be_(j);
        be_prop_me_ini = be_prop(j) - D1_prop / D2_prop;
        be_prop_sd_ini = 2.4 / sqrt(-D2_prop);

        // first_sum = arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );
        first_sum = sumMatProdVec(ind_r_d_, exp_xbeta);
        second_sum = arma::sum((arma::log(D1_2nd_den) % ind_d_).t(), 1);

        loglh_ini = arma::accu(-h_ % first_sum + second_sum);
        // first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1) ;
        first_sum_prop = sumMatProdVec(ind_r_d_, exp_xbeta_prop);
        second_sum_prop = arma::sum((arma::log(D1_2nd_den_prop) % ind_d_).t(), 1);
        loglh_prop = arma::accu(-h_ % first_sum_prop + second_sum_prop);

        /*logprior_prop = R::dnorm( be_prop(j), 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop(j), be_prop_me_ini, be_prop_sd_ini, true);
        logprop_ini = R::dnorm( be_(j), be_prop_me, be_prop_sd, true);*/
        logprior_prop = arma::log_normpdf(be_prop(j), 0.0, sd_be_(j));
        logprior_ini = arma::log_normpdf(be_(j), 0.0, sd_be_(j));
        logprop_prop = arma::log_normpdf(be_prop(j), be_prop_me_ini, be_prop_sd_ini);
        logprop_ini = arma::log_normpdf(be_(j), be_prop_me, be_prop_sd);
        logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;

        // if( log( R::runif(0., 1.) ) < logR )
        if (log(arma::randu()) < logR)
        {
            be_(j) = be_prop(j);
            xbeta_ = xbeta_prop;
            sampleRPg_accept_(j) = sampleRPg_accept_(j) + 1;
        }
    }
    
    return Rcpp::List::create(
        Rcpp::Named("be.ini") = be_,
        Rcpp::Named("acceptl") = sampleRPg_accept_
      );
}
