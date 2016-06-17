function [prox] = prox_sum_ind_epi_conjquad_1d(idx, count, interleaved, a, b, c, alpha, beta)

    dim = 2;
    coeffs = { a, b, c, alpha, beta };
    prox = { 'ind_epi_conjquad_1d', idx, count*dim, false, ...
             { count, dim, interleaved, coeffs } };

end
