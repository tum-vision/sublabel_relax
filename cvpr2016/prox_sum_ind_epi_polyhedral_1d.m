function [prox] = prox_sum_ind_epi_polyhedral_1d(idx, count, interleaved, x, y, alpha, beta, ind_vec, cnt_vec)

    dim = 2;
    coeffs = { x, y, alpha, beta, cnt_vec, ind_vec };
    prox = { 'ind_epi_polyhedral_1d', idx, count*dim, false, ...
             { count, dim, interleaved, coeffs } };

end
