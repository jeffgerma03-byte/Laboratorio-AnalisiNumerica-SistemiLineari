% svd_ex2.m
% Esercizio 2 - Matrice triangolare superiore B e condizionamento
%
% b_{i,j} = 1  se i=j
%         = -1 se i<j
%         = 0  se i>j
%
% Si studiano sigma_max, sigma_min e cond_2(B) per n crescenti.
% Si perturba b_{n,1} della quantità -2^(2-n) e si calcolano gli autovalori.

fprintf(">>> Esercizio 2: Matrice triangolare superiore B\n");

n_list = [5, 10, 20, 40, 80];

sig_max = zeros(size(n_list));
sig_min = zeros(size(n_list));
cond2   = zeros(size(n_list));

for k = 1:length(n_list)
    n = n_list(k);

    B = triu(-ones(n), 1) + eye(n);

    s = svd(B);
    sig_max(k) = max(s);
    sig_min(k) = min(s);
    cond2(k)   = sig_max(k)/sig_min(k);

    fprintf("\nn = %d\n", n);
    fprintf("sigma_max = %.6e\n", sig_max(k));
    fprintf("sigma_min = %.6e\n", sig_min(k));
    fprintf("cond_2    = %.6e\n", cond2(k));

    % Perturbazione di b_{n,1} della quantità -2^(2-n)
    Bp = B;
    Bp(n,1) = Bp(n,1) - 2^(2-n);

    lam = eig(Bp);
    [min_abs, idx] = min(abs(lam));
    fprintf("min |eig(B_pert)| = %.6e (lambda = %.6e)\n", min_abs, lam(idx));

    % Controllo ulteriore: minimo valore singolare della perturbata
    sp = svd(Bp);
    fprintf("min sigma(B_pert) = %.6e\n", min(sp));
end

figure;
semilogy(n_list, sig_max, '-o', n_list, sig_min, '-o', n_list, cond2, '-o');
grid on;
xlabel('n');
ylabel('valore (scala log)');
legend('\\sigma_{max}', '\\sigma_{min}', 'cond_2(B)', 'Location', 'best');
title('Andamento di \\sigma_{max}, \\sigma_{min} e cond_2(B)');
