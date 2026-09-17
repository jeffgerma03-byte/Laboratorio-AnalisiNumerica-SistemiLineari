% autovalori_ex3_potenze.m
% Esercizio 3 - Metodo delle potenze e potenze inverse (shift)

fprintf(">>> Esercizio 3: Metodo delle potenze e potenze inverse\n");

A = [ 1 -1  2;
     -2  0  5;
      6 -3  6];

y0a = [1; 1; 1];
y0b = [3; 10; 4];

maxit = 500;
tol   = 1e-12;

fprintf("Matrice A:\n");
disp(A);

fprintf("\n--- Metodo delle potenze (y0 = [1;1;1]) ---\n");
[lam_a, ~, it_a, res_a] = power_method_residual(A, y0a, maxit, tol);
fprintf("lambda (stima finale) = %.16e\n", lam_a);
fprintf("iterazioni = %d\n", it_a);
fprintf("residuo finale ||A y - lambda y||_2 = %.3e\n", res_a);

fprintf("\n--- Metodo delle potenze (y0 = [3;10;4]) ---\n");
[lam_b, ~, it_b, res_b] = power_method_residual(A, y0b, maxit, tol);
fprintf("lambda (stima finale) = %.16e\n", lam_b);
fprintf("iterazioni = %d\n", it_b);
fprintf("residuo finale ||A y - lambda y||_2 = %.3e\n", res_b);

p = lam_b;

fprintf("\n--- Metodo delle potenze inverse con shift p = %.16e ---\n", p);
[lam_inv, ~, it_inv, res_inv] = inverse_power_shift_residual(A, y0b, p, maxit, tol);
fprintf("lambda (stima finale, inverse shift) = %.16e\n", lam_inv);
fprintf("iterazioni = %d\n", it_inv);
fprintf("residuo finale ||A y - lambda y||_2 = %.3e\n", res_inv);

lam_exact = eig(A);
[~, idxmax] = max(abs(lam_exact));
fprintf("\n(Verifica) autovalore di massimo modulo da eig: %.16e\n", lam_exact(idxmax));

function [lambda, y, it, res] = power_method_residual(A, y0, maxit, tol)
    y = y0 / norm(y0,2);

    for k = 1:maxit
        w = A*y;
        y = w / norm(w,2);

        lambda = (y'*(A*y)) / (y'*y);

        r = A*y - lambda*y;
        res = norm(r,2);

        if res < tol
            it = k;
            return;
        end
    end

    it = maxit;
end

function [lambda, y, it, res] = inverse_power_shift_residual(A, y0, p, maxit, tol)
    n = size(A,1);
    M = A - p*eye(n);

    y = y0 / norm(y0,2);

    for k = 1:maxit
        w = M \ y;
        y = w / norm(w,2);

        lambda = (y'*(A*y)) / (y'*y);

        r = A*y - lambda*y;
        res = norm(r,2);

        if res < tol
            it = k;
            return;
        end
    end

    it = maxit;
end
