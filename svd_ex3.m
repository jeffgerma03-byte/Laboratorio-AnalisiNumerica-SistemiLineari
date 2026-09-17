% svd_ex3.m
% Esercizio 3 - Minimi quadrati: confronto SVD, QR, equazioni normali, backslash
%
% Usa la stessa matrice A dell'esercizio 1 e y_i = sin(x_i).

fprintf(">>> Esercizio 3: Minimi quadrati Ac = y\n");

d0 = 4;  % ultima cifra
d1 = 2;  % penultima cifra
m  = 10*(d0 + 1) + d1;

x = (1:m)'/m;
A = [ones(m,1), x, x.^2];
y = sin(x);

% 1) SVD (pseudoinversa)
[U,S,V] = svd(A, 'econ');
s = diag(S);
tol = max(size(A))*eps(max(s));
s_inv = zeros(size(s));
s_inv(s > tol) = 1./s(s > tol);
A_pinv = V*diag(s_inv)*U';
c_svd = A_pinv * y;

% 2) QR
[Q,R] = qr(A, 0);
c_qr = R \ (Q' * y);

% 3) Equazioni normali
c_ne = (A'*A) \ (A'*y);

% 4) Backslash
c_bs = A \ y;

sol = [c_svd, c_qr, c_ne, c_bs];
names = ["SVD", "QR", "NormalEq", "Backslash"];

fprintf("\nSoluzioni c (colonne: SVD, QR, NormalEq, Backslash):\n");
disp(sol);

fprintf("Confronto ||c_i - c_j||_2:\n");
for i = 1:4
    for j = i+1:4
        fprintf("%s vs %s: %.3e\n", names(i), names(j), norm(sol(:,i)-sol(:,j),2));
    end
end

fprintf("\nResidui ||A c - y||_2:\n");
for i = 1:4
    r = norm(A*sol(:,i) - y, 2);
    fprintf("%s: %.6e\n", names(i), r);
end
