% svd_ex1.m
% Esercizio 1 - SVD di una matrice m x 3
%
% m fissato dalla matricola del primo componente (ordine alfabetico).
% d0 = ultima cifra, d1 = penultima cifra, m = 10(d0+1) + d1.
% Per 5669424: d0=4, d1=2 -> m = 52.

fprintf(">>> Esercizio 1: SVD di A (m x 3) e A^T\n");

d0 = 4;  % ultima cifra
d1 = 2;  % penultima cifra
m  = 10*(d0 + 1) + d1;

x = (1:m)'/m;                 % x_i = i/m
A = [ones(m,1), x, x.^2];     % A(i,:) = [1, x_i, x_i^2]

fprintf("m = %d\n", m);

% SVD di A e A^T
% uso SVD completa per poter descrivere Ker(A^T)
[U,S,V]   = svd(A);
[~,St,~]  = svd(A', 'econ');

sA  = diag(S);
sA  = sA(sA>0);               % valori singolari non nulli (A ha rango 3)
sAt = diag(St);

fprintf("\nValori singolari di A:\n");
disp(sA');

fprintf("Valori singolari di A^T:\n");
disp(sAt');

fprintf("Max |s(A) - s(A^T)| = %.3e\n", max(abs(sA - sAt)));

% Confronto con autovalori di AA^T e A^T A
AA = A*A';
AtA = A'*A;

eig_AA  = sort(eig(AA), 'descend');
eig_AtA = sort(eig(AtA), 'descend');

fprintf("\nAutovalori di A*A^T (primi 6):\n");
disp(eig_AA(1:min(6,end))');

fprintf("Autovalori di A^T*A:\n");
disp(eig_AtA');

fprintf("Quadrati dei valori singolari (s^2):\n");
disp((sA.^2)');

% Confronto immagine (orth) con vettori singolari sinistri/destri
QimA  = orth(A);
QimAt = orth(A');

fprintf("\nDimensione Im(A) (rank) da orth: %d\n", size(QimA,2));
fprintf("Dimensione Im(A^T) (rank) da orth: %d\n", size(QimAt,2));

r = rank(A);
U_r = U(:,1:r);
V_r = V(:,1:r);

P1 = QimA*QimA';
P2 = U_r*U_r';
fprintf("||P_orth(A) - P_U||_2 = %.3e\n", norm(P1-P2,2));

P1t = QimAt*QimAt';
P2t = V_r*V_r';
fprintf("||P_orth(A^T) - P_V||_2 = %.3e\n", norm(P1t-P2t,2));

% Nuclei
NA  = null(A);    % in R^3
NAt = null(A');   % in R^m

fprintf("\nDim Ker(A) = %d\n", size(NA,2));
fprintf("Dim Ker(A^T) = %d\n", size(NAt,2));

% Per SVD completa: Ker(A^T) e' generato dalle colonne di U corrispondenti a valori singolari nulli
U0 = U(:, r+1:end);   % base di Ker(A^T), dimensione m-r = 49

Pn1 = NAt*NAt';
Pn2 = U0*U0';
fprintf("||P_null(A^T) - P_U0||_2 = %.3e\n", norm(Pn1-Pn2,2));
