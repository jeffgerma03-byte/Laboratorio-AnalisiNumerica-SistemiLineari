% autovalori_ex1_jordan.m
% Esercizio 1 - Blocco di Jordan perturbato
% Matricola: 5669424 -> d1=2, d0=4 -> n = 10(d1+1)+d0 = 34.

fprintf(">>> Esercizio 1: Blocco di Jordan perturbato\n");

d0 = 4;  % ultima cifra
d1 = 2;  % penultima cifra
n  = 10*(d1 + 1) + d0;

A = diag(ones(1,n-1), 1) + eye(n);

E = zeros(n);
E(n,1) = 2^(-n);
B = A + E;

VA = eig(A);
VB = eig(B);

rel_mat = norm(B - A) / norm(A);
rel_eig = norm(VB - VA) / norm(VA);

fprintf("n = %d\n", n);
fprintf("E(n,1) = 2^{-n} = %.4e\n", E(n,1));
fprintf("norm(B-A)/norm(A) = %.16e\n", rel_mat);
fprintf("norm(VB-VA)/norm(VA) = %.16e\n", rel_eig);

VA_sorted = sort(VA);
[~, idx] = sortrows([real(VB), imag(VB)]);
VB_sorted = VB(idx);

fprintf("\nAutovalori A (primi 10):\n");
disp(VA_sorted(1:min(10,end)));

fprintf("Autovalori B (primi 10, ordinati per Re/Im):\n");
disp(VB_sorted(1:min(10,end)));

AtA = A'*A;
BtB = B'*B;

VAtA = eig(AtA);
VBtB = eig(BtB);

rel_mat2 = norm(BtB - AtA) / norm(AtA);
rel_eig2 = norm(VBtB - VAtA) / norm(VAtA);

fprintf("\n--- Ripetizione su A^T A e B^T B ---\n");
fprintf("norm(BtB-AtA)/norm(AtA) = %.16e\n", rel_mat2);
fprintf("norm(eig(BtB)-eig(AtA))/norm(eig(AtA)) = %.16e\n", rel_eig2);
