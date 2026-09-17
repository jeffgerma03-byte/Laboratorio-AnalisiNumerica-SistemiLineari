% autovalori_ex2_grafo.m
% Esercizio 2 - Grafo ferrovie lombarde

fprintf(">>> Esercizio 2: Ferrovie lombarde (grafo)\n");

names = { ...
    "Milano", "Pavia", "Lodi", "Brescia", "Bergamo", ...
    "Como", "Varese", "Lecco", "Sondrio", "Cremona", "Mantova" };

n = 11;
A = zeros(n);

% Archi (assunti bidirezionali dalle connessioni in figura)
edges = [ ...
    1 2;  % Milano - Pavia
    1 3;  % Milano - Lodi
    1 4;  % Milano - Brescia
    1 5;  % Milano - Bergamo
    1 6;  % Milano - Como
    1 7;  % Milano - Varese
    6 7;  % Como - Varese
    6 8;  % Como - Lecco
    8 9;  % Lecco - Sondrio
    6 5;  % Como - Bergamo
    8 5;  % Lecco - Bergamo
    5 4;  % Bergamo - Brescia
    4 10; % Brescia - Cremona
    3 10; % Lodi - Cremona
    10 11 % Cremona - Mantova
];

for k = 1:size(edges,1)
    i = edges(k,1); j = edges(k,2);
    A(j,i) = 1; % i -> j
    A(i,j) = 1; % j -> i
end

g = sum(A,1);
D = diag(g);
G = A / D;

fprintf("n = %d nodi\n", n);
fprintf("Out-degree (g_j):\n");
disp(g);

[V, Lambda] = eig(G);
lam = diag(Lambda);

[~, idx1] = min(abs(lam - 1));
lambda1 = lam(idx1);
x = V(:, idx1);

x = real(x);
x = abs(x);
x = x / max(x);

fprintf("\nAutovalore piu' vicino a 1: %.16e\n", lambda1);
fprintf("Massimo modulo altri autovalori: %.16e\n", max(abs(lam([1:idx1-1 idx1+1:end]))));

fprintf("\nVettore x (normalizzato in [0,1]) associato a lambda ~ 1:\n");
for i = 1:n
    fprintf("(%2d) %-8s : %.6f\n", i, names{i}, x(i));
end

fprintf("\nVerifiche (punto c):\n");
fprintf("- |lambda1-1| = %.3e\n", abs(lambda1-1));
fprintf("- Tutti gli altri |lambda| < 1? (tolleranza 1e-12): %d\n", all(abs(lam([1:idx1-1 idx1+1:end])) < 1 - 1e-12));

mixed_sign = 0;
for k = 1:n
    if k == idx1, continue; end
    vk = real(V(:,k));
    if any(vk > 0) && any(vk < 0)
        mixed_sign = mixed_sign + 1;
    end
end
fprintf("- Autovettori degli altri autovalori con segni misti (conteggio): %d su %d\n", mixed_sign, n-1);
