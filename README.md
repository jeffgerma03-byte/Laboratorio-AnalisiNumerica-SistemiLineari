# Analisi Numerica e Sistemi Lineari in C++ e MATLAB

Una suite di esercitazioni sviluppate in C++ e MATLAB per il corso di Algebra Lineare e Analisi Numerica. Il progetto esplora in modo empirico e analitico i limiti dell'aritmetica floating-point (IEEE 754), la stabilità degli algoritmi risolutivi per i sistemi lineari, la sensibilità del problema spettrale e le proprietà della decomposizione ai valori singolari (SVD).

Il repository include le implementazioni sorgente degli algoritmi, le esecuzioni dei test parametrizzati e le relazioni tecniche che documentano e commentano i risultati ottenuti.

## Moduli del Progetto

Il progetto è diviso in quattro macro-aree:

### 1. Calcolo Numerico ed Errori di Arrotondamento (`Esercitazione.cpp`)

Un'analisi pratica della propagazione degli errori macchina:

* **Catastrofi di Cancellazione e Associatività**: Dimostrazione di come l'ordine delle operazioni `(a+b)+c` contro `a+(b+c)` generi risultati completamente diversi a causa dell'arrotondamento in doppia precisione su ordini di grandezza estremi.
* **Sviluppo di Taylor e Metodo di Horner**: Approssimazione polinomiale di `exp(x)`. Valutazione dell'errore relativo al variare del numero di termini (N) e gestione ottimizzata del reciproco per esponenti negativi per prevenire la perdita di cifre significative.
* **Epsilon di Macchina**: Calcolo empirico a runtime della precisione di macchina per i tipi `float` e `double`.

### 2. Risoluzione di Sistemi Lineari (`EsercitazioneSistemi.cpp`)

Implementazione e stress-test di algoritmi per la risoluzione di sistemi `Ax = b`:

* **Eliminazione di Gauss con Pivoting Parziale**: Algoritmo sviluppato da zero per minimizzare l'amplificazione degli errori durante la risoluzione.
* **Analisi del Condizionamento**: Test condotti su matrici strutturate per misurare l'affidabilità delle soluzioni:
   * Matrice di Pascal (10x10): Utilizzata per dimostrare il crollo numerico in scenari di forte mal condizionamento.
   * Matrici Tridiagonali (Diagonale Dominante): Sfruttate per testare scenari ad alta stabilità strutturale.
* **Analisi delle Perturbazioni**: Iniezione di un rumore controllato (perturbazione dell'1%) nel termine noto per misurare empiricamente l'esplosione dell'errore (calcolato tramite norma infinito) in relazione al condizionamento della matrice.

### 3. Autovalori, Autovettori e Metodi Iterativi (MATLAB)

Analisi della sensibilità del problema spettrale e dei metodi iterativi per il calcolo degli autovalori, tramite gli script `autovalori_ex1_jordan.m`, `autovalori_ex2_grafo.m`, `autovalori_ex3_potenze.m` (eseguibili in sequenza tramite `run_all_autovalori.m`):

* **Sensibilità degli Autovalori (Blocco di Jordan)**: Costruzione di una matrice `A` e di una sua perturbazione `B = A + E` per dimostrare come, in presenza di matrici difettive/non diagonalizzabili, perturbazioni minime (norma relativa dell'ordine di 10⁻¹¹) possano generare variazioni relative enormi sugli autovalori (~0.5). Ripetizione dell'analisi su `AᵀA` e `BᵀB` per mostrare la maggiore regolarità numerica delle matrici simmetriche positive.
* **Autovalori di una Matrice Associata a un Grafo**: Costruzione di una matrice di transizione `G = AD⁻¹` a partire dalla matrice di adiacenza di un grafo ferroviario (11 nodi) e dai gradi uscenti. Individuazione dell'autovalore `λ ≈ 1` e interpretazione del relativo autovettore non negativo come indice di importanza relativa dei nodi (logica alla base dei risultati di tipo Perron-Frobenius).
* **Metodi Iterativi**: Implementazione del metodo delle potenze (con confronto tra vettori iniziali diversi) e del metodo delle potenze inverse con shift, valutando velocità di convergenza e residuo `‖Ay − λy‖₂` per verificare l'accuratezza numerica delle stime.

### 4. SVD e Minimi Quadrati (MATLAB)

Studio della decomposizione ai valori singolari (SVD) e delle sue applicazioni numeriche, tramite gli script `svd_ex1.m`, `svd_ex2.m`, `svd_ex3.m` (eseguibili in sequenza tramite `run_all_svd.m`):

* **SVD e Sottospazi Fondamentali**: Confronto tra i valori singolari di una matrice `A` e della sua trasposta `Aᵀ`, verifica del legame tra i quadrati dei valori singolari e gli autovalori non nulli di `AᵀA` e `AAᵀ`, e verifica dei sottospazi fondamentali `Im(A)`, `Ker(A)`, `Im(Aᵀ)`, `Ker(Aᵀ)` tramite confronto con `orth` e `null`.
* **Condizionamento e Quasi-Singolarità**: Analisi di una matrice triangolare superiore `B` al crescere della dimensione `n`, con calcolo di `σ_max`, `σ_min` e `cond₂(B)` (fino a un condizionamento dell'ordine di 10¹⁹ per n = 80). Studio dell'effetto di una perturbazione mirata sulla quasi-singolarità e sul rango numerico, mostrando come la SVD sia più robusta di `eig` nell'individuare la perdita di rango.
* **Minimi Quadrati**: Risoluzione del problema `Ac ≈ y` con quattro metodi differenti (SVD, QR, equazioni normali, operatore `A\y`), confrontando le soluzioni ottenute e il residuo `‖Ac − y‖₂` per valutarne la robustezza numerica.

## Tecnologie Utilizzate

* **Linguaggi**: C++, MATLAB
* **Ambienti**: Compilatore C++ standard, MATLAB (con visualizzazione grafica per l'esercizio sul condizionamento)

## Autore

* Jeffrey Germano (Matricola: 5669424)
