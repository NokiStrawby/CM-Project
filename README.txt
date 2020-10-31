README per il progetto di Computational Mathematics 
Giulia Fois, Nicolas Manini

Il codice scritto è diviso in due directory.

La directory Matlab contiene tutti i codici sorgenti scritti in MATLAB per l'implementazione dell'algoritmo e per i test svolti sui dataset.
L'implementazione dell'algoritmo principale è contenuta nel file LowRankAlgo.m 
L'implementazione della fattorizzazione QR con column pivoting è contenuta in MyQRP.m 
L'implementazione della variante con Random Restart è contenuta in LowRankAlgoRR.m 

Gli script GenerateRandomDataset.m e GenerateRandomDatasetSparse.m permettono di generare dei dataset di matrici casuali, eventualmente sparse.
Lo script ParseDataset.m permette di eseguire l'algoritmo su tutte le matrici salvate dagli script precedenti contenute in una cartella.
Similmente lo script ParseDataset_RR.m esegue la variante Random Restart dell'algoritmo sulle marici salvate dagli script precedenti.

Gli script Benchmark_Test.m e Benchmark_Test_RR.m sono stati utilizzati per eseguire i test sul benchmark dataset.
Lo script Img_Approximation.m è stato utilizzato per eseguire i test di comparazione tra truncated SVD e approssimazione data dal nostro algoritmo.
Lo script Img_Approximation_StopIt.m è stato utilizzato per verificare il numero di iterazioni dell'algoritmo sufficienti ad ottenere una buona approssimazione di immagine.
Gli script ParseDataset_RandomSubset.m e DoublingSize_GradNormAnalysis.m sono stati utilizzati per i test relativi alle stop condition.

Informazioni dettagliate sui parametri delle funzioni sono inserite come commenti nei relativi file.

La directory Python contiene i notebook relativi alla parte di data visualization. Abbiamo inserito tutti i notebook utilizzati per generare i plot contenuti nel report.
