#!/usr/bin/env Rscript

rm(list = ls(all = TRUE)) #Permet de remove l'environnement pour ne pas que les vecteurs gardent ceux du programme précédent.
##############FONCTIONS##############
#Fonction qui va écrire les 100 séquences dans un fichier texte au format fasta
ecriture<-function(taille,proba,filename) {
  for (i in 1:100){
    nom<-paste(i,"_seq_1",sep="")
    echantillon<-sample(c("A","C","G","T"),taille,replace<-T,prob<-proba)
    write.fasta(sequences<-echantillon, name<-nom,file.out<-filename,open<-'a')
  }
}
#Fonction qui calcule la distribution en nucléotide d'une séquence, et renvoi 4 valeurs :
prob<-function(masequence){
  compo<-table(masequence)
  probtot<-table(masequence)/length(masequence)
  return (probtot)
}
#####################################
###INTERACTION UTILISATEUR (input)###
fasta1<-readline(prompt<-"Indiquez le chemin du premier fichier : ")#vecteur qui prend la valeur de la destination du fichier 1 entré par l'utilisateur
fasta2<-readline(prompt<-"Indiquez le chemin du premier fichier : ")#vecteur qui prend la valeur de la destination du fichier 2 entré par l'utilisateur
#####################################

library("seqinr")#On déclare qu'on utilise seqinr pour lire les fichiers fasta

#Puis on lis les fichiers fasta, on les enregistre dans un vecteur, et on selectionne uniquement la séquence avec [[1]] :
#Séquence 1
tseq1<-read.fasta(file<-fasta1)
seq1<-tseq1[[1]]
#Séquence 2
tseq2<-read.fasta(file<-fasta2)
seq2<-tseq2[[1]]

#On détermine la taille des séquence pour créer des echantillons de la même taille :
taille1<-length(seq1)#Nombre de nucléotides dans la séquence 1
taille2<-length(seq2)#Nombre de nucléotides dans la séquence 2

#On utilise la fonction prob pour calculer la distribution en nucleotide de chaque séquence :
proba1<-prob(seq1)#proba1 prend les valeur de probabilitées de la séquence 1
proba2<-prob(seq2)#proba2 prend les valeur de probabilitées de la séquence 2
#On génère ensuite un nombre aléatoire entre 100 et 10000, afin d'obtenir un nom de fichier sortant unique :
nombre<-sample(100:10000,1)
#On additione ce nombre avec "_echantillon1.txt" et "_echantillon2.txt" pour obtenir le nom des deux fichiers :
filename1<-paste(nombre, "_echantillon1.txt", sep="")#Nom du fichier contenant les 100 seq issues de seq1
filename2<-paste(nombre, "_echantillon2.txt", sep="")#Nom du fichier contenant les 100 seq issues de seq2
#Enfin on appelle la fonction ecriture, qui va écrire 2*100 fois une séquence issue de seq1 et seq2, dans deux fichiers séparées
ecriture(taille1,proba1,filename1)#On écris dans un premier fichier, 100 seq de seq1
ecriture(taille2,proba2,filename2)#On écris dans un second fichier, 100 seq de seq2

#On print le nom des deux fichiers à l'utilisateur pour qu'il sache quel nom porte ses deux fichiers :
cat(filename1,' et ',filename2,'produits',sep<-"")
