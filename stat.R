#!/usr/bin/env Rscript
###INTERACTION UTILISATEUR (input)###
fichierscore<-readline(prompt<-"Indiquez le chemin du fichier 'scores.txt' : ")#vecteur qui prend la valeur de la destination du fichier scores.txt
mainscore<-readline(prompt<-"Indiquez le score obtenu lors de l'alignement des deux séquences principales : ")#Enregistre le score des deux séquences d'origine
alpha<-readline(prompt<-"Indiquez une valeur pour alpha (ex : 0.1) : ")
#####################################
#####################################
#########TEST D'INDEPENDANCE#########
#H0 : « les deux sequences ont des distributions de lettres independantes » contre l'alternative
#H1 : « les distributions des deux sequences sont liees ».
###CALCUL DE LA P-VALUE DE L'ALIGNEMENT###
scores<-scan(fichierscore)#Scan le fichier et enregistre les scores dans un vecteur
scoresup<-sum(scores>as.numeric(mainscore))#Calcule le nombre de scores parmis les scores des 100 séquences simulées lesquels sont supérieur au score des séquences d'origine
nombrescore<-length(scores)#Calcule le nombre de scores dans le fichier (normalement 100)

pvalue<-scoresup/nombrescore#Va calculer la p-value entre le score principal et ceux générés par les séquences aléatoires
#Donne le résultat de la p-value du score et le résultat
cat('La p-value observée du score est de',pvalue,'\n')
plot(density(scores), xlab="Valeur des scores", ylab="densité")
if (pvalue>alpha){
  cat('On ne rejette pas H0')
}else{
  cat('On rejette H0')
}