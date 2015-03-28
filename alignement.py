#!/usr/bin/python3
# -*- coding: utf8 -*-
#Programme : alignement.py
#Auteur : Vincent ROCHER
#But du programme : Alignement local(Smith et Waterman), avec ou sans tracé arrière.
#
##############FONCTIONS################
import sys #Module qui permet la vérification de la version de Python de l'utilisateur
if sys.version_info < (3, 0): #Si Version de python inférieur à 3 :
    exit() #Quitter python

else:
    ######MULTIFASTA###################
#Lecture du fichier fasta, récupère le nom et le contenu de la séquence
#de chaque séquence dans deux listes, puis renvoi les deux listes :
    def lireFastaMul(nomFichier):
        fichier=open(nomFichier,"r")#Ouverture en mode lecture
        lignes=fichier.readlines()#Récupère toutes les lignes du fichier
        fichier.close()#Fermeture du fichier
    #Declaration des variables et des listes :
        sequence=""
        nom=""
        sequences=[]
        noms=[]
    #Boucle de lecture
        for ligne in lignes:#Parcours chaque ligne du fichier
            if ligne[0]==">":#Si le premier élément de la ligne contient un '>'
                if sequence != "":#Enregistre seulement si sequence contient quelque chose, évite d'enregistrer le sequence="" du début
                    noms.append(nom)#Enregistre le nom dans la liste des noms
                    sequences.append(sequence)#Enregistre la séquence dans la liste des séquences
                nom=ligne[1:-1]#Enregistre le nom de la séquence sans le '>' et sans le saut de ligne
                sequence=""#Remet la valeur de la variable à 0 pour ne pas enregistrer deux séquences en même temps
            else:#Enregistre ligne dans la variable séquence et sans saut de ligne
                sequence=sequence+ligne[:-1]
        if sequence !="":#Permet d'enregistrer la dernière séquence car n'enregistre qu'après avoir vu un '>'
            noms.append(nom)
            sequences.append(sequence)
        return (noms,sequences)#Renvoi des deux listes contenant toutes les séquences

    ######TABLEAU DE SUBSTITUTION######
    def msubst(L,C):
#L=Nucléotide qui correspondra à la ligne
#C=Nucléotide qui correspondra à la colonne
#Y=score match déterminé par l'utilisateur
#N=score mismatch déterminé par l'utilisateur
#M=score mismatch mais même type de base
        global Y,N,M
        sub=[
            #A C G T W Y R N
            #0 1 2 3 4 5 6 7
            [Y,N,M,N,Y,N,Y,Y],#A/0
            [N,Y,N,M,N,Y,N,Y],#C/1
            [M,N,Y,N,N,N,Y,Y],#G/2
            [N,M,N,Y,Y,Y,N,Y],#T/3
            [Y,N,N,Y,Y,N,N,Y],#W/4
            [N,Y,N,Y,N,Y,N,Y],#Y/5
            [Y,N,Y,N,N,N,Y,Y],#R/6
            [Y,Y,Y,Y,Y,Y,Y,Y]]#N/7
        eq = {'A':0,'C':1,'G':2,'T':3,'W':4,'Y':5,'R':6,'N':7}#Dictionnaire qui remplace les lettres par leur numéro correspondant dans la matrice de substitution
        score=sub[eq[L]][eq[C]]
        return(float(score))

    ######ALIGNEMENT LOCAL#############    
    def Align(seq1,seq2):
    ###Construction Matrice#####
        m=len(seq1)+1 #m = nombre de lignes
        n=len(seq2)+1 #n = nombre de colonnes
        T=[0]*m #Construction de m lignes
        for i in range(m):
            T[i]=[0]*n #Construction de n colonnes
    #########Alignement#########    
        i=j=1
        maxi=0
        while i<m:
            while j<n:
            #Nécessité d'enlever -1 pour obtenir la bonne position de chaque nucleotides
                p1=i-1#position du nucleotide de la séquence 1
                p2=j-1#position du nucleotide de la séquence 2
                score1=T[i-1][j-1]+msubst(seq1[p1],seq2[p2])#Score d'un match/mismatch
                score2=T[i][j-1]+d#Insertion dans la séquence 1, délétion dans séquence 2
                score3=T[i-1][j]+d#Délétion dans la séquence 1, insertion dans séquence 1
                T[i][j]=max(score1,score2,score3,0)#Calcul du score maximum et insertion dans Table
                if T[i][j]>=maxi:
                    maxi=T[i][j]
                    imax=i#Coordonnées du score max au niveau des lignes (i)
                    jmax=j#Coordonnées du score max au niveau des colonnes (j)
                j=j+1
            i=i+1
            j=1#Remet j à 1 pour faire recommencer la boucle a la première colonne
        return (T,maxi,imax,jmax)
    ######TRACE ARRIERE###############
    def backtrac(seq1,seq2,T,imax,jmax):
        Align1="" #Contiendra le premier alignement
        Align2="" #Contiendra le second alignement
        Align3="" #Va servir a illustrer l'alignement dans le html par des . ou des |
        i=len(seq1[0:imax])#i va prendre comme valeur imax
        j=len(seq2[0:jmax])#j va prendre comme valeur jmax
        while (i>0) or (j>0): #Continu tant que i et j ne valent pas 0
            if (T[i][j]==0): #Sert a arrêter la boucle si on croise un score nul (0)
                i=0
            elif (i>0) and (j>0) and (T[i][j]==T[i-1][j-1]+msubst(seq1[i-1],seq2[j-1])): #Vérifie si le score est égal au score de la case diagonale
                Align1=seq1[i-1]+Align1
                Align2=seq2[j-1]+Align2
                if (seq1[i-1]==seq2[j-1]):
                    Align3="|"+Align3
                else:
                    Align3="."+Align3
                i=i-1
                j=j-1
            elif (i>0) and (T[i][j]==T[i-1][j]+d): #Vérifie si le score est égal au score de la case de gauche 
                Align1=seq1[i-1]+Align1
                Align2="-"+Align2
                Align3=" "+Align3
                i=i-1
            elif (j>0) and (T[i][j]==T[i][j-1]+d): #Vérifie si le score est égal au score de la case du haut
                Align1="-"+Align1
                Align2=seq2[j-1]+Align2
                Align3=" "+Align3
                j=j-1
        return(Align1,Align2,Align3)
    #######Print le résultat du tracé arrière dans une page result.html#### 
    #Pas eu le temps de détailler le code, désolé (résultat visible dans result.html)#   
    def printhtml(align1,align2,align3,nom1,nom2,score):
        global Y,N,d,nombreseq
        taille=len(align1)
        html_strdeb='''
        <!doctype html>
        <html>
        <head>
        <meta charset="UTF-8">
        <title>Local Alignment</title>
        <body>'''
        html_strdeb2='''
# Program: Smith & Waterman
# Autor : Rocher Vincent
# Align_format: pair
# Report_file: stdout
########################################
        
#=======================================
#
# Aligned_sequences: '''
        html_strdeb3='''
# 1: '''
        html_strdeb4='''
# 2: '''
        html_strdeb5='''
# Gap_penalty: '''
        html_strdeb6='''
# Match: '''
        html_strdeb7='''
# Mismatch:'''
        html_strdeb8=''' 
# Length: '''
        html_strdeb9='''
# Score: '''
        html_strdeb10='''
# 
#
#======================================='''
        html_strfin='''
        </body>
        </html>'''
        with open("result.html","w") as html: #Ecrit dans le fichier result.html
            n=1
            html.write(html_strdeb)
            html.write('''<pre id="alignmentContent" xml:space="preserve">''')
            html.write(html_strdeb2)
            html.write(str(nombreseq))
            html.write(html_strdeb3)
            html.write(nom1)
            html.write(html_strdeb4)
            html.write(nom2)
            html.write(html_strdeb5)
            html.write(str(d))
            html.write(html_strdeb6)
            html.write(str(Y))
            html.write(html_strdeb7)
            html.write(str(N))
            html.write(html_strdeb8)
            html.write(str(taille))
            html.write(html_strdeb9)
            html.write(str(score))
            html.write(html_strdeb10)
            html.write('<br />')
            while (n<=len(align1)):
                html.write('1: ')
                html.write(align1[n:n+30].replace('',' '))
                html.write('<br />')
                html.write('   ')
                html.write(align3[n:n+30].replace('',' '))
                html.write('<br />')
                html.write('2: ')
                html.write(align2[n:n+30].replace('',' '))
                html.write('<br />')
                html.write('<br />')
                n=n+31
            html.write('''</pre>''')
            html.write(html_strfin)
            print("le tracé arrière est visible dans le fichier result.html")
#####Declaration variables et interaction utilisateur######
    noms1=[]#Liste qui récupère les noms des séquences issue du fichier 1
    noms2=[]#Liste qui récupère les noms des séquences issue du fichier 2
    sequences1=[]#Liste qui va contenir les séquences du premier fichier
    sequences2=[]#Liste qui va contenir les séquences du second fichier

####################INTERACTION UTILISATEUR################
    print('''####################################################''')
    print('''#                                                  #''')
    print('''#  Programme d'alignement local Smith & Waterman   #''')
    print('''#  Réalisé dans le cadre du projet Programmation   #''')
    print('''#                En Bio-Informatique               #''')
    print('''#                Par ROCHER Vincent                #''')
    print('#           Version de python : ',sys.version[0:5],'            #''')
    print('''#                                                  #''')
    print('''####################################################''')
#Récupère la destination des deux fichiers et renvoi les séquences dans deux variables différentes :
    noms1,sequences1=lireFastaMul(input('Destination du fichier de la ou des séquence(s) 1: \n'))
    noms2,sequences2=lireFastaMul(input('Destination du fichier de la ou des séquence(s) 2: \n'))
    d=float(input('Score d\'indel : ')) #valeur du score en cas d'indel
    if d>0:
        d=-d
    Y=float(input('Score match : ')) #valeur du score en cas d'identité
    N=float(input('Score mismatch : ')) #valeur du score en cas de substitution
    M=float(input('Score mismatch même type de base purique/pyrimidique (possibilité de mettre même score que mismatch) : '))
    back=input('Voulez vous effectuer le tracé arrière ?(y/n) :')
#########Alignement de chaque séquence du fichier 1 avec chaque séquence du fichier 2#########
#Vérifie le nombre de séquences dans chaque fichier fasta, prend la plus petite série pour comparer
    if len(sequences1)<len(sequences2):
        nombreseq=len(sequences1)
    else:
        nombreseq=len(sequences2)


########LANCEMENT DES FONCTIONS############


    print('''Début de l'alignement ... Veuillez patienter.''')
    i=0
    scores=[]
    T=[]#Tableau qui va contenir le tableau de score
    while i<nombreseq: #Boucle qui va continuer tant qu'il y a des sequences a aligner
        sequences1[i]=sequences1[i].upper()
        sequences2[i]=sequences2[i].upper()
        print('Comparaison de :\n',noms1[i],'\nAvec :\n',noms2[i],'\nEn cours...',sep='')
        T,score,imax,jmax=Align(sequences1[i],sequences2[i]) #Renvoi le score d'alignement
        print('Sequences (1_',i+1,':2_',i+1,') Aligned. Score : ',score,sep='')
        scores.append(score) #Ajoute le score dans une liste
        if back=='y':
            print('Alignement de :\n',noms1[i],'\nAvec :\n',noms2[i],'\nEn cours...',sep='')
            align1,align2,align3=backtrac(sequences1[i],sequences2[i],T,imax,jmax) #Effectue le tracé arrière si 'y'
            printhtml(align1,align2,align3,noms1[i],noms2[i],scores[i]) #Affiche le résultat dans une page web
        T=[]
        i=i+1

    if nombreseq>1: #Ecrit dans un fichier texte les scores si il y en a plus d'un
        with open("scores.txt", "w") as fichier:
            for element in scores:
                fichier.write(str(element))
                fichier.write(' ')
            print('fichier scores.txt produit.')  


     
    