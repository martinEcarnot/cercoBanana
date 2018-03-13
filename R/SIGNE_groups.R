library(MASS)
# library(mixOmics)
# library(FactoMineR)
library(signal)
library(plyr)
library(caret)

# rm(list = ls())

source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les r?p?titions inutiles)
# set.seed(1)
#
# brb="~/Documents/INRA/Projets/SIGNE/2017/globalmatrix"
# load(file=brb)
#
# # Data Filter
# # Select dates
# dates=list(
#   "0524P"
#   ,"0529P"
#   ,"0606P"
#   # ,"0612P"
#   # ,"0619P"
#   # ,"0626P"
#   # "0703P"
#   # ,"0710P"
#   # ,"0717P"
#   # ,"0724P"
#   # ,"0731P"
#   # ,"0529S"
#   # ,"0606S"
#   # ,"0612S"
#   # ,"0619S"
#   # ,"0626S"
#   # ,"0703S"
#   # ,"0710S"
#   # ,"0717S"
#   # ,"0724S"
#   # ,"0731S"
# )
# iok=substr(rownames(res0)%in% dates  ## To keep only the date(s) you want
# res0=res0[iok,]
#
# titre=rownames(res0)
# ## decommenter pour retirer les jeunes feuilles:
# w=as.numeric(substr(titre,11,12))==4 | as.numeric(substr(titre,11,12))==5 | as.numeric(substr(titre,11,12))==6 | as.numeric(substr(titre,11,12))==10 | as.numeric(substr(titre,11,12))==11 | as.numeric(substr(titre,11,12))==12 | as.numeric(substr(titre,11,12))==16 | as.numeric(substr(titre,11,12))==17 | as.numeric(substr(titre,11,12))==18
# # ## decommenter pour retirer les vieilles feuilles:
# # # w=as.numeric(substr(titre,5,6))==1 | as.numeric(substr(titre,5,6))==2 | as.numeric(substr(titre,5,6))==3 | as.numeric(substr(titre,5,6))==7 | as.numeric(substr(titre,5,6))==8 | as.numeric(substr(titre,5,6))==9 | as.numeric(substr(titre,5,6))==13 | as.numeric(substr(titre,5,6))==14 | as.numeric(substr(titre,5,6))==15
# # ## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix[(complete.cases(w==TRUE)),]
# titre=rownames(globalmatrix)
# ## decommenter pour retirer la syrah:
# z=as.numeric(substr(titre,1,3))==471 | as.numeric(substr(titre,1,3))==525 | as.numeric(substr(titre,1,3))==747 | as.numeric(substr(titre,1,3))==877
## decommenter pour retirer le cabernet sauvignon:
# z=as.numeric(substr(titre,1,3))==015 | as.numeric(substr(titre,1,3))==169 | as.numeric(substr(titre,1,3))==685
## decommenter pour retirer le gamay:
# z=as.numeric(substr(titre,1,3))==787 | as.numeric(substr(titre,1,3))==509 | as.numeric(substr(titre,1,3))==222
## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix[(z==FALSE),]

## FIXATION DES PARAMETRES UTILISES:
# nombre de repetitions de la boucle de FDA:
repet=4
#Parametres du Savitsky-Golay (p=degr? du polynome, n= taille de la fen?tre, m=ordre de d?rivation)
p=2
n=11
m=2
# nombre de DV max autoisees
ncmax=20

#Taille de l'echantillon de validation (1/v):
# v=3
# Nombre de groupes de CV
k=6

sp=resok
ns=nrow(sp)

# creation de la matrice de classes
class=as.factor(substr(idtotok,17,19))
# variable qui mesure le nombre de classes
c=length(levels(class))

## Pretraitements
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
# sp=adj_asd(sp,c(651,1451))
# # # Reduction des variables (extremites bruitees)
# sp=sp[,seq(100,2100,1)]
# # # SNV
# sp=t(scale(t(sp)))
# # # Derivation Savitsky Golay
# sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

## Boucle pour effectuer plusieurs PLSDA (reduire impact tirage aleatoire)

# initialisation vecteur de % de bons classments par DV
perok=vector(mode='numeric',length=ncmax)
# creation matrice de % de mauvais classements par clone
mc=matrix(nrow = ncmax,ncol = c)

## definition des matrices de resultat final
# creation matrice des perok finale
perok_finalm0=matrix(nrow = repet, ncol = ncmax)

# creation de la matrice des dv et perok maximaux
maxi_final=matrix(nrow= repet, ncol = 2)
# creation de la matrice de % de mauvais classements
mc_final=matrix(nrow= repet, ncol = length(levels(class)))
# creation d'un matrice cubique pour enregistrer les tables de contingence
t_final=array(dim=c(c,c,repet))
# noms des colonnes et lignes
colnames(t_final)=c(basename(levels(class)))
rownames(t_final)=c(basename(levels(class)))
colnames(maxi_final)= c("maxi.id","perok max")
colnames(mc_final)= c(basename(levels(class)))

for(j in 1:repet) {

  # creation des jeux d'apprentissage et validation
  flds <- createFolds(1:ns, k = k)
  predm0=as.data.frame(matrix(nrow = ns, ncol = ncmax))

  # Boucle sur les groupes de CV
  for (i in 1:k) {
    id_val=sort(unlist(flds[i]))
    sp_val=sp[id_val,]
    class_val=class[id_val]
    sp_cal=sp[-id_val,]
    class_cal=class[-id_val]

    # PLSDA and application to have loadings and scores
    rplsda=caret::plsda(sp_cal, class_cal,ncomp=ncmax)
    sc_cal=rplsda$scores
    sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
    sc_val=sp_val_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    for (ii in 2:ncmax) {
      # Validation
      predm0[id_val,ii]=SIGNE_maha0(sc_cal[,1:ii], class_cal, sc_val[,1:ii])$class
    }
  }
    # Table de contingence
    tsm0=lapply(as.list(predm0), class, FUN = table)

    # Matrice mauvais classements par clone
    # mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
    # Sumpred=lapply(ts, FUN = rowSums)
    diagsm0=lapply(tsm0, FUN = diag)

    # Poucentage de bien classes
    perokm0=100*unlist(lapply(diagsm0, FUN = sum))/length(class)

    # pred[id_val,]=predict(rplsda ,sp_val,ncomp = 1:ncmax)

  # # Table de contingence
  # ts=lapply(as.list(pred), class, FUN = table)
  # # Matrice mauvais classements par clone
  # # mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
  # # Sumpred=lapply(ts, FUN = rowSums)
  # diags=lapply(ts, FUN = diag)
  # # Poucentage de bien classes
  # perok=100*unlist(lapply(diags, FUN = sum))/length(class)

  maxi=max(perokm0)
  maxi.id=which.max(perokm0)

  ## Enregistrement des matrices de resultat final
  # remplissage de la matrice des perok finale

  perok_finalm0[j,]=perokm0

  # remplissage de la dv max et de son % de bon classements globaux
  maxi_final[j,1]=maxi.id
  maxi_final[j,2]=maxi
  # remplissage de la matrice de mauvais classements par clone
  # mc_final[j,]=mc[maxi.id,]
  # t_final[,,j]=ts
}

# affichage matrices de resultat final
cat("RESULTATS FINAUX SUR", repet, "TIRAGES", "\n")
# print("Ensemble des tables de contingence:")
# print(t_final)
# print("Mauvais classements par clone:")
# print(mc_final)
# print("moyenne de mauvais classements par clone:")
# mc_final=mc_final[complete.cases(mc_final),]
# print(colMeans(mc_final))
# print(perok_final)
# print("Nombre de DV au max et maximum de bon classements:")
# print("Pourcentages d'erreur:")
# print(maxi_final)
# # Tracage des spectres
# # matplot(t(sp),pch='.')
# # Tracage de l'evolution des perok en fct du nb de DV utilisees
plot(colMeans(perok_finalm0))

# # # Tracage des moyennes de mauvais classements par clone
# # plot(colMeans(mc_final))
print("Moyenne de DV max et de perok:")
print(colMeans(maxi_final))

mc_final[!rowSums(!is.finite(mc_final)),]
mc_final[!is.finite(mc_final)] <- 0

mc_cep=vector(mode="logical", length = 3)
mc_cep[1]=mean(colMeans(mc_final)[1],colMeans(mc_final)[2],colMeans(mc_final)[7])
mc_cep[2]=mean(colMeans(mc_final)[3],colMeans(mc_final)[5],colMeans(mc_final)[9])
mc_cep[3]=mean(colMeans(mc_final)[4],colMeans(mc_final)[6],colMeans(mc_final)[8],colMeans(mc_final)[10])
names(mc_cep)=c("cabernet sauv","gamay","syrah")
# print("pourcentage d'erreur moyenne de classement des clones par c?page:")
# print(mc_cep)
############ END ############

## Global model
rplsdag=caret::plsda(sp, class,ncomp=ncmax)
scg=rplsdag$scores

## PLot
# plot(sc_cal[,1],sc_cal[,2], col=class_cal,pch=substr(rownames(sp_cal),3,3))
# plot(sc_cal[,1],sc_cal[,2], pch='')
# plot(sc_cal[,1],sc_cal[,2], col=dat_cal)
# text(sc_val[,1],sc_val[,2], labels=substr(rownames(sp_val),11,12), col = "red")
# plot(sc_cal[,1],sc_cal[,2], col=as.factor(v))

# plot(scg[,1],scg[,2], col=class)
# points(sc2[,1],sc2[,2],pch=3)
