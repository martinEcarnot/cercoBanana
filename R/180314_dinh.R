library(MASS)
# library(mixOmics)
# library(FactoMineR)
library(signal)
library(plyr)
library(caret)

source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha0.R')

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d='D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 06 09 rge'
res0=SIGNE_load(d)




# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
res=adj_asd(res0,c(651,1451))
# # Reduction des variables (extremites bruitees)
res=res[,seq(100,2100,1)]
# Normalisation SNV
res=t(scale(t(res)))

# nombre de repetitions de la boucle de FDA:
repet=4

# # Derivation Savitsky Golay
#Parametres du Savitsky-Golay (p=degr? du polynome, n= taille de la fen?tre, m=ordre de d?rivation)
p=2
n=11
m=2
res=t(apply(res,1,sgolayfilt,p=p,n=n,m=m))
ncmax=20
# Nombre de groupes de CV
k=6

nboit=rownames(res)
idtot=vector()

for (i in 1:nrow(res))  {
  b=substr(nboit[i],1,nchar(nboit[i])-5)
  ib=which(b==id$Pastille.boite)

  id1=sprintf("%s_%9s_%s",substr(basename(d),9,18),id$Variete[ib],id$Boite.ds.var[ib])
  idtot=rbind(idtot,id1)
}

iok = which(substr(idtot,17,19) == "T2." |substr(idtot,17,19) == "T3."|substr(idtot,17,19) == "T4."|substr(idtot,17,19) == "T5."|substr(idtot,17,19) == "T6."|substr(idtot,17,19) == "T7."|substr(idtot,17,19) == "T8."|substr(idtot,17,19) == "T9.")

resok=res[iok,]
idtotok=idtot[iok]

# creation de la matrice de classes
class=as.factor(substr(idtotok,17,19))
# variable qui mesure le nombre de classes
c=length(levels(class))




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

  #maxi=max(perokm0)
  #maxi.id=which.max(perokm0)

  ## Enregistrement des matrices de resultat final
  # remplissage de la matrice des perok finale

  #perok_finalm0[j,]=perokm0

  # remplissage de la dv max et de son % de bon classements globaux
  #maxi_final[j,1]=maxi.id
  #maxi_final[j,2]=maxi
  # remplissage de la matrice de mauvais classements par clone
  # mc_final[j,]=mc[maxi.id,]
  # t_final[,,j]=ts
}

# affichage matrices de resultat final
#cat("RESULTATS FINAUX SUR", repet, "TIRAGES", "\n")

# # Tracage des spectres
# # matplot(t(sp),pch='.')
# # Tracage de l'evolution des perok en fct du nb de DV utilisees
#plot(colMeans(perok_finalm0))

# # # Tracage des moyennes de mauvais classements par clone
# # plot(colMeans(mc_final))
#print("Moyenne de DV max et de perok:")
#print(colMeans(maxi_final))

#mc_final[!rowSums(!is.finite(mc_final)),]
#mc_final[!is.finite(mc_final)] <- 0

#mc_cep=vector(mode="logical", length = 3)
#mc_cep[1]=mean(colMeans(mc_final)[1],colMeans(mc_final)[2],colMeans(mc_final)[7])
#mc_cep[2]=mean(colMeans(mc_final)[3],colMeans(mc_final)[5],colMeans(mc_final)[9])
#mc_cep[3]=mean(colMeans(mc_final)[4],colMeans(mc_final)[6],colMeans(mc_final)[8],colMeans(mc_final)[10])
#names(mc_cep)=c("cabernet sauv","gamay","syrah")
# print("pourcentage d'erreur moyenne de classement des clones par c?page:")
# print(mc_cep)
