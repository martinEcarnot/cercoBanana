library(MASS)
library(signal)
library(plyr)
library(caret)
library(ggplot2)
library(lattice)


source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha.R')
source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_maha0.R')

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d="D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 07 18"
res0=SIGNE_load(d)

## Pretraitements
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
res=adj_asd(res0,c(651,1451))
# # Reduction des variables (extremites bruitees)
res=res[,seq(100,2100,1)]
# Normalisation SNV
res=t(scale(t(res)))

nboit=rownames(res)
idtot=vector()

for (i in 1:nrow(res))
  {
  b=substr(nboit[i],1,nchar(nboit[i])-5)
  ib=which(b==id$Pastille.boite)

  id1=sprintf("%s_%9s_%s",substr(basename(d),9,18),id$Variete[ib],id$Boite.ds.var[ib])
  idtot=rbind(idtot,id1)
}

iok = which(substr(idtot,17,19) == "T2." |substr(idtot,17,19) == "T3."|substr(idtot,17,19) == "T4."|substr(idtot,17,19) == "T5."|substr(idtot,17,19) == "T6."|substr(idtot,17,19) == "T7."|substr(idtot,17,19) == "T8."|substr(idtot,17,19) == "T9.")

resok=res[iok,]
idtotok=idtot[iok]
resok=aggregate(resok,list(idtotok),mean)
idtotok=resok$Group.1
resok=as.matrix(resok[,2:ncol(resok)])

sp=resok
ns=nrow(sp)


## Boucle pour effectuer plusieurs PLSDA (reduire impact tirage aleatoire)
## FIXATION DES PARAMETRES UTILISES:
repet=4
# nombre de DV max autoisees
ncmax=20
# # Derivation Savitsky Golay
#Parametres du Savitsky-Golay
p=2
n=11
m=2
res=t(apply(res,1,sgolayfilt,p=p,n=n,m=m))
# Nombre de groupes de CV
k=6

# creation de la matrice de classes
iso= substr(idtotok, 23, 24)
vectoriso = as.factor(iso)
class=vectoriso
# variable qui mesure le nombre de classes
c=length(levels(class))

# initialisation vecteur de % de bons classments par DV
perok=vector(mode='numeric',length=ncmax)

## definition des matrices de resultat final
# creation matrice des perok finale
perok_finalm0=matrix(nrow = repet, ncol = ncmax)

for(j in 1:repet)
  {
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

  # Sumpred=lapply(ts, FUN = rowSums)
  diagsm0=lapply(tsm0, FUN = diag)

  # Poucentage de bien classes
  perokm0=100*unlist(lapply(diagsm0, FUN = sum))/length(class)

  maxi=max(perokm0)
  maxi.id=which.max(perokm0)

  ## Enregistrement des matrices de resultat final
  # remplissage de la matrice des perok finale

  perok_finalm0[j,]=perokm0
}

# # Tracage de l'evolution des perok en fct du nb de DV utilisees
plot(colMeans(perok_finalm0))


