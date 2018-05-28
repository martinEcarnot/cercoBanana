library(ggplot2)
library(signal)
library(hyperSpec)

source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d="D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 06 09 rge"
res0=SIGNE_load(d)

# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
res=adj_asd(res0,c(651,1451))
# # Reduction des variables (extremites bruitees)
res=res[,seq(100,2100,1)]
# Normalisation SNV
res=t(scale(t(res)))
# # Derivation Savitsky Golay
#Parametres du Savitsky-Golay (p=degr? du polynome, n= taille de la fen?tre, m=ordre de d?rivation)
# p=2
# n=11
# m=2
# res=t(apply(res,1,sgolayfilt,p=p,n=n,m=m))


nboit=rownames(res)
idtot=vector()

for (i in 1:nrow(res))  {
  b=substr(nboit[i],1,nchar(nboit[i])-5)
  ib=which(b==id$Pastille.boite)


  id1=sprintf("%s_%9s_%s",substr(basename(d),9,18),id$Variete[ib],id$Boite.ds.var[ib])
  if (length(id1)==0) {id1=nboit[i]}
  idtot=rbind(idtot,id1)
}

iok = which(substr(idtot,17,19) == "T2."|substr(idtot,17,19) == "T3."|substr(idtot,17,19) == "T4."|substr(idtot,17,19) == "T5."|substr(idtot,17,19) == "T6."|substr(idtot,17,19) == "T7."|substr(idtot,17,19) == "T8."|substr(idtot,17,19) == "T9.")

resok=res[iok,]
idtotok=idtot[iok]

#To do the average of the 3 spectra
# resok=aggregate(resok,list(idtotok),mean)
# idtotok=resok$Group.1
# resok=as.matrix(resok[,2:ncol(resok)])


geno=substr(idtotok,17,19)
# geno[geno %in% c("T2.2_","T2.3_","T2.4_")]="T2"
# geno[geno %in% c("T3.1","T3.3","T3.4")]= "T3"
# geno[geno %in% c("T4.2","T4.3","T4.4")] = "T4"
# geno[geno %in% c("T5.1","T5.2","T5.3")] = "T5"
# geno[geno %in% c("T6.1","T6.2","T6.4")] = "T6"
# geno[geno %in% c("T7.1","T7.3","T7.4")] = "T7"
# geno[geno %in% c("T8.3","T8.4","T8.5")] = "T8"
# geno[geno %in% c("T9.1","T9.2","T9.5")] = "T9"


#combine levels geno
geno = as.character(geno)
#
vectorgeno = as.factor(geno)

rownames(resok)=idtotok
write.table(resok, file=paste(trimws(basename(d)),"_pre_sans der_san moyen.csv"), sep=";", dec = ".",row.names = idtotok)
write.table(as.numeric(vectorgeno), file=paste(trimws(basename(d)),"_class_sans der_san moyen.csv"),sep=";", dec = ".",row.names = idtotok)

