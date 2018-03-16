library(ggplot2)
library(signal)

source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d="D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 07 18"
res0=SIGNE_load(d)

# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
res=adj_asd(res0,c(651,1451))
# # Reduction des variables (extremites bruitees)
res=res[,seq(100,2100,1)]
# Normalisation SNV
res=t(scale(t(res)))
# # Derivation Savitsky Golay
#Parametres du Savitsky-Golay (p=degr? du polynome, n= taille de la fen?tre, m=ordre de d?rivation)
p=2
n=11
m=2
res=t(apply(res,1,sgolayfilt,p=p,n=n,m=m))


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
iso= substr(idtotok, 23, 24)
#geno=substr(idtotok,17,19)
rpca=prcomp(resok)

#iout=[]
#if (substr(basename(d),9,18)=="2017 06 21") {iout=which(rpca$x > 0.05)}
#if (substr(basename(d),9,18)=="2017 07 06") {iout=which(rpca$x > 0.2)}
if (substr(basename(d),9,18)=="2017 07 18") {iout=which(rpca$x < -0.02)}

resok=resok[-iout,]
iso=iso[-iout]
idtotok=idtotok[-iout]

rpca=prcomp(resok)

vectoriso = as.factor(iso)
vectoriso[which(vectoriso==".2")]=".1"
vectoriso[which(vectoriso==".4")]=".3"
vectoriso[which(vectoriso==".6")]=".5"

df=data.frame(rpcax=rpca$x,vectoriso=vectoriso)
p=ggplot(data = df, aes(rpcax.PC1,rpcax.PC2, colour=vectoriso)) + geom_point()
plot(p)


