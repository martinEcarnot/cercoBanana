library(ggplot2)
library(signal)

source('D:\\Intern 2018\\Dinh\\scripts R\\SIGNE_load.R')  
source('D:\\Intern 2018\\Dinh\\scripts R\\adj_asd.R')  

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d="D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 06 21"
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
  ib=which(b==id$Pastille.bo�te)
  
  id1=sprintf("%s_%9s_%s",substr(basename(d),9,18),id$Vari�t�[ib],id$Bo�te.ds.var[ib])
  idtot=rbind(idtot,id1)
}

iok = which(substr(idtot,17,19) == "T2." |substr(idtot,17,19) == "T3."|substr(idtot,17,19) == "T4."|substr(idtot,17,19) == "T5."|substr(idtot,17,19) == "T6."|substr(idtot,17,19) == "T7."|substr(idtot,17,19) == "T8."|substr(idtot,17,19) == "T9.")

resok=res[iok,]
geno=substr(idtot[iok],17,19)
rpca=prcomp(resok)

#iout=[]
#if (substr(basename(d),9,18)=="2017 06 21") {iout=which(rpca$x < -0.5)}
#if (substr(basename(d),9,18)=="2017 07 06") {iout=which(rpca$x < -0.3)}
iout=which(rpca$x < -0.5)
resok=resok[-iout,]
geno=geno[-iout]
rpca=prcomp(resok)

vectorgeno = as.factor(geno)


df=data.frame(rpcax=rpca$x,vectorgeno=vectorgeno)
p=ggplot(data = df, aes(rpcax.PC1,rpcax.PC2, colour=vectorgeno)) + geom_point()
plot(p)

