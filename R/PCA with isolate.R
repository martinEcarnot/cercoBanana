library(ggplot2)
library(signal)

source('D:\\Intern 2018\\cercoBanana\\R\\SIGNE_load.R')
source('D:\\Intern 2018\\cercoBanana\\R\\adj_asd.R')

id=read.table("D:\\Intern 2018\\Dinh\\Classeur fs5.csv",sep=",",header =TRUE)

d="D:\\Intern 2018\\Dinh\\biov tous asd\\biov1 fs5 asd\\fs5 asd 2017 07 18"
res0=SIGNE_load(d)

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

isolate=substr(idtotok,23,24)
rpca=prcomp(resok)


rpca=prcomp(resok)

vectorisolate = as.factor(isolate)


df=data.frame(rpcax=rpca$x,vectorisolate=vectorisolate)
p=ggplot(data = df, aes(rpcax.PC1,rpcax.PC2, colour=vectortisolate)) + geom_point()
plot(p)


