library(ssmarina)
## expression matrix
exp_f = 'TPM_log_qn.PC.symbol.txt' 
## network file from aracne
net = 'network_0805.txt' 

df=read.table(exp_f,header = T, row.names = 1)  
df=data.matrix(df)
## 
day1_cols= grep("Int_1",colnames(df),value=T)  ## grep("Int_1|Int_8",colnames(df),value=T)
day8_cols= grep("Int_8",colnames(df),value=T) 
signature = rowTtest(df[, day8_cols], df[, day1_cols])
signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
signature[is.na(signature)] <- 0

nullmodel <- ttestNull(df[, day8_cols], df[, day1_cols], per=1000, repos=T)
regulon <- aracne2regulon(net, df)  
mrs <- marina(signature, regulon, nullmodel)
num_tf = length(mrs$regulon)
s_mrs=summary(mrs,num_tf)
## sort by pvalue
s_mrs=s_mrs[order(s_mrs[,4]),]
## rm NA
s_mrs=na.omit(s_mrs)
## Leading-edge analysis
mrs <- ledge(mrs)
## Shadow analysis
mrshadow <- shadow(mrs,)
tmp = summary(mrshadow,num_tf)
s_mrshadow = tmp$MARINA.results
s_mrshadow=s_mrshadow[order(s_mrshadow[,4]),]
s_mrshadow=na.omit(s_mrshadow)
write.table(s_mrshadow,'Int_mrshadow.txt',sep='\t',row.names =F)
pair_mrshadow = tmp$Shadow.pairs
pair = list()
for (i in 1:length(pair_mrshadow)){
pair[[i]] = pair_mrshadow[i]
}
write.table(pair,'Int_shadow_pair.txt',sep='\n',row.names=F,col.names=F)
## write 
tfs=as.character(s_mrs[['Regulon']])
n=length(tfs)
lst=list()
for (i in 1:n){
tf=tfs[i]
targets=names(mrs$regulon[[tf]]$'tfmode')
targets=paste(targets,collapse=",")
lst[[i]]=targets
}
s_mrs$targets=as.character(lst)
lst2=list()
for (i in 1:n){
tf=tfs[i]
targets=mrs$ledge[[tf]]
targets=paste(targets,collapse=",")
lst2[[i]]=targets
}
s_mrs$ledge=as.character(lst2)
write.table(s_mrs,'Int_mrs.txt',sep='\t',row.names =F)
pdf('Int_mrs.pdf')
plot(mrs,20)
dev.off()

