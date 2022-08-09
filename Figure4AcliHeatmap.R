######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)
cliFile="clinical.txt"      #临床数据文件
riskFile="risk.txt"         #风险文件
setwd("D:\\biowolf\\FerrLnc\\24.cliHeatmap")

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65, ">65","<=65") )

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
cli=cli[samSample,,drop=F]
risk=risk[samSample,,drop=F]
data=cbind(risk, cli)
data=data[order(data$riskScore),,drop=F]      #根据风险打分对样品排序
Type=data[,(ncol(risk):ncol(data))]           #提取临床信息，作为热图注释文件
exp=data[,(3:(ncol(risk)-2))]                 #提取lncRNA表达量

#对临床性状进行循环，得到显著性标记
sigVec=c("risk")
for(clinical in colnames(Type[,2:ncol(Type)])){
	data=Type[c("risk", clinical)]
	colnames(data)=c("risk", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	tableStat=table(data)
	stat=chisq.test(tableStat)
	pvalue=stat$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(clinical, Sig))
	#print(paste(clinical, pvalue, Sig, sep="\t"))
}
colnames(Type)=sigVec

#定义热图注释的颜色
colorList=list()
#Type=Type[apply(Type,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#FF0000","#0066FF","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", 
         "#6ad157", "#f7aa5d", "#9ed84e", "#39ba30", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#ddd53e", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502", "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
	cliLength=length(levels(factor(Type[,cli])))
	cliCol=bioCol[(j+1):(j+cliLength)]
	j=j+cliLength
	names(cliCol)=levels(factor(Type[,cli]))
	if("unknow" %in% levels(factor(Type[,cli]))){
		cliCol["unknow"]="grey75"}
	colorList[[cli]]=cliCol
}

#热图可视化
pdf("heatmap.pdf", width=9, height=6)
pheatmap(t(exp),
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("green",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
