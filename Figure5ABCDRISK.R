#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")

#引用包
library(survival)
library(survminer)
library(timeROC)

inputFile="risk.txt"        #输入文件
survFile="survival.pdf"         #生存曲线文件
rocFile="ROC.pdf"               #ROC曲线文件
setwd("C:\\Users\\11111\\Desktop\\坏死性凋亡necroptosis\\17生存曲线及ROC曲线")      #设置工作目录


#读取输入文件
rt=read.table(inputFile,header=T,sep="\t")

#比较高低风险组生存差异，得到显著性p值
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
#绘制生存曲线
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=5,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()


###ROC曲线
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=5,height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
        c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
          paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
          paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
        col=c("green",'blue','red'),lwd=2,bty = 'n')
dev.off()

