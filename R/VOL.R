#加载包
library(ggpubr)
library(ggplot2)
#读取数据
Dat <- read.csv("C:/Users/runru/Desktop/seq/VOL/all.FDR.csv",header = T)
#作图
ggplot(Dat,aes(x=log2FoldChange,y=-log10(FDR)))+ geom_point()
#加载包
library(ggrepel)
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$FDR < 0.05 & abs(Dat$log2FoldChange) >= 1.5, ifelse(Dat$log2FoldChange>= 1.5 ,'Up','Down'),'No'),levels=c('Up','Down','No'))
#新加一列label
#Dat$Label=""

#对差异表达基因的FDR值进行从小到大排top 20

#Dat<-Dat[order(Dat$FDR)]
# 高表达的基因中，选择FDR最小的30
#Dat <- Dat[order(abs(Dat$FDR),decreasing = T), ]
#fdr.genes <- head(Dat$锘縄D, 30)
#DEG.top30<-c(as.character(fdr.genes))
#Dat$Label[match(DEG.top30,Dat$锘縄D)]<-DEG.top30

p<-ggplot(Dat,aes(x=log2FoldChange,y=-log10(FDR),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#BC3C28","MidnightBlue","grey"))+#确定点的颜色
geom_text_repel(
    data = Dat,
    label = Dat$Label,
    size = 3,
    segment.color = "black", show.legend = FALSE,family="sans" )+#添加关注的点的基因名  

  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-Log10(P)')+#修改y轴名称
  xlab('Log2(FoldChange)')+#修改x轴名称
  ggtitle("Volcano Plot of DEGs")+#标题
  geom_vline(xintercept=c(-1.5,1.5),lty=3,col="black",lwd=0.5) +#添加横线log|FoldChange|>1.5
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p<-p + xlim(-8,8) +ylim(0,25)#设置横纵坐标范围
 
p<-p+theme(axis.text.x =element_text(family="sans",size=12,color="black",face="bold"), axis.text.y=element_text(family="sans",size=12,color="black",face="bold"))
p<-p+theme(title =element_text(family="sans",size=12,color="black",face="bold",hjust=0.5))
p<-p+theme(legend.text = element_text(family="sans",face="bold",size = 12))#图例内容
p<-p+theme(plot.title = element_text(family="sans",face="bold",size = 14,hjust=0.5))#标题居中
p<-p  + theme(panel.grid=element_blank())#除去背景网格
print(p)
ggsave("VOL-M-sans.png",width=6,height=5.5)# 设定画布大小


