setwd('H:/causal/radar/')

# plot radar
library(fmsb)
tmp <- read.csv('D3_Insilico50_hetero_AUROC.csv',row.names = 1)
tmp=rbind(rep(0.7,5),rep(0.5,5),tmp)
radarchart(tmp,axistype = 4,plwd = 4,plty = 1,caxislabels = c(0.5,0.55,0.60,0.65,0.70),
           pcol = c("#FC4E07","#00AFBB", "#E7B800"),cglcol = 'grey',cglty = 1,
           axislabcol = 'black',cglwd = 1,vlcex = 1,title = 'DREAM3_Insilico100_trajectory')
legend(x = 'right',title = 'AUPR',col = c( "#FC4E07","#00AFBB", "#E7B800"),pt.cex=0,legend = rownames(tmp[-c(1,2),]),lty = 1,bty = 'n',cex = 1.2,lwd =3)



tmp <- read.csv('d3_aupr_heatmap_data.csv',row.names = 1)
library(pheatmap)
library(RColorBrewer)
pheatmap(tmp,scale = 'none',cluster_cols = F,cluster_rows = F,
         color = brewer.pal(9,'OrRd'),cellheight = 16,cellwidth = 12,main = 'AUPR')

tmp <- read.table('tmp.txt',header = T)
library(ggplot2)

ggplot(tmp,aes(AUROC,AUPR))+
  geom_point(aes(shape = DataSet,color = Method),size=5)+
  #scale_shape_manual(values = c())+
  scale_color_manual(values = c("#377EB8","#E41A1C" , "#4DAF4A"))+
  theme_bw()+
  scale_shape_manual(values = c(15,16,17,18,25))+
  scale_fill_manual(values = c(NA, "red"))





