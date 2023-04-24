###################################################################################################
# Code to reproduce Figure 8 in paper
###################################################################################################


library(ggplot2)

set.seed(12)
X1 = c(0.3,0.7,1,1.4,1.5,2.1,2.3,1.5,1.4)
X2 = 6
Y1 = 1.5 + rnorm(length(X1),mean=0,sd=1)
Y2 = X2*2 + rnorm(1)
X = c(X1,X2)
Y = c(Y1,Y2)

Z = rnorm(length(X),mean=0,sd=1)
g_Y = Y+ Z
f_Y = Y -Z
split = c(rep("Split 1",length(X)/2),rep("Split 2",length(X)/2))


df <- data.frame(X = X, Y = Y, f_Y= f_Y, g_Y= g_Y,split = split[c(sample(1:length(X)))])

ggplot(subset(df,split == "Split 1"),aes(x=X,y=Y)) +
  geom_point(size = 8,color="blue",shape="circle") + 
  geom_smooth(method = "lm", fill = NA,fullrange=TRUE,linetype="dashed",size=2.5) +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position =c(0.2,0.80), legend.text = element_text(size = 30)) +
  ggtitle("Split 1") + 
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,15))
ggsave("example_highleverage_split1.pdf")

ggplot(subset(df,split == "Split 2"),aes(x=X,y=Y)) +
  geom_point(size = 8,color="red",shape="square") + 
  geom_smooth(method = "lm", color= "red",fill = NA,fullrange=TRUE,linetype="dashed",size=2.5) +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position =c(0.2,0.80), legend.text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,15)) +
  ggtitle("Split 2")
ggsave("example_highleverage_split2.pdf")

ggplot(df,aes(x=X,y=Y)) +
  geom_point(size = 8) + 
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position =c(0.2,0.80), legend.text = element_text(size = 30)) +
  scale_y_continuous(limits = c(0,15))
ggsave("example_highleverage_fullds.pdf")

ggplot(df,aes(x=X,y=Y,color=split,shape=split)) +
  geom_point(size = 8) +
  geom_smooth(method = "lm", fill = NA,fullrange=TRUE,linetype="dashed",size=2.5) +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position =c(0.2,0.80), legend.text = element_text(size = 30)) +
  scale_y_continuous(limits = c(0,15))
ggsave("example_highleverage_1.pdf")

ggplot(df,aes(x=X,y=f_Y)) +
  geom_point(color="blue",size=8) +
  geom_smooth(method = "lm", color = "blue", linetype ="dashed",fill = NA,fullrange=TRUE,size=2.5) +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position = "bottom", legend.text = element_text(size =30)) + 
  ylab("f(Y)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,15))
ggsave("example_highleverage_2.pdf")


ggplot(df,aes(x=X,y=g_Y)) +
  geom_point(color="darkgreen",size=8) +
  geom_smooth(method = "lm", color= "darkgreen", linetype="dashed", fill = NA,fullrange=TRUE,size=2.5) +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 30),
        legend.position = "bottom", legend.text = element_text(size = 30)) +
  ylab("g(Y)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,15))
ggsave("example_highleverage_3.pdf")


