install.packages("ggplot2")

install.packages("tikzDevice")
library(tikzDevice)

library(ggplot2)


###############################################
#Evolution of Lambda Values
###############################################
#Define relevant dates
dates <- seq(as.Date("2001-07-01"), as.Date("2021-06-30"), by="months")

###############################################
#Lasso Paramters
###############################################

###############################################
#Standard Lasso
###############################################

tikz(file = "Lambda_Lasso.tex", width = 6, height = 2)

plot <- ggplot(data.frame(dates[-c(1:13)],Lasso$Evaluation$Lambda),aes(x=dates[-c(1:13)],y=Lasso$Evaluation$Lambda) ) +
          geom_segment(aes(x=dates[-c(1:13)], xend=dates[-c(1:13)], y=0, yend=Lasso$Evaluation$Lambda), size = 0.5, linetype = "solid", color = "grey") +
          geom_smooth(se = F) +
          theme_bw() +
          scale_x_date(date_breaks = "1 year", date_labels = "%Y",expand = c(0,0)) +
          theme(axis.text.x = element_text(size=rel(0.7),angle = 90, vjust=0.5),
                axis.title.x = element_text(size=rel(0.8)),
                axis.title.y = element_text(size=rel(0.8))) +  
          geom_point(shape = 21, size = 1, color="black",fill = "orange") +
          ylab("Optimal $\\lambda$") +
          xlab("Time")
#Preview of plot
print(plot)


dev.off()

###############################################
#peLasso (Average)
###############################################
tikz(file = "Lambda_peAverage.tex", width = 6, height = 2)

plot <- ggplot(data.frame(dates[-c(1:13)],peLassoAverage$Evaluation$Lambda),aes(x=dates[-c(1:13)],y=peLassoAverage$Evaluation$Lambda) ) +
  geom_segment(aes(x=dates[-c(1:13)], xend=dates[-c(1:13)], y=0, yend=peLassoAverage$Evaluation$Lambda), size = 0.5, linetype = "solid", color = "grey") +
  geom_smooth(se = F) +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",expand = c(0,0)) +
  theme(axis.text.x = element_text(size=rel(0.7),angle = 90, vjust=0.5),
        axis.title.x = element_text(size=rel(0.8)),
        axis.title.y = element_text(size=rel(0.8))) +  
  geom_point(shape = 21, size = 1, color="black",fill = "orange") +
  ylab("Optimal $\\lambda$") +
  xlab("Time")

print(plot)

dev.off()

###############################################
#peLasso (Ridge)
###############################################
tikz(file = "Lambda_peRidgeStep1.tex", width = 6, height = 2)

plot <- ggplot(data.frame(dates[-c(1:13)],peLassoRidge$Evaluation$Lambda..Lasso.),aes(x=dates[-c(1:13)],y=peLassoRidge$Evaluation$Lambda..Lasso.) ) +
  geom_segment(aes(x=dates[-c(1:13)], xend=dates[-c(1:13)], y=0, yend=peLassoRidge$Evaluation$Lambda..Lasso.), size = 0.5, linetype = "solid", color = "grey") +
  geom_smooth(se = F) +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",expand = c(0,0)) +
  theme(axis.text.x = element_text(size=rel(0.7),angle = 90, vjust=0.5),
        axis.title.x = element_text(size=rel(0.8)),
        axis.title.y = element_text(size=rel(0.8))) +  
  geom_point(shape = 21, size = 1, color="black",fill = "orange") +
  ylab("Optimal $\\lambda$") +
  xlab("Time")

print(plot)

dev.off



###############################################
#Plots of Ridge Parameters
###############################################

###############################################
#Standard Ridge
###############################################
tikz(file = "Lambda_Ridge.tex", width = 6, height = 2)

plot <- ggplot(data.frame(dates[-c(1:13)],Ridge$Evaluation$Lambda),aes(x=dates[-c(1:13)],y=Ridge$Evaluation$Lambda) ) +
  geom_segment(aes(x=dates[-c(1:13)], xend=dates[-c(1:13)], y=0, yend=Ridge$Evaluation$Lambda), size = 0.5, linetype = "solid", color = "grey") +
  geom_smooth(se = F) +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",expand = c(0,0)) +
  theme(axis.text.x = element_text(size=rel(0.7),angle = 90, vjust=0.5),
        axis.title.x = element_text(size=rel(0.8)),
        axis.title.y = element_text(size=rel(0.8))) +  
  geom_point(shape = 21, size = 1, color="black",fill = "orange") +
  ylab("Optimal $\\lambda$") +
  xlab("Time")

print(plot)

dev.off()

###############################################
#peLasso (Ridge)
###############################################
tikz(file = "Lambda_peRidgeStep2.tex", width = 6, height = 2)

plot <- ggplot(data.frame(dates[-c(1:13)],peLassoRidge$Evaluation$Lambda..Ridge.),aes(x=dates[-c(1:13)],y=peLassoRidge$Evaluation$Lambda..Ridge.) ) +
  geom_segment(aes(x=dates[-c(1:13)], xend=dates[-c(1:13)], y=0, yend=peLassoRidge$Evaluation$Lambda..Ridge.), size = 0.5, linetype = "solid", color = "grey") +
  geom_smooth(se = F) +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",expand = c(0,0)) +
  theme(axis.text.x = element_text(size=rel(0.7),angle = 90, vjust=0.5),
        axis.title.x = element_text(size=rel(0.8)),
        axis.title.y = element_text(size=rel(0.8))) +  
  geom_point(shape = 21, size = 1, color="black",fill = "orange") +
  ylab("Optimal $\\lambda$") +
  xlab("Time")

print(plot)

dev.off()

###############################################
#Plot of Barrier Function
###############################################
#Logarithmic funcitons for different values of mu
eq1 = function(x){ -0.5 * log(-x)}
eq2 = function(x){ -1 * log(-x)}
eq3 = function(x){ -2 * log(-x)}

#Plot
tikz(file = "Barrier.tex", width = 4.5, height = 3.25)

plot <- ggplot(data.frame(x=c(-5, -0.0001)), aes(x=x)) +
          geom_path(aes(colour="0.5"), stat="function", fun=eq1) +
          geom_path(aes(colour="1"), stat="function", fun=eq2) +
          geom_path(aes(colour="2"), stat="function", fun=eq3) +
          geom_segment(aes(x=0, y=0 , xend = 0,  yend = 5), linetype = "dashed", colour = "#666666") +  
          geom_segment(aes(x=0, y=0 , xend = -5,  yend = 0), linetype = "dashed", colour = "#666666")  +
          labs(color="$\\mu$") +
          ylab("- $\\mu$ log(x)") +
          xlab("x") +
          theme_bw() +
          coord_cartesian(
            xlim = c(-5,2.5),
            ylim = c(-2.5,5),
            expand = FALSE,
            default = FALSE,
            clip = "on")
  
print(plot)

dev.off()
