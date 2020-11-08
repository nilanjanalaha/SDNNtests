#_____________________________________________________________________#

###################### ecdf of d1 and d2 ####################################
#' Descriptive statistics
#' @param d1 first dataset
#' @param d2 second dataset
#' @return Makes plots of the ECDF, the Boxplot, the KDE, and the
#'    histograms of the datasets d1 and d2
#' @export
do_description <- function(d1, d2)
{
# Pre-processing the data
Trial <- as.factor(c(rep("HVTN 097",length(d1)), rep("HVTN 100",length(d2))))
dat.boxplot <- data.frame(c(d1,d2),Trial)
colnames(dat.boxplot) <- c("Immune_response", "Trial")
dat.boxplot[,1] <- c(d1,d2)

#------------- Common width and height of hist, boxplot, and ecdf-------

wd <- 5
hi <- 3
sz1 <- 1.5
sz2 <- 2
al <- 0.5
ts <- 15 #  text size
tts <- 15 # axis label text size
tl <- 15 # legend text size

##################### Histogram  #############################################

p <- ggplot(dat.boxplot, aes(Immune_response,fill=Trial)) +geom_histogram(alpha=.5,position="identity")
p <- p +xlab("log(net MFI)")+theme_bw()+ylab("Frequency")
p <- p+ theme(legend.position = "top",
              axis.title= element_text(size= tts),
              axis.text.y=element_text(size= ts),
              axis.text.x=element_text(size= ts),
              legend.text=element_text(size= tl),
              strip.text.x=element_text(size= ts),
              legend.title=element_text(size= tts))
p <- p+scale_fill_manual(breaks = c("HVTN 097", "HVTN 100"),
                         values=c("red", "blue"))
print(p)
#ggsave("hist.pdf", p, width= wd, height= hi, units = "in")

#------------------ ECDF -------------------------------------------

p <- ggplot(dat.boxplot, aes(Immune_response,colour=Trial)) +stat_ecdf()
p <- p +xlab("log(net MFI)")+theme_bw()+ ylab(" ")
p <- p + theme(legend.position = "top",
               axis.title= element_text(size= tts),
               axis.text.y=element_text(size= ts),
               axis.text.x=element_text(size= ts),
               legend.text=element_text(size= tl),
               strip.text.x=element_text(size= ts),
               legend.title=element_text(size= tts))
p <- p+scale_color_manual(breaks = c("HVTN 097", "HVTN 100"),
                          values=c("red", "blue"))
print(p)
#ggsave("ecdf.pdf", p, width=wd, height=hi, units = "in")

##################### Boxplot of d1 and d2 ################################################

p <- ggplot(dat.boxplot, aes(x=Trial, y=Immune_response, color=Trial)) +geom_boxplot()
p <- p+theme_bw()+xlab("log(net MFI)")+ylab(" ")+ theme(legend.position = "top",
                                                        axis.title= element_text(size= tts),
                                                        axis.text.y=element_text(size= ts),
                                                        axis.text.x=element_text(size= ts),
                                                        legend.text=element_text(size= tl),
                                                        strip.text.x=element_text(size= ts),
                                                        legend.title=element_text(size= tts))
p <- p+scale_color_manual(breaks = c("HVTN 097", "HVTN 100"),
                          values=c("red", "blue"))
print(p)
#ggsave("box.pdf", p, width= wd, height= hi, units = "in")

################## KDE #####################################################
kd1 <-  kde(sort(d1), h=hlscv(d1))
kd2 <-  kde(sort(d2), h=hlscv(d2))




p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))+ylim(0, 0.25) + xlim(-2, 12)
p <- p + stat_function(fun = function(x)  dkde(x, kd1), aes(color="F"))
p <- p + stat_function(fun = function(x)  dkde(x, kd2), aes(color="G"))
p <- p+ theme_bw() + theme(legend.position = "top",
                           axis.title= element_text(size= tts),
                           axis.text.y=element_text(size= ts),
                           axis.text.x=element_text(size= ts),
                           legend.text=element_text(size= tl),
                           strip.text.x=element_text(size= ts),
                           legend.title=element_text(size= tts))
p <- p+labs(x= "log(net MFI)", y= " ")
#p <- p + theme( axis.text=element_blank())
p <- p +scale_color_manual(name = " ",
                           values = c( "red", "blue"), # Color specification
                           labels = c("HVTN 097", "HVTN 100"))
print(p)
#ggsave("kde.pdf", p, width= wd, height= hi, units = "in")
}
