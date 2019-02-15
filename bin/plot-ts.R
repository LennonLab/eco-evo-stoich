#Plot function for each chemostat community
plot.ts <- function(data, cid = "", day.start, day.end){
  par(mar = c(5,6,5,1))
  abd <- data[data$cID == cid,]
  pt <- plot(abd$day[abd$microbe == "Syn"],abd$abd[abd$microbe == "Syn"],
             xlim = c(day.start,day.end),ylim = c(1*10^4,10^9), log = "y", xaxt = 'n',yaxt = 'n',
             xlab = expression(bold("Time (d)")), ylab = expression(bold(paste("Abundance (mL"^"-1",")"))), cex.lab = 1.5, main = cid,
             lty = 1, lwd = 2, type = "l", col = "black", 
             pch = 21, cex = 2.5)
  
  axis(1, cex.axis = 1.25,c(-120,-80,-40,0,40,80,120,160))
  ticks <- seq(4, 9, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(2, cex.axis = 1.25, at = c(10^4,10^5, 10^6, 10^7, 10^8, 10^9), labels = labels, las = 1)
  points(abd$day[abd$microbe == "Phage"],abd$abd[abd$microbe == "Phage"],
         col = "green",type = 'l',pch = 21,lwd = 2,cex = 2.5,lty = 1)
  abline(v = 0, col = "grey", lty = 1, lwd = 1)
  box(lwd=1)
  return(abd)
}