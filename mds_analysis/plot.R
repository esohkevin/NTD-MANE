#!/usr/bin/Rscript

args <- commandArgs(TRUE)
p <- args[1]
t <- args[2]
p <- args[3]
d <- args[4]

mds_plot <- function(pe=potential_file,
		     te=temperature_file,
		     pr=pressure_file,
		     de=density_file) {
   
   e <- read.table(pe, h=F)
   colnames(e) <- c("Time_ps","Energy_kjmol_1")
   
   t <- read.table(te, h=F)
   colnames(t) <- c("Time_ps","K")
   t_avg <- mean(t$K)
   
   p <- read.table(pr, h=F)
   colnames(p) <- c("Time_ps","bar")
   #p$bar_10ps <- (p$bar)/10
   
   d <- read.table(de, h=F)
   colnames(d) <- c("Time_ps","den")
   #d$den_10ps <- (d$den)/10
   
   png("grom_energies.png", height=10, width=10, units = "in", points=10, res=200)
   par(mfrow=c(2,2))
   plot(e$Time_ps, e$Energy_kjmol_1, type="l", 
        main="Potential Energy: 1AKI, NVT Equilibration", 
        xlab="Time (ps)", 
        ylab="Potential Energy (kJ/mol)",
        col="red")
   #dev.off()
   
   plot(t$Time_ps, t$K, type="l", 
        main="Temperature: 1AKI, NVT Equilibration", 
        xlab="Time (ps)", 
        ylab="Temperature (K)",
        col = "black")
   abline(h=t_avg, lty=2)
   
   plot(p$Time_ps, p$bar, type="l",
        main="Pressure: 1AKI, NPT Equilibration",
        xlab="Time (ps)",
        ylab="Pressure (bar)",
        col="black")
   #lines(p$Time_ps, p$bar_10ps, col="red")
   
   plot(d$Time_ps, d$den, type="l",
        main="Density: 1AKI, NPT Equilibration",
        xlab="Time (ps)",
        ylab="Density (kg/m^3)",
        col="black")
   #lines(d$Time_ps, d$den_10ps, col="red")
   
   dev.off()

}

mds_plot(pe=p,te=t,pr=p,de=d)
