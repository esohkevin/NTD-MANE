#!/usr/bin/Rscript

args <- commandArgs(TRUE)
pe <- args[1]
te <- args[2]
pr <- args[3]
de <- args[4]
ra <- args[5]
rb <- args[6]
ga <- args[7]
gb <- args[8]

engy_plot <- function(pe=potential_file,
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

ansis_plot <- function(rmsda=first_rmsd_file,
		       rmsdb=second_rmsd_file, 
		       gyra=first_gyration_file,
		       gyrb=second_gyration_file) {

   ra <- read.table(rmsda, h=F)
   colnames(ra) <- c("time","rmsd")

   rb <- read.table(rmsdb, h=F)
   colnames(rb) <- c("time","rmsd")

   ga <- read.table(gyra, h=F)
   colnames(ga) <- c("time","Rt","Rx","Ry","Rz")
   ga$t_ns <- ga$time/1000

   gb <- read.table(gyrb, h=F)
   colnames(gb) <- c("time","Rt","Rx","Ry","Rz")
   gb$t_ns <- gb$time/1000


   png("grom_analysis.png", height=10, width=6, units = "in", points=10, res=200)
   par(mfrow=c(2,1))
      plot(ra$time, ra$rmsd, 
           main="RMSD",
      	type="l",
      	xlab="Time (ns)",
      	ylab="RMSD (nm)")
      lines(rb$time, rb$rmsd,col="red")
      legend("bottomright", legend=c("Mutant","Wild Type"), col=c("black","red"),pch="-")

      plot(ga$t_ns, ga$Rt,
           main="Radius of Gyration",
           type="l",
           xlab="Time (ns)",
           ylab="Rg (nm)",
	   ylim=c(1.30, 1.50))
      lines(gb$t_ns, gb$Rt,col="red")
   dev.off()
}


engy_plot(pe=pe,te=te,pr=pr,de=de)
ansis_plot(rmsda=ra, rmsdb=rb, gyra=ga, gyrb=gb)
