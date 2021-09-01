library(bio3d)
# set file [open dipl_phase2_native_smd.log r]
# set output [open dipl_phase2_native_smd.dat w]
# set file [open dipl_phase2_native_smd.log r]
# while { [gets $file line] != -1 } {
#   ### Determine if a line contains SMD output. If so, write the
#   ### timestep followed by f(dot)n to the output file
#   if {[lindex $line 0] == "SMD"} {
#     puts $output "[lindex $line 1] [lindex $line 3] [lindex $line 6]]"
#   }
# }
# ### Close the log file and the output .dat file
# close $file
# close $output
# 
# set file [open tauro_phase2_native_smd.log r]
# set output [open tauro_phase2_native_smd.log.dat w]
# set file [open tauro_phase2_native_smd.log r]
# while { [gets $file line] != -1 } {
#   ### Determine if a line contains SMD output. If so, write the
#   ### timestep followed by f(dot)n to the output file
#   if {[lindex $line 0] == "SMD"} {
#     puts $output "[lindex $line 1] [lindex $line 4] [lindex $line 7]]"
#   }
# }
# ### Close the log file and the output .dat file
# close $file
# close $output
#diploptene dragged by y-axis
#tauro dragged in x-axis

#well-compiled log file till 1 ns
setwd('/Users/sayaneshome/Desktop/prelims_stuffs/data/data_analysis_dipl/')
f117r_dipl <- read.table('f117r_dipl_smd01.dat')
f541r_dipl <- read.table('f541r_dipl_smd01.dat')
l48f_dipl <- read.table('l48f_dipl_smd01.dat')
l826f_dipl <- read.table('l826f_Dipl_smd01.dat')
native_dipl <- read.table('smd01_native.dat')
#force-position graph 
plot(f117r_dipl$V2,f117r_dipl$V3,type = 'line')
lines(f541r_dipl$V2,f541r_dipl$V3,col = "red")
lines(l826f_dipl$V2,l826f_dipl$V3,col = "green")
lines(native_dipl$V2,native_dipl$V3,col = "purple")
#time-position graph;get the frames from this graph specially native one

plot((native_dipl$V1)/100000,smooth(native_dipl$V3),type = 'line',ylab="Force", xlab="Time(ns)",xlim = c(0,10))
lines((f117r_dipl$V1)/100000,smooth(f117r_dipl$V3),col = "purple")
lines((f541r_dipl$V1)/100000,smooth(f541r_dipl$V3),col = "red")
lines((l826f_dipl$V1)/100000,smooth(l826f_dipl$V3),col = "green")


#diploptene pdbs

f117r_dipl_pdb <- read.pdb('f117r_dipl_stride100_full.pdb',multi = TRUE)
native_dipl_pdb <- read.pdb('hpnn_dipl_stride100_full.pdb',multi = TRUE)
l826f_dipl_pdb <- read.pdb('l826f_dipl_stride100_full.pdb',multi = TRUE)
f541r_dipl_pdb <- read.pdb('f541r_dipl_stride100_full.pdb',multi = TRUE)


ca.inds_native <- atom.select(native_dipl_pdb, elety="CA")
xyz_native <- fit.xyz(fixed=native_dipl_pdb$xyz[1,], mobile=native_dipl_pdb$xyz[2:11,])
rd_native <- rmsd(xyz_native[1,ca.inds_native$xyz], xyz_native[,ca.inds_native$xyz])
ca.inds_l826f <- atom.select(l826f_dipl_pdb, elety="CA")
xyz_l826f <- fit.xyz(fixed=l826f_dipl_pdb$xyz[1,], mobile=l826f_dipl_pdb$xyz[2:11,])
rd_l826f <- rmsd(xyz_l826f[1,ca.inds_l826f$xyz], xyz_l826f[,ca.inds_l826f$xyz])
ca.inds_f117r <- atom.select(f117r_dipl_pdb, elety="CA")
xyz_f117r <- fit.xyz(fixed=f117r_dipl_pdb$xyz[1,], mobile=f117r_dipl_pdb$xyz[2:11,])
rd_f117r <- rmsd(xyz_f117r[1,ca.inds_f117r$xyz], xyz_f117r[,ca.inds_f117r$xyz])
ca.inds_f541r <- atom.select(f541r_dipl_pdb, elety="CA")
xyz_f541r <- fit.xyz(fixed=f541r_dipl_pdb$xyz[1,], mobile=f541r_dipl_pdb$xyz[2:11,])
rd_f541r <- rmsd(xyz_f541r[1,ca.inds_f541r$xyz], xyz_f541r[,ca.inds_f541r$xyz])

rf_native <- rmsf(xyz_native[,ca.inds_native$xyz])
plot(rf_native, ylab="RMSF", xlab="Residue Position", typ="l",col="pink")
rf_l826f <- rmsf(xyz_l826f[,ca.inds_l826f$xyz])
lines(rf_l826f, ylab="RMSF", xlab="Residue Position", typ="l",col="red")
rf_f117r <- rmsf(xyz_f117r[,ca.inds_f117r$xyz])
lines(rf_f117r, ylab="RMSF", xlab="Residue Position", typ="l",col="blue")
rf_f541r <- rmsf(xyz_f541r[,ca.inds_f541r$xyz])
lines(rf_f541r, ylab="RMSF", xlab="Residue Position", typ="l",col="green")
rf_dipl <- data.frame(cbind(rf_f117r,rf_native,rf_f541r,rf_l826f))
legend(700,5,legend=c("Native","l826f","f117r","f541r"),col=c("black","red","blue","green"), lty=1:2, cex=0.7)

plot(rd_l826f, typ="l", ylab="RMSD", xlab="Frame No.",col="red")
#points(lowess(rd_l826f), typ="l", col="red", lty=2, lwd=2)
lines(rd_native, typ="l", ylab="RMSD", xlab="Frame No.")
#points(lowess(rd_native), typ="l", col="black", lty=2, lwd=2)
lines(rd_f541r, typ="l", ylab="RMSD", xlab="Frame No.",col="blue")
#points(lowess(rd_f541r), typ="l", col="blue", lty=2, lwd=2)
lines(rd_f117r, typ="l", ylab="RMSD", xlab="Frame No.",col="orange")
#points(lowess(rd_f117r), typ="l", col="orange", lty=2, lwd=2)

legend(8,1.5,legend=c("Native", "l826f","f117r","f541r"),col=c("black","red","orange","blue"), lty=1:2, cex=0.7)

setwd('/Users/sayaneshome/Desktop/prelims_stuffs/data/data_analysis_tauro/')
native_tauro <- read.table('hpnn_t_native_phase2_trial2.dat')
f117r_tauro <- read.table('f117r_smd_tauro.dat')
f541r_tauro <- read.table('f541r_tauro_smd.dat')
l826f_tauro <- read.table('l826f_smd_tauro.dat')

#force-position graph 
plot(f117r_tauro$V2,f117r_tauro$V3,type = 'line', ylab="Force", xlab="Time(ns)")
lines(f541r_tauro$V2,f541r_tauro$V3,col = "red")
lines(l826f_tauro$V2,l826f_tauro$V3,col = "green")
lines(native_tauro$V2,native_tauro$V3,col = "purple")
#time-position graph;get the frames from this graph specially native one
plot(native_tauro$V1/100000,smooth(native_tauro$V3),type = 'line',ylab="Force", xlab="Time(ns)",xlim = c(0,10))
lines((f541r_tauro$V1)/100000,smooth(f541r_tauro$V3),col = "orange")
lines((l826f_tauro$V1)/100000,smooth(l826f_tauro$V3),col = "blue")
lines((f117r_tauro$V1)/100000,smooth(f117r_tauro$V3),col = "pink")


f117r_tauro_pdb <- read.pdb('hpnn_f117r_tauro_stride100_full.pdb',multi = TRUE)
native_tauro_pdb <- read.pdb('hpnn_tauro_phase2_trial2_stride100_full.pdb',multi = TRUE)
l826f_tauro_pdb <- read.pdb('hpnn_l826f_tauro_stride100_full.pdb',multi = TRUE)
f541r_tauro_pdb <- read.pdb('hpnn_f541r_tauro_stride100_full.pdb',multi = TRUE)

ca.inds_native <- atom.select(native_tauro_pdb, elety="CA")
xyz_native <- fit.xyz(fixed=native_tauro_pdb$xyz[1,], mobile=native_tauro_pdb$xyz[2:11,])
rd_native <- rmsd(xyz_native[1,ca.inds_native$xyz], xyz_native[,ca.inds_native$xyz])
ca.inds_l826f <- atom.select(l826f_tauro_pdb, elety="CA")
xyz_l826f <- fit.xyz(fixed=l826f_tauro_pdb$xyz[1,], mobile=l826f_tauro_pdb$xyz[2:11,])
rd_l826f <- rmsd(xyz_l826f[1,ca.inds_l826f$xyz], xyz_l826f[,ca.inds_l826f$xyz])
ca.inds_f117r <- atom.select(f117r_tauro_pdb, elety="CA")
xyz_f117r <- fit.xyz(fixed=f117r_tauro_pdb$xyz[1,], mobile=f117r_tauro_pdb$xyz[2:11,])
rd_f117r <- rmsd(xyz_f117r[1,ca.inds_f117r$xyz], xyz_f117r[,ca.inds_f117r$xyz])

rf_native <- rmsf(xyz_native[,ca.inds_native$xyz])
plot(rf_native, ylab="RMSF", xlab="Residue Position", typ="l",col="black")
rf_l826f <- rmsf(xyz_l826f[,ca.inds_l826f$xyz])
lines(rf_l826f, ylab="RMSF", xlab="Residue Position", typ="l",col="red")
rf_f117r <- rmsf(xyz_f117r[,ca.inds_f117r$xyz])
lines(rf_f117r, ylab="RMSF", xlab="Residue Position", typ="l",col="blue")
rf_f541r <- rmsf(xyz_f541r[,ca.inds_f541r$xyz])
lines(rf_f541r, ylab="RMSF", xlab="Residue Position", typ="l",col="green")
rf_tauro <- data.frame(cbind(rf_f117r,rf_native,rf_f541r,rf_l826f))

legend(700,4.3,legend=c("Native","l826f","f117r","f541r"),col=c("black","red","blue","green"), lty=1:2, cex=0.7)


plot(rd_l826f, typ="l", ylab="RMSD", xlab="Frame No.",col="red")
#points(lowess(rd_l826f), typ="l", col="red", lty=2, lwd=2)
lines(rd_native, typ="l", ylab="RMSD", xlab="Frame No.")
#points(lowess(rd_native), typ="l", col="black", lty=2, lwd=2)
lines(rd_f541r, typ="l", ylab="RMSD", xlab="Frame No.",col="blue")
#points(lowess(rd_f541r), typ="l", col="blue", lty=2, lwd=2)
lines(rd_f117r, typ="l", ylab="RMSD", xlab="Frame No.",col="orange")
#points(lowess(rd_f117r), typ="l", col="orange", lty=2, lwd=2)
legend(8,1.5,legend=c("Native","l826f","f117r","f541r"),col=c("black","red","orange","blue"), lty=1:2, cex=0.7)


setwd('/Users/sayaneshome/Desktop/prelims_stuffs/data/Diploptene/aug_18/')
library(bio3d)
native_dipl2 <- read.table('dipl_phase2_native_smd.dat')
native_tauro2 <- read.table('tauro_phase2_native_smd.log.dat')

native_dipl1 <- read.table('dipl_phase1_native_smd.dat')
native_tauro1 <- read.table('tauro_phase1_native_smd.log.dat')




