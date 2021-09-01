library(bio3d)
pdb1 <- read.pdb('hpnn_tauro_phase2_transport_stride30.pdb',multi = TRUE)
ca.inds_pdb1 <- atom.select(pdb1,elety = 'CA')
xyz_pdb1 <- pdbfit(pdb1,inds = ca.inds_pdb1)
rf1 <- data.frame(rmsf(xyz_pdb1[,ca.inds_pdb1$xyz]))
rf1$residue_no <- rownames(rf1)
colnames(rf1)[1] <- "RMSF_values"
plot(rf1$residue_no,rf1$RMSF_values, ylab="RMSF", xlab="RMSF versus residues", typ="l",col = 'red',lwd = 2)

pdb2 <- read.pdb('hpnn_tauro_f541r_stride30.pdb',multi = TRUE)
ca.inds_pdb2 <- atom.select(pdb2,elety = 'CA')
xyz_pdb2 <- pdbfit(pdb2,inds = ca.inds_pdb2)
rf2 <- data.frame(rmsf(xyz_pdb2[,ca.inds_pdb2$xyz]))
rf2$residue_no <- rownames(rf2)
colnames(rf2)[1] <- "RMSF_values"
lines(rf2$residue_no,rf2$RMSF_values, typ="l",col = 'blue',lwd = 2)

pdb3 <- read.pdb('hpnn_tauro_l826f_stride30.pdb',multi = TRUE)
ca.inds_pdb3 <- atom.select(pdb3,elety = 'CA')
xyz_pdb3 <- pdbfit(pdb3,inds = ca.inds_pdb3)
rf3 <- data.frame(rmsf(xyz_pdb3[,ca.inds_pdb3$xyz]))
rf3$residue_no <- rownames(rf3)
colnames(rf3)[1] <- "RMSF_values"
lines(rf3$residue_no,rf3$RMSF_values, typ="l",col = 'orange',lwd = 2)

pdb4 <- read.pdb('hpnn_tauro_f117r_stride30.pdb',multi = TRUE)
ca.inds_pdb4 <- atom.select(pdb4,elety = 'CA')
xyz_pdb4 <- pdbfit(pdb4,inds = ca.inds_pdb4)
rf4 <- data.frame(rmsf(xyz_pdb4[,ca.inds_pdb4$xyz]))
rf4$residue_no <- rownames(rf4)
colnames(rf4)[1] <- "RMSF_values"
lines(rf4$residue_no,rf4$RMSF_values, typ="l",col = 'purple',lwd = 2)
text(locator(), labels = c("native", "f541r","l826f","f117r")))

pdb1a <- read.pdb('dipl_2ndphase_native_stride30.pdb',multi = TRUE)
ca.inds_pdb1a <- atom.select(pdb1a,elety = 'CA')
xyz_pdb1a <- pdbfit(pdb1a,inds = ca.inds_pdb1a)
rf1a <- data.frame(rmsf(xyz_pdb1a[,ca.inds_pdb1a$xyz]))
rf1a$residue_no <- rownames(rf1a)
colnames(rf1a)[1] <- "RMSF_values"
plot(rf1a$residue_no,rf1a$RMSF_values, ylab="RMSF", xlab="RMSF versus residues", typ="l",col = 'orange')

pdb2a <- read.pdb('hpnn_dipl_l826_stride30_edit.pdb',multi = TRUE)
ca.inds_pdb2a <- atom.select(pdb2a,elety = 'CA')
xyz_pdb2a <- pdbfit(pdb2a,inds = ca.inds_pdb2a)
rf2a <- data.frame(rmsf(xyz_pdb2a[,ca.inds_pdb2a$xyz]))
rf2a$residue_no <- rownames(rf2a)
colnames(rf2a)[1] <- "RMSF_values"
lines(rf2a$residue_no,rf2a$RMSF_values, typ="l",col = 'blue')

pdb3a <- read.pdb('hpnn_dipl_f117r_stride30.pdb',multi = TRUE)
ca.inds_pdb3a <- atom.select(pdb3a,elety = 'CA')
xyz_pdb3a <- pdbfit(pdb3a,inds = ca.inds_pdb3a)
rf3a <- data.frame(rmsf(xyz_pdb3a[,ca.inds_pdb3a$xyz]))
rf3a$residue_no <- rownames(rf3a)
colnames(rf3a)[1] <- "RMSF_values"
lines(rf3a$residue_no,rf3a$RMSF_values, typ="l",col = 'green')

pdb4a <- read.pdb('hpnn_dipl_f541r_stride30_edit.pdb',multi = TRUE) #needs change
ca.inds_pdb4a <- atom.select(pdb4a,elety = 'CA')
xyz_pdb4a <- pdbfit(pdb4a,inds = ca.inds_pdb4a)
rf4a <- data.frame(rmsf(xyz_pdb4a[,ca.inds_pdb4a$xyz]))
rf4a$residue_no <- rownames(rf4a)
colnames(rf4a)[1] <- "RMSF_values"
lines(rf4a$residue_no,rf4a$RMSF_values, typ="l",col = 'purple')

cij1<-dccm(xyz_pdb1[,ca.inds_pdb1$xyz])
plot(cij1)
cij2<-dccm(xyz_pdb2[,ca.inds_pdb2$xyz])
plot(cij2)
cij3<-dccm(xyz_pdb3[,ca.inds_pdb3$xyz])
plot(cij3)
cij4<-dccm(xyz_pdb4[,ca.inds_pdb4$xyz])
plot(cij4)

cij1a<-dccm(xyz_pdb1a[,ca.inds_pdb1a$xyz])
plot(cij1a)
cij2a<-dccm(xyz_pdb2a[,ca.inds_pdb2a$xyz])
plot(cij2a)
cij3a<-dccm(xyz_pdb3a[,ca.inds_pdb3a$xyz])
plot(cij3a)
cij4a<-dccm(xyz_pdb4a[,ca.inds_pdb4a$xyz])
plot(cij4a)


pdb_d <- read.pdb('hpnn_dipl_fulltransport.pdb',multi = TRUE) #needs change
ca.inds_pdb_d <- atom.select(pdb_d,elety = 'CA')
xyz_pdb_d <- pdbfit(pdb_d,inds = ca.inds_pdb_d)
rf_d <- data.frame(rmsf(xyz_pdb_d[,ca.inds_pdb_d$xyz]))
rf_d$residue_no <- rownames(rf_d)
colnames(rf_d)[1] <- "RMSF_values"
plot(rf_d$residue_no,rf_d$RMSF_values, typ="l",col = 'purple')

pdb_t <- read.pdb('hpnn_tauro_fulltransport.pdb',multi = TRUE)
ca.inds_pdb_t <- atom.select(pdb_t,elety = 'CA')
xyz_pdb_t <- pdbfit(pdb_t,inds = ca.inds_pdb_t)
rf_t <- data.frame(rmsf(xyz_pdb_t[,ca.inds_pdb_t$xyz]))
rf_t$residue_no <- rownames(rf_t)
colnames(rf_t)[1] <- "RMSF_values"
lines(rf_t$residue_no,rf_t$RMSF_values,typ="l",col = 'blue')

#community clustering

cij1<-dccm(xyz_pdb1[,ca.inds_pdb1$xyz],cutoff.cij=0.3)
test1<-cna(cij1,cm=cm)
cm1<- cmap(xyz_pdb1[,ca.inds_pdb1$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test1,pdb1)
