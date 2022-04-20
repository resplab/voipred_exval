#GRconnect("voipredex")
#GRupload("./R/server.R")

for(id in 1:10)
{
  system(paste("Rscript.exe -e \" library('GRcomp'); GRconnect('voipred_exval','meta_sim'); GRrun('meta_sim.R',",id,") \" "), invisible = F, wait=F)
}



