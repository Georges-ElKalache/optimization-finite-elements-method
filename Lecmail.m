function [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=Lecmail(nomfile)
fid = fopen(sprintf('%s.amdba',nomfile),'r') ;
N = fscanf(fid,'%i') ;
Nbpt = N(1) ;
Nbtri = N(2) ;
line = fgets(fid) ;
tmp = fscanf(fid,'%f',[4,Nbpt]) ;
Coorneu = tmp(2:3,:)' ;
Refneu = tmp(4,:)' ;
tmp = fscanf(fid,'%i',[5,Nbtri]) ;
Numtri=tmp(2:4,:)' ;
Reftri=tmp(5,:)' ;