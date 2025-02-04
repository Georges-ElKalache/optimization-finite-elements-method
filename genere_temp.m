function Temp=genere_temp(nomfile,n)
[Nbpt,~,~,~,~,~]=Lecmail(nomfile);
R=genere_R_alea(n);
R(1,2) = 0.75;
R(2,2) = 0.75;
R(3,3) =-0.75;
R(2,3) =-0.75;
R(1,3) = 0.75;
R(3,2) = -0.75;
Temp=sparse(Nbpt,n);
for i=1:n
    Temp(:,i)=prob_direct(R(i,2),R(i,3),nomfile,0,1e4);
end
end