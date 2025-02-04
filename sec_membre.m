function f = sec_membre(x,y,nomfile,A0)
[Nbpt,Nbtri,Coorneu,~,Numtri,~]=Lecmail(nomfile);
%construction de f du second membre
f=sparse(Nbpt,1);
for l=1:Nbtri
    %Coordonnees des trois sommets du triangle l
    x1=Coorneu(Numtri(l,1),1);y1=Coorneu(Numtri(l,1),2);
    x2=Coorneu(Numtri(l,2),1);y2=Coorneu(Numtri(l,2),2);
    x3=Coorneu(Numtri(l,3),1);y3=Coorneu(Numtri(l,3),2);
    delta=(x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
    %calcul des trois coordonnees barycentriques
    lam1=(1/delta)*((y2-y3)*(x-x3)-(x2-x3)*(y-y3));
    lam2=(1/delta)*((y3-y1)*(x-x1)-(x3-x1)*(y-y1));
    lam3=(1/delta)*((y1-y2)*(x-x2)-(x1-x2)*(y-y2));
    %voir ou se trouve la resistance pour utiliser la bonne formule
    if (0<=lam1) && (0<=lam2) && (0<=lam3) && (1>=lam1) && (1>=lam2) && (1>=lam3)
        for k=1:3
            if k==1
                f(Numtri(l,k))=A0*lam1;
            end
            if k==2
                f(Numtri(l,k))=A0*lam2;
            end
            if k==3
                f(Numtri(l,k))=A0*lam3;
            end
        end
    end
end


            
        
    

