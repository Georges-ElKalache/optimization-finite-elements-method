function fig = tri_contour(elm_som,som_x,som_y,som_u,hautr)
% TRI_CONTOUR(ELM_SOM,SOM_X,SOM_Y,SOM_U,HAUTR) trace les lignes de
% contours (iso-lignes ou iso-valeurs definies par le tableau hautr) 
% de la fonction som_u (vecteur des valeurs aux sommets du maillage)
% 
% elm_som = topologie du maillage
% som_x = vecteur ligne des 1e coord des sommets
% som_y = vecteur ligne des 2e coord des sommets
%
% Copyright (c) 2005 by Frederic Pascal, ENS de CACHAN
  
if nargin < 5
    error('Nombres arguments insuffisants')
end

% Taille du pb
nbe = size(elm_som,1);

% Transformation eventuelle du vecteur ligne en vecteur colonne
if size(som_u,1)==1 
    v = som_u';
else
    v = som_u;
end

% hautr doit etre un vecteur ligne
if size(hautr,2)== 1
    hautr = hautr';
end
% Rien a tracer
if size(hautr,2)== 0
    error('Rien a tracer');
end

% Initialisation
plt_n =zeros(1,size(hautr,2));

map=jet(128);

i=0;
nb=length(hautr);

for h=hautr
     
    i = i + 1;
    
    clear plt_x
    clear plt_y
    
    diff = v(elm_som(:,:))-h;
    
    sdiff=sign(diff);
    
    ind=find(abs(sum(sdiff,2))==1 & prod(sdiff,2)~=0);

    edges1=sdiff(ind,[2,3]);
    edges2=sdiff(ind,[3,1]);
    edges3=sdiff(ind,[1,2]);
    
    ind1=find(sum(edges1,2)==0);
    
    aux(ind1,1)=diff(ind(ind1),2)./(diff(ind(ind1),2)-diff(ind(ind1),3));
    
    plt_x(1,ind1) = som_x(elm_som(ind(ind1),2)) - aux(ind1).*(som_x(elm_som(ind(ind1),2))-som_x(elm_som(ind(ind1),3)));
    plt_y(1,ind1) = som_y(elm_som(ind(ind1),2)) - aux(ind1).*(som_y(elm_som(ind(ind1),2))-som_y(elm_som(ind(ind1),3)));
     
    ind2=find(sum(edges1,2)==0 & sum(edges2,2)==0);
    
    aux(ind2,1)=diff(ind(ind2),3)./(diff(ind(ind2),3)-diff(ind(ind2),1));
    
    plt_x(2,ind2) = som_x(elm_som(ind(ind2),3)) - aux(ind2).*(som_x(elm_som(ind(ind2),3))-som_x(elm_som(ind(ind2),1)));
    plt_y(2,ind2) = som_y(elm_som(ind(ind2),3)) - aux(ind2).*(som_y(elm_som(ind(ind2),3))-som_y(elm_som(ind(ind2),1)));
    
    ind2=find(sum(edges3,2)==0 & sum(edges2,2)==0);
    
    aux(ind2,1)=diff(ind(ind2),3)./(diff(ind(ind2),3)-diff(ind(ind2),1));
    
    plt_x(1,ind2) = som_x(elm_som(ind(ind2),3)) - aux(ind2).*(som_x(elm_som(ind(ind2),3))-som_x(elm_som(ind(ind2),1)));
    plt_y(1,ind2) = som_y(elm_som(ind(ind2),3)) - aux(ind2).*(som_y(elm_som(ind(ind2),3))-som_y(elm_som(ind(ind2),1)));
    
    ind3=find(sum(edges3,2)==0);
    
    aux(ind3,1)=diff(ind(ind3),1)./(diff(ind(ind3),1)-diff(ind(ind3),2));
    
    plt_x(2,ind3) = som_x(elm_som(ind(ind3),1)) - aux(ind3).*(som_x(elm_som(ind(ind3),1))-som_x(elm_som(ind(ind3),2)));
    plt_y(2,ind3) = som_y(elm_som(ind(ind3),1)) - aux(ind3).*(som_y(elm_som(ind(ind3),1))-som_y(elm_som(ind(ind3),2)));
    
    a=plot(plt_x,plt_y);
    set(a,'Color',map(round(i*128/nb),:))
    hold on
    
end
hold off
