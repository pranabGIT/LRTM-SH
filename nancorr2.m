function A=nancorr2(U,V)
    id1=find(isnan(U));
     id2=find(isnan(U));
     id=union(id1,id2);
     U(id)=[];
     V(id)=[];
     A=corr2(U,V);
end