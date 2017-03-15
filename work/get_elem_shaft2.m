function [NODES ELM ww] =get_elem_shaft2(W,D)
% same  as get_elem_shaft but splits elements
% so that they are not longer than D

n=size(W,1);
ww=[];
for q=1:n
     n1 = ceil(W(q,1)/D);
     dx=W(q,1)/n1;
     ww=[ww; [dx*ones(n1,1) W(q,2)*ones(n1,1) W(q,3)*ones(n1,1)]];
end    

n=size(ww,1);

NODES=[];
ELM=[];
x=0;

 
for q=1:n-1

    w=ww(q,1);
    y=ww(q,2)/2;
    
    NODES=[NODES x];
        
    ELM=[ELM ;q q+1 ww(q,2) 0  ww(q,3)];
     
    
    x=x+w;
end
    NODES=[NODES x x+ww(end,1)];
    ELM=[ELM ;n n+1 ww(end,2) 0  1];

shg
hold off
axis equal
 