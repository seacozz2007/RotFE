function [NODES ELM]=get_elem_shaft(W)

n=size(W,1);
NODES=[];
ELM=[];
x=0;
for q=1:n-1
    
    
    w=W(q,1);
    
    y=W(q,2)/2;
    
    NODES=[NODES x];
    ELM=[ELM ;q q+1 W(q,2) 0  1];
    
    plot([x x+w x+w x x],[-y -y y y -y],'b')
    
    hold on
   
    x=x+w;
end
    NODES=[NODES x];


shg
hold off
axis equal
 