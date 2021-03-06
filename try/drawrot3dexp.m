function drawrot3dexp(Rot,full)
%drawrot3dexp(Rot)
% or
%drawrot3dexp(Rot,'full')
%
%Description
% draws (3D)  the model of a rotor 
%  contained in a Rot structure
%   drawrot3dexp(Rot) - draws  the mode
%and
%   drawrot3dexp(Rot,'anything') - draws  the model with a1/4 section
%
% Copyright by
% I. Bucher September 1998
% under the GNU liecense 

if nargin<3, clf, end

cQ=0; nQ=0;
nF=(nargin>1); % full or section flag
if nargin<2, CallBack='    ';end 

[ne m]=size(Rot.ELEMENTS);
Lt=norm(Rot.NODES)/ne;			% typical length
Dt=norm(Rot.ELEMENTS(:,3))/(ne); %typical diameter

% >>>>> Draw elements
for q=1:ne,
   do=Rot.ELEMENTS(q,3);   di=Rot.ELEMENTS(q,4); 
   x1=Rot.NODES(Rot.ELEMENTS(q,1)); 
   x2=Rot.NODES(Rot.ELEMENTS(q,2)); 
   
   x=[x1 x2 x2 x1];
   y1=[di/2 di/2 do/2 do/2];
   set(gcf,'rende','zbuffer')
   c=Rot.ELEMENTS(q,5);
   set(gca,'drawmode','fast');   
   h2=draw_cyl([x1;x2],[do/2;do/2],24,c,di/2,nF);
   set(gca,'drawmode','fast');   
   h1=draw_cyl([x1;x2],[di/2;di/2],24,c,di/2,nF);
end


% >>>>> Draw Discs
if isfield(Rot,'DISCS')
[nd tmp]=size(Rot.DISCS);
Do_max=0;
for q=1:nd,
   x1=Rot.NODES(Rot.DISCS(q,1));	% center coords
   Do=Rot.DISCS(q,2);   Di=Rot.DISCS(q,3);
   W=Rot.DISCS(q,4);
   c=Rot.DISCS(q,5);
   x=[x1-W/2 x1+W/2 x1+W/2 x1-W/2];
   y=[Di/2 Di/2 Do/2 Do/2];
   h=draw_cyl([x1-W/2;x1+W/2 ],[Do/2 Do/2],32,c,Di/2,nF);
   Do_max=max(Do_max,Do);
end
end
% >>>>> Draw Springs
if isfield(Rot,'SPRINGS')
[ns tmp]=size(Rot.SPRINGS);
th=linspace(0,2*pi,1000); 
z=linspace(0,1,length(th));
ncoils=7; D3=Dt;
for q=1:ns,
   x1=Rot.NODES(Rot.SPRINGS(q,1));	% center coords
   v=Rot.SPRINGS(q,2:5); v=[v(1) v(3);v(4) v(2)];
   % find principle axes of stiffness    
   [T dd]=eig(v); dd=sqrt(diag(dd)); dd=dd/max(dd);
   h1=line([x1 x1],[-Do_max Do_max]/2 ,[0 0]);
   h2=line([x1 x1],[0 0],[-Do_max Do_max] /2);
   
   draw_spring(x1+D3*cos(ncoils*th),z*Lt,D3*sin(ncoils*th),T);
   draw_spring(x1+D3*cos(ncoils*th),-z*Lt,D3*sin(ncoils*th),T);
   draw_spring(x1+D3*sin(ncoils*th),D3*cos(ncoils*th),z*Lt,T);
   draw_spring(x1+D3*sin(ncoils*th),D3*cos(ncoils*th),-z*Lt,T);
   
   set([h1 h2],'linewidth',3,'linestyle','-.');
end
end
% >>>>> Draw BC
if isfield(Rot,'BCNodeDir')
[nb tmp]=size(Rot.BCNodeDir(:));
d_shaft_max=max(max(Rot.ELEMENTS(:,3:4)));
for q=1:nb,
   x1=Rot.NODES(fix(Rot.BCNodeDir(q)));	% center coords
   c=[1 1 1];
   z=[x1 x1-Dt*2 x1+Dt*2  x1];
   y=[0  Dt*4 Dt*4 0];
   x=[0 0 0 0];
   h1=plot3(z,x,y);
   h2=plot3(z,x,-y);
   h3=plot3(z,y,x);
   h4=plot3(z,-y,x);   
end
end
% >>>>> Draw Point mass
if isfield(Rot,'POINT_MASS')
[np tmp]=size(Rot.POINT_MASS);
%d1=max(d_shaft_max/2,Dt/2);
for q=1:np,
   x1=Rot.NODES(Rot.POINT_MASS(q,1));	% center coords
   jj=union(find(Rot.ELEMENTS(:,1)==q),...
      find(Rot.ELEMENTS(:,2)==q));
   d1=max(Rot.ELEMENTS(jj,3))/1.7;
   c=[1 0 0];
   theta=[0:11]'*2*pi/12;
   x=[Dt*2*cos(theta) ];
   y=[+Dt*2*sin(theta)];
   %   h1=patch(x,y,3*x1*ones(size(x)));
   [x y z]=sphere(12);
   h1=surf(x1+z*d1,x*d1,y*d1);
   set(h1,'facecolor','r'); 
end
end

axis equal
axis off
figure(gcf)
hold off


set(findobj('type','surface'),'FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5,'EdgeColor','interp')
light('Position',[0.9 -10 -30],'Style','infinite');
light('Position',[10 10 30],'Style','infinite');
light('Position',[0.9 10 -30],'Style','infinite');



 
function h=draw_cyl(z,r,n,c,rmin,nF)
%
m = length(r); 
n1=fix(n*.75);
if ~nF, n1=n;end
theta = (0:n1)/n*2*pi+1e-3+pi;
sintheta = sin(theta); %sintheta(n1) = 0;

x = r(:) * cos(theta);
y = r(:) * sintheta;
z = z(:) * ones(1,n1+1);

% draw outer cylinder
h1=surf(z,x,y,c*ones(size(z)));
hold on
%   h=draw_cyl([x1-W/2;x1+W/2 ],[Do/2 Do/2],32,c,Di/2);

if nargin<5
   rmin=0;
end
rmax=max(r);
for q=1:n1
   x1=[rmin rmax rmax rmin].'*cos(theta(q:q+1)); 
   y1=[rmin rmax rmax rmin].'*sintheta(q:q+1);
   h2=surf(z(1)*ones(size(x1)),x1,y1,c*ones(size(x1)));
   h3=surf(z(2)*ones(size(x1)),x1,y1,c*ones(size(x1)));
   h1=[h1(:) ; h2(:); h3(:)];
end

h=h1;        

r=[rmin rmax rmax rmin]; t1=theta(1); t2=theta(n1+1);
h4=fill3([z(1) z(1) z(2) z(2)]',r*cos(t1),r*sin(t1),c*ones(size(r)));

h5=fill3([z(1) z(1) z(2) z(2)]',r*cos(t2),r*sin(t2),c*ones(size(r)));

set(h4,'tag','in1')
set(h5,'tag','in1')
%set(h4,'facecolor','r','edgecolor','w')             
%set(h5,'facecolor','r','edgecolor','w')
set(h4,'edgecolor','w')             
set(h5,'edgecolor','w')


function draw_spring(x,y,z,T)
%
 yz=[y(:) z(:)]; yz=yz*T.';
 h=plot3(x,yz(:,1),yz(:,2));
 set(h,'linewidth',2)
