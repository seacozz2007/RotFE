function tF=rot2templ(Rot,File,qFile)
%rot2templ(Rot,File)
%
% convert a ROTOR struct (output of rotfe) into a
% template file callable by rotfe modeplot etc.
%
%Used as a utility by various functions
%
% by I. bucher 8.1998
% revised 6/99 ignores struct sub-fileds of Rot (caused error)
% revised 8/99 now includes struct sub-fields so that reduced order models can be handeled
% 
tF=[];
if nargin<2, disp('rot2templ(Rot,File),  => 2 arguments are needed')
   error('existing')
end
if nargin<=2, qFile=0; end
if nargin>2, qFile=1; end

if isa(Rot,'struct');
   Rot=rm_mats(Rot,{'M' 'K' 'G' 'D' 'KH' 'Kst'});
   f=fieldnames(Rot);
   
   n=size(f,1);  

   fid=fopen(File,'w');
   if fid==-1, error(' could not open file in rot2templ'); end
        
 if (Rot.dim>25) & (qFile==0) % too large for an ascii file
	 tFILE=['tempmat.mat']
        fprintf(' variables very large, stored in temp. file %s \n',tFILE);
        tF=tFILE;
     for q=1:n,      
        eval([f{q} '=Rot.' f{q} ';' ]);
        if q==1,	
		save(tFILE,f{q});
	else,
		save(tFILE,f{q},'-APPEND');
	end
      end
     fprintf(fid,'%c temporary file \n %c description of a Rot structure \n','%','%');    
     fprintf(fid,'%c Due to large size data stored in %s \n ','%',tFILE)    ;
      fprintf(fid,'\n load %s \n',tFILE);
      fprintf(fid,'\n %c%s \n','%',' end');

else
      for q=1:n,      
        x=eval(['Rot.' f{q} ]);
        if isa(x,'struct') 
          write_struct(fid,x,f{q});
         elseif isa(x,'cell'),
         fprintf(' ignoring struct sub-field %s\n',f{q})
      else
         if ~isstr(x)            
            fprintf(fid,'%c ==> %s \n',37,f{q});
            cmd=[f{q} ' = ... ' ];
            fprintf(fid,'%c-------------------------------- \n',37);
            fprintf(fid,' %s  \n %s ; \n',cmd,mat2str(x));
            fprintf(fid,'%c-------------------------------- \n',37);            
         end
	      end
   end
 end
   fclose(fid);
else
   error(' Rot is not a valid struct')
end

function r=rm_mats(r,NAMEs)
   f=fieldnames(r);
   n=size(f,1);
   for q=1:n, 
      if ismember(f{q},NAMEs)
         r=rmfield(r,f{q});
      end
   end
   
    
%<<<<<<<<<<<<<<<<<<<<<   
function write_struct(fid,x,NAME)
% printf a structure x recursively
   f=fieldnames(x);
   n=size(f,1);
   for q=1:n, 
 if isa(x,'struct') 
   f=fieldnames(x);
   n=size(f,1);
   for q=1:n, 
%      fprintf(' evaluating %s q=%d n=%d \n',['x.' f{q}],q,n)
      z=eval(['x.' f{q}]);
      write_struct(fid,z,[NAME '.' f{q} ]);
   end
elseif isa(x,'cell')
    fprintf(' ignoring %s . %s \n',NAME,x)
else
   fprintf(fid,'%c ==> %s \n',37,[NAME ]);
   if isstr(x), cmd=[NAME  ' =' char(39)  x  char(39)];
   else,cmd=[NAME  ' =' mat2str( x )];
      end
   fprintf(fid,'%c-------------------------------- \n',37);
   fprintf(fid,' %s ; \n',cmd);
   fprintf(fid,'%c-------------------------------- \n',37);            
end
  
end
