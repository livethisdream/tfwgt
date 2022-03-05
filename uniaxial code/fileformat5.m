function varargout=fileformat5(varargin)
% usage:  [varargout]=fileformat(varargin)
%         
% inputs: 
%         - filename is an optional argument, which is a string of the name of the 
%           file to be opened.  fileopen asks the user to select a file via a gui 
%           if a filename is not given in the input arguments.
%         
%  outputs: the vectors S-parameters and f
%  functions used:  fileopen (requires user input)
%--------------------------------------------------------------------------

if nargin==0
    fid=fileopen();  % use the fileopen function to open the file - asks the user for a file
else
    fid=fopen(varargin{1});  % use the filename input by the user
end

fseek(fid,0,'eof');
fileend=ftell(fid);
fseek(fid,0,'bof');

line=fgetl(fid);
while strncmp(line,'DATA',4)==0
    line=fgetl(fid);
end

% count the number of variables we'll have
numvars=0;
varnames={};
while strncmp(line,'DATA',4)==1
    numvars=numvars+1;
    varnames(numvars)={regexprep(regexprep(line,'DATA ',''),' RI','')};
    line=fgetl(fid);
end

% frequencies will always come first
while strcmp(line,'VAR_LIST_BEGIN')==0  % throw out the header stuff
    line=fgetl(fid);
end
line=fgetl(fid);
f=[];
while strcmp(line,'VAR_LIST_END')==0  % do this until we get to the end of f data
    f=vertcat(f,str2double(line));
    line=fgetl(fid);
end

% now, get the variables
A=zeros(length(f),numvars*2+1); % pre-allocate (freq, then real, imag for each s-param)
A(:,1)=f; % assign the frequency vals
for n=2:2:size(A,2) % we already assigned freqency vals
    line=fgetl(fid);
    while strcmp(line,'BEGIN')==0 % wait until we find the beginning of the data
        line=fgetl(fid);
    end
    line=fgetl(fid);
    idx=1;
    while strcmp(line,'END')==0 % loop until we find the end of the data
        A(idx,n:n+1)=sscanf(strrep(line,',',' '),'%f')';
        line=fgetl(fid);
        idx=idx+1;
    end
    % we now have the data in an array where A(:,1) is the frequency and 
    % A(:,even) is the real part and A(:,odd) is imaginary part
end
% for m=1:size(A,2)
%     varargout(m)=A(:,m);
% end
fclose(fid);

for ii=1:length(varnames)
   name=cell2mat(varnames(ii));
   rmatch=strfind(name,'R');
   if isempty(rmatch)==0
      varnames(ii)={[regexprep(name,'/','') name(end)]}; 
   end
end

if nargout==1
    varargout={A};
elseif nargout==2
    varargout(1)={A};
    varargout(2)={varnames};
end

end % end of function

