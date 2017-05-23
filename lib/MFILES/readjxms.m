function [hmap,lo,la,age,Plm,hmapref]=readjxms(fnames,minage,maxage,ageref,degres,L3,Plm)
% [hmap,lo,la,age,Plm,hmapref]=READJXMS(fnames,minage,maxage,ageref,degres,L3)
% 
% INPUT:
%
% fnames        String with path of files (wildcards allowed)
% minage		minimum age of files to load
% maxage		maximum age of files to load
% ageref		age of reference datum (pass NaN for none)
% degres        The degree reslution of the output grid OR
%               a matrix with [lon(:) lat(:)] in degrees
% L3            An optional extra truncation for the expansions
%
% OUTPUT:
%
% hmap          The data maps or vectors appropriately referenced
% lo            Longitudes of these data, in degrees
% la            Latitudes of these data, in degrees
% age			Ages of the data
% Plm			Legendre polynomial expansion
% hmapref		Reference map or vectors
% 
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Jul 9 21:14:08 EDT 2013

% Define defaults
defval('fnames','seagl_*.mn')
defval('minage',0)
defval('maxage',200)
defval('L3',256)
defval('ageref',0)
defval('degres',0.7)

fnames2=[fnames '.gz'];
gzfiles=dir(fnames2);
if length(gzfiles)>0
	disp('Unzipping files...');
	for i=1:length(gzfiles)
		system(['gzip -d ' gzfiles(i).name]);
	end
end

u=strfind(fnames,'/');
if length(u)>0
	fpath = fnames(1:u(end));
else
	fpath = '';
end

ffiles = dir(fnames);
for i=1:length(ffiles)
	spec = readheader(fullfile(fpath,ffiles(i).name));
	L1s(i) = spec(1); age(i) = spec(2);
end

fname2 = 0;
if isfinite(ageref)
	sub = find(age==ageref);
	if length(sub)>0
		fname2 = fullfile(fpath,ffiles(sub(1)).name);
		L2 = L1s(i);
		disp(['      AGE ' num2str(ageref)]); 
	end
end

sub=find((age>=minage).*(age<=maxage));
age=age(sub); ffiles=ffiles(sub); L1s=L1s(sub);
[agesort,agesorti] = sort(age);
age=agesort(end:-1:1); ffiles=ffiles(agesorti(end:-1:1)); L1s=L1s(agesorti(end:-1:1));

% Load the reference model
if ~isstr(fname2)
  hmapref=fname2;
else
  [hlm2,dels2,dems2]=readit(fname2);
  if L3<L2
    % Restrict?
    hlm2=hlm2(1:2*addmup(L3),:); 
    dels2=dels2(1:addmup(L3));
    dems2=dems2(1:addmup(L3));
    L2=L3;
  end

  if prod(size(degres))==1
    % It's a map!
    [hmapref,lo,la,Plm]=plm2xyz(...
		[dels2 dems2 reshape(hlm2,2,addmup(L2))'], degres);
	hmapref=fliplr(hmapref);
  elseif size(degres,2)==2
    % It's an irregular grid!
    if exist('Plm','var')
     [hmapref,lo,la,Plm]=plm2xyz([dels2 dems2 reshape(hlm2,2,addmup(L2))'],... 
			degres(:,2),360-degres(:,1),[],[],Plm);
   else
  	  [hmapref,lo,la,Plm]=plm2xyz([dels2 dems2 reshape(hlm2,2,addmup(L2))'], ...
			degres(:,2),360-degres(:,1));
	end
  end
end

if length(ffiles)>0
	if exist('Plm','var')
		startind = 1;
	else
		startind =2 ;
		if prod(size(degres))==1
		  % It's a map!
		  % Initialize the data maps
		  nlon=ceil(360/degres+1);
		  nlat=ceil(180/degres+1);
		  hmap=nan(nlat,nlon,length(ffiles));
		elseif size(degres,2)==2
		  hmap=nan(size(degres,1),length(ffiles));
		end

		ind=1;
		disp(['      AGE ' num2str(age(ind))]); 
		fname1=fullfile(fpath,ffiles(ind).name);
		L1=L1s(ind);
		
		% Load the data model as coefficients
		[hlm1,dels1,dems1]=readit(fname1);
		
		% Restrict?
		if L3<L1
		  hlm1=hlm1(1:2*addmup(L3),:);
		  dels1=dels1(1:addmup(L3));
		  dems1=dems1(1:addmup(L3));
		  L1=L3;
		end
			
		% Get the data maps
		  if prod(size(degres))==1
			% It's a map!
		    [hmap(:,:,ind),lo2,la2,Plm2]=plm2xyz(...
		    	[dels1 dems1 reshape(hlm1,2,addmup(L1))'],degres);
			hmap(:,:,ind) = fliplr(hmap(:,:,ind))-hmapref;
		  elseif size(degres,2)==2
			% It's an irregular grid!
			
			[hmap(:,ind),lo2,la2,Plm]=plm2xyz([dels1 dems1 reshape(hlm1,2,addmup(L1))'],... 
						degres(:,2),360-degres(:,1));
			hmap(:,ind) = hmap(:,ind)-hmapref;
		  end
	end


	if prod(size(degres))==1
		%It's a map
		
		parfor ind=startind:length(ffiles)
			disp(['      AGE ' num2str(age(ind))]); 
			fname1=fullfile(fpath,ffiles(ind).name);
			L1=L1s(ind);
			
			% Load the data model as coefficients
			[hlm1,dels1,dems1]=readit(fname1);
			
			% Restrict?
			if L3<L1
			  hlm1=hlm1(1:2*addmup(L3),:);
			  dels1=dels1(1:addmup(L3));
			  dems1=dems1(1:addmup(L3));
			  L1=L3;
			end
				
			% Get the data maps
		    [hmap(:,:,ind),lo2,la2,Plm2]=plm2xyz(...
		    	[dels1 dems1 reshape(hlm1,2,addmup(L1))'],degres);
			hmap(:,:,ind) = fliplr(hmap(:,:,ind))-hmapref;
		end
	elseif size(degres,2)==2
		% It's an irregular grid!
		parfor ind=2:length(ffiles)
			disp(['      AGE ' num2str(age(ind))]); 
			fname1=fullfile(fpath,ffiles(ind).name);
			L1=L1s(ind);
			
			% Load the data model as coefficients
			[hlm1,dels1,dems1]=readit(fname1);
			
			% Restrict?
			if L3<L1
			  hlm1=hlm1(1:2*addmup(L3),:);
			  dels1=dels1(1:addmup(L3));
			  dems1=dems1(1:addmup(L3));
			  L1=L3;
			end
				
			% Get the data maps
			[hmap(:,ind),lo2,la2]=plm2xyz([dels1 dems1 reshape(hlm1,2,addmup(L1))'],... 
					degres(:,2),360-degres(:,1),[],[],Plm);
			hmap(:,ind) = hmap(:,ind)-hmapref;
		end
	end
else
	hmap = hmapref;
end

if ~isstr(fname2)
	lo = lo2; la = la2;
	if ~exist('Plm','var')
		if exist('Plm2','var')
			Plm = Plm2;
		end
	end
end


%%%%%%%
function [specs]=readheader(fname)

fid = fopen(fname);
tline = fgetl(fid);
specs = sscanf(tline,'%f');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hlm,dels,dems]=readit(fname)
% Load the file
fid=fopen(fname,'r');
% Read the header
head=fgetl(fid);
LT=sscanf(head,'%f',2);
% Spherical harmonic bandwidth
L=LT(1);
% Time
T=LT(2);
% Read the coefficients
COF=fscanf(fid,'%f');
fclose(fid);
% Check that the length is as I expect it
difer(length(COF)-2*addmup(L)-2*(L/2+1),[],[],NaN)

% Now rearrange
mods=(L-[0:L]+1);
% The indices of the blank pairs
modu=cumsum(mods+mod(mods,2));
modb=modu(~~mod(mods,2));
blanx=sort([2*modb-1 2*modb]);
% Check these are indeed blanks
difer(COF(blanx),[],[],NaN)
% Remaining are the non-blanks
COF=skip(COF,blanx);
% Which should be the same number as expected
difer(length(COF)-2*addmup(L),[],[],NaN)

% Now split into the real and imaginary parts
C=COF(1:2:end);
S=COF(2:2:end);
% The first L+1 should be zero
difer(S(1:L+1),[],[],NaN)

% Create a blank array with FJS format
[dems,dels,mz,lmcosi]=addmon(L);
% Fake arrival to align with readshds
hlm=lmcosi(:,3:4);
modm=[cumsum([0 mods(1:end-1)])+1 addmup(L)+1];
% Populate this with the JXM coefficients
for m=0:L
  % Stick in the "cosine" coefficients
  hlm(addmup(m-1:L-1)+m+1,1)=C(modm(m+1):modm(m+2)-1);
  % Stick in the "sine" coefficients
  hlm(addmup(m-1:L-1)+m+1,2)=S(modm(m+1):modm(m+2)-1);
end

% Some conventions need to be realigned
CSC=(-1).^dems;
dom=sqrt(2-(dems==0));
hlm(:,1)=hlm(:,1).*CSC.*dom;
hlm(:,2)=hlm(:,2).*CSC.*dom;
% Fake arrival to align with readshds
hlm=hlm';
hlm=hlm(:);

