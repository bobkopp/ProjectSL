function [hmap,lo,la,age,hmapref]=readjxms_ssht(fnames,minage,maxage,ageref,L3)
% [hmap,lo,la,age,hmapref]=readjxms_ssht(fnames,minage,maxage,ageref,L3)
% 
% INPUT:
%
% fnames        String with path of files (wildcards allowed)
% minage		minimum age of files to load
% maxage		maximum age of files to load
% ageref		age of reference datum (pass NaN for none)
% L3            An optional extra truncation for the expansions
%
% OUTPUT:
%
% hmap          The data maps or vectors appropriately referenced
% lo            Longitudes of these data, in degrees
% la            Latitudes of these data, in degrees
% age			Ages of the data
% hmapref		Reference map or vectors
% 
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu 11 Jul 2013 00:52:06 EDT

% Define defaults
defval('fnames','seagl_*.mn')
defval('minage',0)
defval('maxage',200)
defval('L3',256)
defval('ageref',0)

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
  [hlm2,dels2]=readit(fname2,L3);
  L2=max(dels2);
  [hmapref,lo,la] = sshtinverse(hlm2,L2);
end

if length(ffiles)>0

		ind=1;
		disp(['      AGE ' num2str(age(ind))]); 
		fname1=fullfile(fpath,ffiles(ind).name);
		L1=L1s(ind);
		
		% Load the data model as coefficients
		[hlm1,dels1]=readit(fname1,L3);
		L1 = max(dels1);
		[hmap1,lo2,la2] = sshtinverse(hlm1,L1);
		hmap = nan([size(hmap1) length(ffiles)]);
		hmap(:,:,ind) = hmap1-hmapref;

		parfor ind=2:length(ffiles)

			disp(['      AGE ' num2str(age(ind))]); 
			fname1=fullfile(fpath,ffiles(ind).name);
			L1=L1s(ind);
		
			% Load the data model as coefficients
			[hlm1,dels1]=readit(fname1,L3);
			L1 = max(dels1);
			hmap(:,:,ind) = sshtinverse(hlm1,L1)-hmapref;
		end
else
	hmap = hmapref;
end

if ~isstr(fname2)
	lo = lo2; la = la2;
end


%%%%%%%
function [specs]=readheader(fname)

fid = fopen(fname);
tline = fgetl(fid);
specs = sscanf(tline,'%f');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hlm,dels,dems]=readit0(fname)
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
%dom=sqrt(2-(dems==0));
dom=1;
hlm(:,1)=hlm(:,1).*CSC.*dom;
hlm(:,2)=hlm(:,2).*CSC.*dom;
% Fake arrival to align with readshds
%hlm=hlm';
%hlm=hlm(:);

%%%%%%%%
function [hlm,dels,dems]=readit(fname,Ltrunc)

defval('Ltrunc',[]);

[hlm0,dels0,dems0]=readit0(fname);

if length(Ltrunc)>0
	sub=find(dels0<=Ltrunc);
	hlm0=hlm0(sub,:); dels0=dels0(sub); dems0=dems0(sub);
end

Lmax = max(dels0);

hlm=zeros((Lmax+1)^2,1);
counter=0;
for ii=0:Lmax
	sub=find(dels0==ii);
	if length(sub)>1
		sub2=[sub(end:-1:2) ; sub(1:end)];
	else
		sub2=sub;
	end
	hlm((counter+1):(counter+length(sub2))) = hlm0(sub2,1) + i*hlm0(sub2,2);
	counter=counter+length(sub2);
end
dels=floor(sqrt([1:length(hlm)]-1));
dems=[1:length(hlm)] - 1 - dels.*dels - dels;
sub=find(dems<0);
hlm(sub)=(-1).^dems(sub)' .* conj(hlm(sub));

%%%%%
function [f,lo,la] = sshtinverse(hlm,L)

%tic;
f=ssht_inverse(hlm,L+1,'Method','MWSS','Reality',true);
%toc;

f=f*sqrt(4*pi);

[thetas,phis]=ssht_sampling(L+1,'Method','MWSS');
la=90-rad2deg(thetas);
lo=mod(rad2deg(phis)+180,360);
[lo,ind]=sort(lo);
f=f(:,ind);