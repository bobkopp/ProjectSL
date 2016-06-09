function [fp,fpname,lo,la] = readFingerprint(subdir)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Nov 22 09:24:39 EST 2013

defval('subdir','IFILES/slr/FPRINT');

L3=512;

pd=pwd;
cd(subdir);
files=dir('*.mn.gz');
for i=1:length(files)
	system(['cp ' files(i).name ' /tmp']);
end
cd('/tmp');
clear fp fpname;
for i=1:length(files)
	s=regexp(files(i).name,'\_(.+)\.mn','tokens');
	fpname{i}=s{1}{1};
	disp(fpname(i));
	system(['gzip -d ' files(i).name ]);
	[fp0,lo,la] = readjxms_ssht([files(i).name(1:end-3)],0,0,NaN,L3);
	fp(:,:,i)=fp0;
	system(['rm /tmp/' files(i).name(1:end-3)]);	
end
cd(pd);
