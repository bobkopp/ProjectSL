function [fp,fpname,lo,la] = readFingerprint(subdir,usessht)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Jun 15 22:08:21 EDT 2014

defval('subdir','IFILES/slr/FPRINT');
defval('usessht',1);

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
        if usessht
            [fp0,lo,la] = readjxms_ssht([files(i).name(1:end-3)],0,0,NaN,L3);
        else
            [fp0,lo,la] = readjxms([files(i).name(1:end-3)],0,0,NaN,360/(2*L3),L3);
        end
        
	fp(:,:,i)=fp0;
	system(['rm /tmp/' files(i).name(1:end-3)]);	
end
cd(pd);
