% outputSLRProjections_table_GIAhyperparameters
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 24 08:27:49 EST 2014

defval('PARAMDIR','PARAMS/');
defval('COASTLINESFILE',[PARAMDIR 'coastlines.txt']);
coastlines=importdata(COASTLINESFILE);

for i=1:size(coastlines.data,1)
	fid=fopen(fullfile(PARAMDIR,['SLModel_' coastlines.textdata{i+1} '.txt']),'r');
	thetGLRtemp=fscanf(fid,'%f');
	thetGLR(i,:)=thetGLRtemp';
	fclose(fid); 
end
M=thetGLR;
M(:,end+1) = 50*sqrt(thetGLR(:,1).^2+thetGLR(:,5).^2);

fid=fopen('hyperparameters.tsv','w');
fprintf(fid,'\t\\theta_{%0.0f}',1:size(thetGLR,2));
fprintf(fid,['\t\\theta_\\Delta']);
if exist('finescale','var')
    M(:,end+1) = finescale;
    fprintf(fid,['\t\\theta_6' char(39)]);
end
fprintf(fid,'\n');

for i=1:size(coastlines.data,1)
	fprintf(fid,[coastlines.textdata{i+1,1}]);
	fprintf(fid,'\t%0.2f',M(i,:));
	fprintf(fid,'\n');
end
fclose(fid);
