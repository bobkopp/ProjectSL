function thetGLR=runTrainGPSLModel(coastlinestorun,thetG,ROOTDIR);

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 10 23:27:38 EST 2014

    mf=mfilename('fullpath'); mfsl=strfind(mf,'/'); mypath=mf(1:mfsl(end)-1);
    addpath(mypath,fullfile(mypath,'../lib/MFILES'),fullfile(mypath,'MFILES'),fullfile(mypath,'MFILES/scripts'));

    defval('ROOTDIR',mypath);
    IFILES=fullfile(ROOTDIR,'../../IFILES/slr/');
    PARAMDIR=fullfile(ROOTDIR,'PARAMS/');
    OUTDIR=[pwd '/'];

    coastlines=importdata(fullfile(PARAMDIR,'coastlines.txt'));
    giafile=fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc');

    psmsldir=fullfile(IFILES,'rlr_annual');
    gslfile=fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv');

    defval('coastlinestorun',1:size(coastlines.data,1));
    defval('thetG',[]);

    %thetG=[1.621 11.966 40.225 0.549];

    giamodel.gia=ncread(giafile,'Dsea_250');
    giamodel.lat=ncread(giafile,'Lat');
    giamodel.long=ncread(giafile,'Lon');


    for i=coastlinestorun
        disp(coastlines.textdata{i+1});
        [X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen]=ReadPSMSLData(coastlines.data(i,1),coastlines.data(i,2),15,psmsldir,gslfile,[],coastlines.data(i,3:end));

        sub = find(regionsu>0);
        if length(sub)>0
            if length(sub)>5
                minlength=quantile(sitelen(sub),.5);
            else
                minlength=15;
            end

            [thetGLRtemp,cvfuncGLR,thetGLRA] = TrainGPSLModel(X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen,minlength,giamodel,thetG);
            thetGLR(i,:) = thetGLRtemp;
            thetG = thetGLRtemp(1:4);
        end
    
        fid=fopen(fullfile(OUTDIR,['SLModel_' coastlines.textdata{i+1} '.txt']),'w');
        fprintf(fid,'%0.4f ',thetGLR(i,:));
        fclose(fid);
        
    end

end


