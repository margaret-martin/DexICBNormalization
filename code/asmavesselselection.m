function asmavesselselection( rootPath,action )


%Channel 1 is CD3 and channel 2 is aSMA
modCD3=1;  modaSMA=2;  moddapi=3;
modCt = 3;

reset1=0;
if exist('actionParam','var'),  eval(actionParam);  end

warning off all;
visuPath = 'visu'; visuPath=fullfile(rootPath,visuPath);
if ~exist(visuPath,'dir'), mkdir(visuPath);   end
csvPath = 'csv'; csvPath=fullfile(rootPath,csvPath);
if ~exist(csvPath,'dir'), mkdir(csvPath);   end
dataPath = 'Data'; dataPath=fullfile(rootPath,dataPath);
temp=what;  code_path=temp.path;
matPath = 'mat'; matPath=fullfile(rootPath,matPath);
if ~exist(matPath,'dir'), mkdir(matPath);   end

%get cd31 bw 
if strcmp(action,'step1')
    % aSMA_areafraction_nonvessel=[];
    paths=dirPaths(dataPath);
    aSMA_vesseldensity=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        %         im = [];
        %         running_dist = [];
        %         tot_avgdist = [];
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' cd31 ') && ~contains(data_name,'asma')
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,type));
                if ~exist(filename,'file') || reset1
                    tic;
                    %                 data=readStack_mod_modCt(fullfile(path,files(f).name),modCt,modaSMA);
                    %                 data=data(:,:,4);
                    data = imread(fullfile(path,files(f).name));
                    data = data(:,:,modCD3);
                    
                    %VESSEL SElECTION BEGINS
                    %ACTIVE CONTOUR
                    mask = zeros(size(data));
                    locs = find(data>120);
                    mask(locs) = 1;
                    %mask(200:end-200,200:end-200) = 1;
                    acbw = activecontour(data,mask,1000);
                    
                    %DILATE AND ERODE
                    %se1 = strel('line',4,0);
                    %se2 = strel('line',4,90);
                    %                 se = strel('disk',4,4);
                    %                 bw = imdilate(acbw, se);
                    %                 bw = imerode(bw,se);
                    %                 %ERODE AND DILATE
                    %                 se2 = strel('disk',3,4);
                    %                 bw = imerode(bw,se2);
                    %                 bw = imdilate(bw,se2);
                    %
                    %                 %CREATE SKEL
                    %                 bwt5 = bwmorph(bw,'thin',5);
                    %                 bwtk5 = bwmorph(bwt5,'thicken',5);
                    bwsI = bwmorph(acbw,'skel',Inf);
                    %thin lines and then thicken to get smoother lines
                    %then get skel
                    
                    bw=acbw;
                    bw_skel = bwsI;
                    %BRANCH POINTS
                    %                 bwbp = bwmorph(bwsI,'branchpoints');
                    %
                    %                 bwsubtract = logical(bw_skel-bwbp);
                    %                 bw_l = bwlabel(bwsubtract);
                    %                 for i=1:max(max(bw_l))
                    %                     vessel = find(bw_l==i);
                    %                     if length(vessel)<20
                    %                         for j=1:length(vessel)
                    %                             bw_skel(vessel(j))=0;
                    %                         end
                    %                     end
                    %                 end
                    
                    bw_skel = logical(bwmorph(bw_skel,'clean'));
                    
                    
                    %get vessel density and non-aSMA area fraction
                    bwbp = bwmorph(bw_skel,'branchpoints');
                    bwsubtract = logical(bw_skel-bwbp);
                    bw_l = bwlabel(bwsubtract);
                    num_vessels = max(max(bw_l));
                    
                    nonaSMA_size=length(find(bw==0));
                    [m,n]=size(bw);
                    tot_pixels = m*n;
                    vessel_d=num_vessels/tot_pixels;
                    nonaSMA_areafrac=length(nonaSMA_size)/tot_pixels;
                    
                    %                 threshlocs = find(data<85);
                    %                 for i=1:length(threshlocs)
                    %                     data(threshlocs(i))=0;
                    %                 end
                    %                 bw = imbinarize(data);
                    %                 bw_l = bwlabel(bw);
                    %                 for i=1:max(max(bw_l))
                    %                     area = find(bw_l==i);
                    %                     if length(area)<40
                    %                         for j=1:length(area)
                    %                             bw(area(j))=0;
                    %                         end
                    %                     end
                    %                 end
                    
                    %VESSEL SELECTION ENDS
                    %SAVE
                    save(filename,'bw','bw_skel');   clear im bw;
                    tt=toc;  display(sprintf('%s : %g/%g {%g sec}',data_name,f,numel(files),tt));
                end
                %             spaces = strfind(data_name,' ');
                %
                %             curr_name = char(data_name);
                %             aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
                %             load(aSMAmatname);
                %
                %             [m,n] = size(bw);
                %             total_size=m*n;
                %             aSMA_size=length(find(bw==1));
                %
                %             fraction=aSMA_size/total_size;
                aSMA_vesseldensity=[aSMA_vesseldensity vessel_d];
                curr_name = char(data_name);
                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'aSMAvesseldensity.csv');
    aSMA_vesseldensity=transpose(aSMA_vesseldensity);
    % CD3_areafraction=transpose(CD3_areafraction);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,aSMA_vesseldensity);
    writetable(T,csvname);
    
end

%take out small parts of cd31 bw ang get total vessel density
if strcmp(action,'step2')
    paths=dirPaths(dataPath);
    totalvesseldensity=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' cd31 ') && ~contains(data_name,'asma')
                tic;
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,type));
                load(filename)
                bw_l=bwlabel(bw);
                for i=1:max(max(bw_l))
                    vessel = find(bw_l==i);
                    if length(vessel)<50
                        for j=1:length(vessel)
                            bw(vessel(j))=0;
                            bw_skel(vessel(j))=0;
                        end
                    end
                end
                
                bw_skel = logical(bwmorph(bw_skel,'clean'));
                
                
                %get vessel density 
                bwbp = bwmorph(bw_skel,'branchpoints');
                bwsubtract = logical(bw_skel-bwbp);
                bw_l = bwlabel(bwsubtract);
                num_vessels = max(max(bw_l));
                
                [m,n]=size(bw);
                tot_pixels = m*n;
                vessel_d=num_vessels/tot_pixels;
                new_filename=fullfile(matPath,sprintf('%s_%s_final.mat',data_name,type));
                
                %VESSEL SELECTION ENDS
                %SAVE
                save(new_filename,'bw','bw_skel');   clear im bw;
                tt=toc;  display(sprintf('%s : %g/%g {%g sec}',data_name,f,numel(files),tt));
                totalvesseldensity=[totalvesseldensity vessel_d];
                curr_name = char(data_name);
                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'totalvesseldensity.csv');
    totalvesseldensity=transpose(totalvesseldensity);
    % CD3_areafraction=transpose(CD3_areafraction);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,totalvesseldensity);
    writetable(T,csvname);
    
end

%get just vessel density
if strcmp(action,'totalvesseldensity')
    paths=dirPaths(dataPath);
    totalvesseldensity=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' cd31 ') && ~contains(data_name,'asma')
                tic;
                filename=fullfile(matPath,sprintf('%s_%s_final.mat',data_name,type));
                load(filename)
                
                %get vessel density 
                bwbp = bwmorph(bw_skel,'branchpoints');
                bwsubtract = logical(bw_skel-bwbp);
                bw_l = bwlabel(bwsubtract);
                num_vessels = max(max(bw_l));
                
                [m,n]=size(bw);
                tot_pixels = m*n;
                vessel_d=num_vessels/tot_pixels;
                
                %VESSEL SELECTION ENDS
                %SAVE
                curr_name = strcat(data_name,'_',type);
                curr_name = char(curr_name);
                tt=toc;  display(sprintf('%s : %g/%g {%g sec}',curr_name,f,numel(files),tt));
                totalvesseldensity=[totalvesseldensity vessel_d];
                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'totalvesseldensity.csv');
    totalvesseldensity=transpose(totalvesseldensity);
    % CD3_areafraction=transpose(CD3_areafraction);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,totalvesseldensity);
    writetable(T,csvname);
    
end

%get asma bw 
if strcmp(action,'step3')
    paths=dirPaths(dataPath);
    nonvesselareafrac_asma=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        %         im = [];
        %         running_dist = [];
        %         tot_avgdist = [];
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' asma ') && ~contains(data_name,'cd31')
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,type));
                a=strfind(data_name,'asma');
                loadname=strcat('mat/',data_name(1:a-1),' cd31',data_name(a+4:end),'_',type,'_final.mat');
                load(loadname);
                vesselbw=bw;
                data = imread(fullfile(path,files(f).name));
                data = data(:,:,modaSMA);
                
                %VESSEL SElECTION BEGINS
                %ACTIVE CONTOUR
                mask = zeros(size(data));
                locs = find(data>120);
                mask(locs) = 1;
                %mask(200:end-200,200:end-200) = 1;
                acbw = activecontour(data,mask,1000);
                bw_l=bwlabel(acbw);
                for i=1:max(max(bw_l))
                    vessel = find(bw_l==i);
                    if length(vessel)<50
                        for j=1:length(vessel)
                            acbw(vessel(j))=0;
                        end
                    end
                end
                
                bw_skel = bwmorph(acbw,'skel',Inf);
                %thin lines and then thicken to get smoother lines
                %then get skel
                
                bw=acbw;
                
                bw_skel = logical(bwmorph(bw_skel,'clean'));
                
                
                %get non-aSMA area fraction
                vessellocs=find(vesselbw==1);
                nonvesselarea=bw;
                nonvesselarea(vessellocs)=0;
                
                
                nonvessel_size=length(find(nonvesselarea==1));
                [m,n]=size(bw);
                tot_pixels = m*n;
                nonvessel_areafrac=length(nonvessel_size)/tot_pixels;
                
                %                 threshlocs = find(data<85);
                %                 for i=1:length(threshlocs)
                %                     data(threshlocs(i))=0;
                %                 end
                %                 bw = imbinarize(data);
                %                 bw_l = bwlabel(bw);
                %                 for i=1:max(max(bw_l))
                %                     area = find(bw_l==i);
                %                     if length(area)<40
                %                         for j=1:length(area)
                %                             bw(area(j))=0;
                %                         end
                %                     end
                %                 end
                
                %VESSEL SELECTION ENDS
                %SAVE
                save(filename,'bw','bw_skel');   clear im bw;
                display(sprintf('%s : %g/%g',data_name,f,numel(files)));
                
                %             spaces = strfind(data_name,' ');
                %
                %             curr_name = char(data_name);
                %             aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
                %             load(aSMAmatname);
                %
                %             [m,n] = size(bw);
                %             total_size=m*n;
                %             aSMA_size=length(find(bw==1));
                %
                %             fraction=aSMA_size/total_size;
                nonvesselareafrac_asma=[nonvesselareafrac_asma nonvessel_areafrac];
                curr_name = char(data_name);
                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'nonvessel_aSMAareafraction.csv');
    nonvesselareafrac_asma=transpose(nonvesselareafrac_asma);
    % CD3_areafraction=transpose(CD3_areafraction);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,nonvesselareafrac_asma);
    writetable(T,csvname);

end

%get nonvessel asma area fraction
if strcmp(action,'step4')
    paths=dirPaths(dataPath);
    nonvesselareafrac_asma=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        %         im = [];
        %         running_dist = [];
        %         tot_avgdist = [];
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' asma ') && ~contains(data_name,'cd31')
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,type));
                a=strfind(data_name,'asma');
                loadname=strcat('mat/',data_name(1:a-1),' cd31',data_name(a+4:end),'_',type,'_final.mat');
                load(loadname);
                vesselbw=bw;
                load(filename);
                nonvesselarea=bw;
                
                
                %get non-aSMA area fraction
                vessellocs=find(vesselbw==1);
                for i=1:length(vessellocs)
                    nonvesselarea(vessellocs(i))=0;
                end
                
                nonvessel_size=length(find(nonvesselarea==1));
                [m,n]=size(bw);
                tot_pixels = m*n;
                nonvessel_areafrac=nonvessel_size/tot_pixels;
                
                
                %VESSEL SELECTION ENDS
                %SAVE
                %                 save(filename,'bw','bw_skel');   clear im bw;
                curr_name = strcat(data_name,'_',type);
                curr_name = char(curr_name);
                display(sprintf('%s : %g/%g',curr_name,f,numel(files)));
                
                nonvesselareafrac_asma=[nonvesselareafrac_asma nonvessel_areafrac];
                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'nonvessel_aSMAareafraction.csv');
    nonvesselareafrac_asma=transpose(nonvesselareafrac_asma);
    % CD3_areafraction=transpose(CD3_areafraction);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,nonvesselareafrac_asma);
    writetable(T,csvname);

end


%aSMA+ vessel density
if strcmp(action,'step5')
    thresh=0.1;
    paths=dirPaths(dataPath);
    asma_vesseldensity=[];
    asma_vesselfrac=[];
    name = [];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        %         im = [];
        %         running_dist = [];
        %         tot_avgdist = [];
        cut=strfind(paths(p).name,'_');
        type=paths(p).name(1:cut(1)-1);
        
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' asma ') && ~contains(data_name,'cd31')
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,type));
                a=strfind(data_name,'asma');
                loadname=strcat('mat/',data_name(1:a-1),' cd31',data_name(a+4:end),'_',type,'_final.mat');
                load(loadname);
                vesselbw=bw;
                vesselskel=bw_skel;
                load(filename);
                asmabw=bw;
                
                bwbp = bwmorph(vesselskel,'branchpoints');
                bwsubtract = logical(vesselskel-bwbp);
                bw_l = bwlabel(bwsubtract);
                num_vessels = max(max(bw_l));
                
                vessel_count=0;
                for i=1:num_vessels
                    locs=find(bw_l==i);
                    check=zeros(length(locs),1);
                    for j=1:length(locs)
                        check(j)=asmabw(locs(j));
                    end
                    ratio=length(find(check==1))/length(check);
                    if ratio>=thresh
                        vessel_count=vessel_count+1;
                    end
                    
                end
                
                
                [m,n]=size(bw);
                tot_pixels = m*n;
                curr_density=vessel_count/tot_pixels;
                curr_frac=vessel_count/num_vessels;
                
                
%                 %SAVE
%                 save(filename,'bw','bw_skel');   clear im bw;
                
                asma_vesseldensity=[asma_vesseldensity curr_density];
                asma_vesselfrac=[asma_vesselfrac curr_frac];
                curr_name = strcat(data_name,'_',type);
                curr_name = char(curr_name);
                display(sprintf('%s : %g/%g',curr_name,f,numel(files)));

                %             filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                %             name = [name;string(filename)];
                name = [name;string(curr_name)];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'aSMA+vessels.csv');
    asma_vesseldensity=transpose(asma_vesseldensity);
    asma_vesselfrac=transpose(asma_vesselfrac);
    % aSMA_celldensity=transpose(aSMA_celldensity);
    % CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,asma_vesseldensity,asma_vesselfrac);
    writetable(T,csvname);

end



end
function paths=dirPaths(path)
paths=dir(path);    paths=paths(3:end);
paths=paths([paths.isdir]);
end

function im = readStack_mod_modCt(fn,modCt,mod2read)
inf = imfinfo(fn);
zCt=numel(inf)./modCt;
for (i=1:zCt)
    im(:,:,i) = imread(fn, i + ((mod2read-1).*zCt));
end
end

function files = dirTiff(path)
files=dir(fullfile(path,'*.TIF'));
if isempty(files),    files=dir(fullfile(path,'*.tif'));    end
if isempty(files),    files=dir(fullfile(path,'*.tiff'));   end
if isempty(files),    files=dir(fullfile(path,'*.TIFF'));   end
end

