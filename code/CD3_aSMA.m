function CD3_aSMAv20190103(rootPath,action,actionParam)


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

paths=dirPaths(dataPath);

if strcmp(action,'area fraction')
%     csvname = fullfile(csvPath,'aSMA_area_fraction.csv');
%     fid = fopen(csvname, 'wt');
%     fprintf(fid,'Name,aSMA Area Fraction\n');
    %visu(:,:,3)=aSMAbw
    %visu(:,:,2)=CD3bw
    aSMA_areafraction=[];
    CD3_areafraction=[];
    aSMA_celldensity=[];
    CD3_celldensity=[];
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
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,' aSMA ') && ~contains(data_name,'CD3-aSMA')
                spaces = strfind(data_name,' ');
                curr_name = char(data_name);
                curr_name = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end));
                aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
                load(aSMAmatname);
                
                [m,n] = size(bw);
                total_size=m*n;
                aSMA_size=length(find(bw==1));
                
                fraction=aSMA_size/total_size;
                aSMA_areafraction=[aSMA_areafraction fraction];
                curr_name = char(data_name);
                filename = strcat(curr_name(1:spaces(1)),curr_name(spaces(end-1):end));
                name = [name;string(filename)];
                %fprintf(fid,'%s,%g\n',files(f).name,fraction);
%             elseif contains(data_name,' a-SMA ') && ~contains(data_name,'CD3-aSMA')
%                 spaces = strfind(data_name,' ');
%                 curr_name = char(data_name);
%                 curr_name = strcat(curr_name(1:spaces(1)),' a-SMA',curr_name(spaces(2):end));
%                 aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
%                 load(aSMAmatname);
%                 
%                 [m,n] = size(bw);
%                 total_size=m*n;
%                 aSMA_size=length(find(bw==1));
%                 
%                 fraction=aSMA_size/total_size;
%                 aSMA_areafraction=[aSMA_areafraction fraction];
                %fprintf(fid,'%s,%g\n',files(f).name,fraction);
            elseif contains(data_name,' CD3 ') && ~contains(data_name,'CD3-aSMA')
                spaces = strfind(data_name,' ');
                curr_name = char(data_name);
                curr_name = strcat(curr_name(1:spaces(1)),' CD3',curr_name(spaces(2):end));
                CD3matname = fullfile(matPath,sprintf('%s_CD3.mat',curr_name));
                load(CD3matname);
                
                [m,n] = size(bw);
                total_size=m*n;
                CD3_size=length(find(bw==1));
                
                fraction=CD3_size/total_size;
                CD3_areafraction=[CD3_areafraction fraction];
            elseif contains(data_name, 'dapi')
                visuname = fullfile(visuPath,data_name);
                curr=load(visuname);
                curr=curr.visu;
                [L,num]=bwlabel(curr(:,:,2));
                [m,n] = size(curr(:,:,2));
                total_size=m*n;
                d= num/total_size;
                CD3_celldensity=[CD3_celldensity d];
                [L,num]=bwlabel(curr(:,:,3));
                d= num/total_size;
                aSMA_celldensity=[aSMA_celldensity d];
            end
            
        end
    end
    
    csvname = fullfile(csvPath,'areafraction_celldensity.csv');
    aSMA_areafraction=transpose(aSMA_areafraction);
    CD3_areafraction=transpose(CD3_areafraction);
    aSMA_celldensity=transpose(aSMA_celldensity);
    CD3_celldensity=transpose(CD3_celldensity);
    T = table(name,aSMA_areafraction,CD3_areafraction,aSMA_celldensity,CD3_celldensity);
    writetable(T,csvname);
    
end


%aSMA
if strcmp(action,'aSMA')
    param.thresh=70;   param.holes=1;  param.minSize=50;   cTable=[200 1450];
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
%         files=dirTiff(path);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,'aSMA') && ~contains(data_name,'CD3-aSMA')
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,action));
                if ~exist(filename,'file') || reset1
                    tic;
                    %                 data=readStack_mod_modCt(fullfile(path,files(f).name),modCt,modaSMA);
                    %                 data=data(:,:,4);
                    data = imread(fullfile(path,files(f).name));
                    data = data(:,:,modaSMA);
                    %                 visu=zeros(size(data,1),size(data,2),3);
                    %                 visu(:,:,2)=data./cTable(2);
                    %                 visu(visu>1)=1;
                    %                 [bw,param]=manualProcess_2D_1(data,visu,param);
                    threshlocs = find(data<85);
                    for i=1:length(threshlocs)
                        data(threshlocs(i))=0;
                    end
                    bw = imbinarize(data);
                    bw_l = bwlabel(bw);
                    for i=1:max(max(bw_l))
                        area = find(bw_l==i);
                        if length(area)<40
                            for j=1:length(area)
                                bw(area(j))=0;
                            end
                        end
                    end
                    save(filename,'bw','param');   clear im bw;
                    tt=toc;  display(sprintf('%s %s : %g/%g {%g sec}',action,data_name,f,numel(files),tt));
                end
            end
        end
    end
end


if strcmp(action,'CD3')
    param.thresh=190;   param.holes=1;  param.minSize=2;   cTable=[200 1450]; %20170317 updated threshold from 1000
    for p=1:numel(paths)
        path=fullfile(dataPath,paths(p).name);
        %files=dirTiff(path);
%         files = path(end-1:end);
%         if files(1)=='_'
%             files = files(end);
%         end
%         files = str2num(files);
        files=dir(fullfile(path));
        t = length(files);
        num_files = path(end-1:end);
        if num_files(1)=='_'
            num_files = num_files(end);
        end
        num_files = str2num(num_files);
        extra = t - num_files +1;
        files=files(extra:end);
        %MyFolderInfo = dir('myfolder')
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if ~isempty(strfind(data_name,'CD3')) && isempty(strfind(data_name,'CD3-aSMA'))
                filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,action));
                if ~exist(filename,'file') || reset1
                    tic;
                    %data=readStack_mod_modCt(fullfile(path,files(f).name),modCt,modCD3);
                    %data=data(:,:,4);
                    data = imread(fullfile(path,files(f).name));
                    data = data(:,:,modCD3);
                    %visu=zeros(size(data,1),size(data,2),3);
                    %visu(:,:,2)=data./cTable(2);
                    %visu(visu>1)=1;
                    %[bw,param]=manualProcess_2D_1(data,visu,param);
                    threshlocs = find(data<175);
                    for i=1:length(threshlocs)
                        data(threshlocs(i))=0;
                    end
                    bw = imbinarize(data);
                    
%                     se = strel('disk',3,6);
%                     bw = imerode(bw,se);
%                     bw = imdilate(bw,se);
                    bw_l = bwlabel(bw);
                    for i=1:max(max(bw_l))
                        area = find(bw_l==i);
                        if length(area)<15
                            for j=1:length(area)
                                bw(area(j))=0;
                            end
                        end
                    end
                    % Create masked image.
                    %maskedImage = X;
                    %maskedImage(~bw) = 0;
                    se = strel('disk',2,4);
                    bw = imerode(bw,se);
                    bw = imdilate(bw,se);
                    save(filename,'bw');   clear im bw;
                    tt=toc;  display(sprintf('%s %s : %g/%g {%g sec}',action,data_name,f,numel(files),tt));
                end
            end
        end
    end
end

if strcmp(action,'dapi')
    param.thresh=60;   param.holes=1;  param.minSize=50;   cTable=[200 1450];
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
%         files=dirTiff(path);
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            filename=fullfile(matPath,sprintf('%s_%s.mat',data_name,action));
            if ~exist(filename,'file') || reset1
                tic;
%                 data=readStack_mod_modCt(fullfile(path,files(f).name),modCt,moddapi);
%                 data=data(:,:,4);
%                 visu=zeros(size(data,1),size(data,2),3);
%                 visu(:,:,2)=data./cTable(2);
%                 visu(visu>1)=1;
%                 [bw,param]=manualProcess_2D_1(data,visu,param);
                data = imread(fullfile(path,files(f).name));
                data = data(:,:,moddapi);
                bw = imbinarize(data);
                save(filename,'bw','param');   clear im bw;
                tt=toc;  display(sprintf('%s %s : %g/%g {%g sec}',action,data_name,f,numel(files),tt));
            end
        end
    end
end

if strcmp(action,'overlay')
    csvname = fullfile(csvPath,'analyze.csv');
    fid = fopen(csvname, 'wt');
    fprintf(fid,'Name,Average Distance\n');
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
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            %if contains(data_name,'CD3-aSMA') || contains(data_name,'CD3_a-SMA')
            %if contains(data_name,'CD3-aSMA')
            if contains(data_name,'dapi')
                spaces = strfind(data_name,' ');
                curr_name = char(data_name);
                curr_name = strcat(curr_name(1:spaces(1)),' CD3',curr_name(spaces(2):end));
                CD3matname = fullfile(matPath,sprintf('%s_CD3.mat',curr_name));
                load(CD3matname);
                CD3bw = bw;
                spaces = strfind(data_name,' ');
                curr_name = char(data_name);
                if contains(data_name,'CD3-aSMA')
                    aSMAbmpname = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end),'.bmp');
                    curr_name = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end));
                    
                elseif contains(data_name,'CD3_a-SMA')
                    aSMAbmpname = strcat(curr_name(1:spaces(1)),' a-SMA',curr_name(spaces(2):end),'.bmp');
                    curr_name = strcat(curr_name(1:spaces(1)),' a-SMA',curr_name(spaces(2):end));
                    
                end
                %curr_name = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end));
                aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
                load(aSMAmatname);
                aSMAbw = bw;
                visu=zeros(size(aSMAbw,1),size(aSMAbw,2),3);
                %visu(:,:,1)=CD31bw.*255;
                %visu(:,:,2)=CD3bw.*255;
                %visu(:,:,3)=CA9bw.*255;
                
                visu(:,:,3)=aSMAbw.*255;
                visu(:,:,2)=CD3bw.*255;
                %visu(:,:,1)=aSMAbw.*255;
                
                visufile = fullfile(visuPath,data_name);
                save(visufile,'visu');
                
                %changed to visu 2
                [L,num]=bwlabel(visu(:,:,2));
                %find centroid
                s = regionprops(L,'centroid');
                centroids = cat(1, s.Centroid);
                centroids = int64(centroids);
                
                %find number of pixels in region
                j = regionprops(L,'Area');
                CD3_areas = cat(1, j.Area);
                CD3_areas = int64(CD3_areas);
                
                %             CD31_distmap = bwdist(visu(:,:,2));
                %             CD3_distmap = bwdist(visu(:,:,3));
                CD3_distmap = bwdist(visu(:,:,2));
                aSMA_distmap = bwdist(visu(:,:,3));
                
                dist=[];
                for x = 1:num
                    dist(x) = aSMA_distmap(centroids(x,2),centroids(x,1));
                end
                
                
                csvname_cell = fullfile(csvPath,sprintf('%s_analyze_percell.csv', files(f).name));  %sprintf('%s_CD31.mat',data_name)
                fid_cell = fopen(csvname_cell, 'wt');
                
                %T = table(number,dist_to_vessel,dist_to_hyp, size_pix);
                T = table(num,dist);
                %T = table(num_hyp_cell,avg_hypcell_vessel_dist_array,avg_hypcell_hypoxia_dist_array,hyp_cell_size);%,num_norm_cell,avg_normcell_vessel_dist_array,avg_normcell_hypoxia_dist_array,norm_cell_size);
                writetable(T,csvname_cell);
                
                %fprintf(fid_cell,'%d',round(num_hyp_cell));
                %fprintf(fid_cell,'%d',avg_hypcell_vessel_dist_array);
                %fprintf(fid_cell,'%d',avg_hypcell_hypoxia_dist_array);
                %fprintf(fid_cell,'%d',hyp_cell_size);
                
                %fprintf(fid_cell,'%d',round(num_norm_cell));
                %fprintf(fid_cell,'%d',avg_normcell_vessel_dist_array);
                %fprintf(fid_cell,'%d',avg_normcell_hypoxia_dist_array);
                %fprintf(fid_cell,'%d',norm_cell_size);
                
                %             clear csvname_cell fid_cell num_hyp_cell num_norm_cell T
                %             clear number dist_to_vessel dist_to_hyp size_pix
                %
                %             clear dist_to_norm cell_location cell_location_hyp cell_location_norm %%%
                %             clear dist_to_Tcell dist_to_noTcell%%%
                
                %             avg_hypcell_vessel_dist = mean(avg_hypcell_vessel_dist_array);
                %avg_hypcell_hypoxia_dist = mean(avg_hypcell_hypoxia_dist_array);
                %             avg_normcell_vessel_dist = mean(avg_normcell_vessel_dist_array);
                %             avg_normcell_hypoxia_dist = mean(avg_normcell_hypoxia_dist_array);
                
                
                
                avgdist = mean(dist);
                fprintf(fid,'%s,%g\n',files(f).name,avgdist);
                %RUNNING AVG ATTEMPT
%                 space = strfind(data_name,' ');
%                 space = space(1);
%                 im_name = char(data_name);
%                 im_name = im_name(1:space-1);
%                 im = [im im_name];
%                 tot_avgdist = [tot_avgdist avgdist];
                
                
            end
%             clear avg_hypcell_vessel_dist_array avg_normcell_vessel_dist_array avg_normcell_hypoxia_dist_array avg_hypcell_hypoxia_dist_array hyp_cell_size norm_cell_size
%             clear avg_hypcell_normoxia_dist_array avg_normcell_hypoxia_dist_array avg_normcell_normoxia_dist_array normoxia_bw normoxia_distmap%%%
%             space = strfind(data_name,' ');
%             space = space(1);
%             im_name = char(data_name);
%             im_name = im_name(1:space-1);
%             running_dist = [running_dist avgdist];
%             if f+1 <= numel(files)
%                 new_data_name = files(f+1).name;
%                 space = strfind(a,' ');
%                 space = space(1);
%                 new_im_name = char(new_data_name);
%                 new_im_name = new_im_name(1:space-1);
%                 if im_name ~=new_im_name
%                     im = [im im_name];
%                     tot_avgdist = [tot_avgdist mean(running_dist)];
%                     running_dist = [];
%                 end
%             else
%                 im = [im im_name];
%                 tot_avgdist = [tot_avgdist mean(running_dist)];
%             end
            
        end
        
        
        
    end
    
    %RUNNING AVG ATTEMPT
%     fprintf(fid,'\n\nType,Average Distance\n');
%     running_avg = [];
%     curr=im(1);
%     for k=1:length(im)
%         if im(k) == curr
%             running_avg = [running_avg tot_avgdist(k)];
%         else
%             fprintf(fid,'%s,%g\n',im(k-1),mean(running_avg));
%             curr = im(k);
%             running_avg = [running_avg];
%         end
%         
%     end
%     fprintf(fid,'%s,%g\n',im(end),mean(running_avg));
    
end


if strcmp(action,'thresh_overlay')
    %aSMAthresh=0.0435;
    aSMAthresh=0.02;
    %aSMA_fracs=csvread('csv/aSMA_area_fraction.csv');
    %filenames=aSMA_fracs(:,1);
    readcsvname = fullfile(csvPath,'aSMA_area_fraction.csv');
    fid2 = fopen(readcsvname, 'r');
    C=textscan(fid2,'%s%s',1,'delimiter',',');
    C=textscan(fid2,'%s%f','delimiter',',');
    fclose(fid2);
    aSMA_names=C{:,1};
    aSMA_fracs=C{:,2};
    csvname = fullfile(csvPath,'analyze.csv');
    fid = fopen(csvname, 'wt');
    fprintf(fid,'Name,Average Distance\n');
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
        for f=1:numel(files)
            data_name = strrep(files(f).name,'.bmp','');
            if contains(data_name,'dapi')
                spaces = strfind(data_name,' ');
                curr_name = char(data_name);
                if contains(data_name,'CD3-aSMA')
                    aSMAbmpname = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end),'.bmp');
                    curr_name = strcat(curr_name(1:spaces(1)),' aSMA',curr_name(spaces(2):end));
                    
                elseif contains(data_name,'CD3_a-SMA')
                    aSMAbmpname = strcat(curr_name(1:spaces(1)),' a-SMA',curr_name(spaces(2):end),'.bmp');
                    curr_name = strcat(curr_name(1:spaces(1)),' a-SMA',curr_name(spaces(2):end));
                    
                end
                idx = find(strcmp(aSMA_names(:), aSMAbmpname));
                %idx = strcmp(aSMA_names(:), aSMAbmpname);
                if aSMA_fracs(idx)>=aSMAthresh
                    
                    aSMAmatname = fullfile(matPath,sprintf('%s_aSMA.mat',curr_name));
                    load(aSMAmatname);
                    aSMAbw = bw;
                    
                    spaces = strfind(data_name,' ');
                    curr_name = char(data_name);
                    curr_name = strcat(curr_name(1:spaces(1)),' CD3',curr_name(spaces(2):end));
                    CD3matname = fullfile(matPath,sprintf('%s_CD3.mat',curr_name));
                    load(CD3matname);
                    CD3bw = bw;
                    
                    visu=zeros(size(aSMAbw,1),size(aSMAbw,2),3);
                    
                    visu(:,:,3)=aSMAbw.*255;
                    visu(:,:,2)=CD3bw.*255;
                    
                    visufile = fullfile(visuPath,data_name);
                    save(visufile,'visu');
                    
                    %changed to visu 2
                    [L,num]=bwlabel(visu(:,:,2));
                    %find centroid
                    s = regionprops(L,'centroid');
                    centroids = cat(1, s.Centroid);
                    centroids = int64(centroids);
                    
                    %find number of pixels in region
                    j = regionprops(L,'Area');
                    CD3_areas = cat(1, j.Area);
                    CD3_areas = int64(CD3_areas);
                    
                    %             CD31_distmap = bwdist(visu(:,:,2));
                    %             CD3_distmap = bwdist(visu(:,:,3));
                    CD3_distmap = bwdist(visu(:,:,2));
                    aSMA_distmap = bwdist(visu(:,:,3));
                    
                    dist=[];
                    for x = 1:num
                        dist(x) = aSMA_distmap(centroids(x,2),centroids(x,1));
                    end
                    
                    
                    csvname_cell = fullfile(csvPath,sprintf('%s_analyze_percell.csv', files(f).name));  %sprintf('%s_CD31.mat',data_name)
                    fid_cell = fopen(csvname_cell, 'wt');
                    
                    %T = table(number,dist_to_vessel,dist_to_hyp, size_pix);
                    T = table(num,dist);
                    %T = table(num_hyp_cell,avg_hypcell_vessel_dist_array,avg_hypcell_hypoxia_dist_array,hyp_cell_size);%,num_norm_cell,avg_normcell_vessel_dist_array,avg_normcell_hypoxia_dist_array,norm_cell_size);
                    writetable(T,csvname_cell);
                    
                    
                    avgdist = mean(dist);
                    fprintf(fid,'%s,%g\n',files(f).name,avgdist);
                else
                    
                    fprintf(fid,'%s\n',files(f).name);
                end
            end
   
            
        end
        
        
        
    end
    
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