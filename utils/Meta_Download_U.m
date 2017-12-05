%  This program is for downloading the universal metabolic reactions from BIGG and preprocess the data.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
datapath = strcat('../data/');       % path of the data folder
for runthis = 1:1
    %% to download universal reaction names
    system('python download_unr.py');             %first run "data/download_unr.py" to download the universal reaction names
    unr = loadjson('unr.txt');         %load the universal reactions downloaded from BIGG, which is in json format
    %% to save all reaction names of unr into unrnames
    fid=fopen(strcat(datapath,'unrnames.txt'),'w+');
    for i = 1:length(unr)-1
        rname = unr{i}.bigg_id;
        fprintf(fid,strcat(rname,'\r\n'));
    end
    rname = unr{i+1}.bigg_id;
    fprintf(fid,rname);
    fclose(fid);
    %% to download all universal reactions' detailed information from BIGG
    system('python data/process_unrnames.py');        % Be Careful. This will take several hours! It loads unrnames.txt and download corresponding reaction's information into unrmet.txt
    unrmet = loadjson(strcat(datapath,'unrmet.txt'));  % load universal reactions' metabolites from "unrmet.txt"
    save unrmet.mat unrmet;     % save to mat format to save reading time in the future
end

