%  This program is for downloading the metabolic models from BIGG and preprocess the data.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
datapath = strcat(pwd,'/data/');       %path of the data folder
for runthis = 1:0
    %% to download universal reaction names
    system('python data/download_unr.py');             %first run "data/download_unr.py" to download the universal reaction names
    unr = loadjson('data/unr.txt');         %load the universal reactions downloaded from BIGG, which is in json format
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
    system('python data/process_unrnames.py');        %Be Careful! This will take several hours! It loads unrnames.txt and download corresponding reaction's information into unrmet.txt
    unrmet = loadjson(strcat(datapath,'unrmet.txt'));  %load universal reactions' metabolites from "unrmet.txt"
    save data/unrmet.mat unrmet;     %save to mat format to save reading time in the future
    delete('data/unrmet.txt');
end

