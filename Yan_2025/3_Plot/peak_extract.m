% 用峰值index读每个站点的值
result_path = '';
peak_path = '';
site_path = ''
save_path = ''
sub = 'AdultSubjects'
LR_all = ["LH","RH"];
% model_all = ["M1_Dx_Meta","M4_Recur_PostHoc_FirstEpisodeVsHC','M4_Recur_PostHoc_RecurrentVsHC','M5_Rem_PostHoc_AcutelyDepressedVsHC','M9_HDRS','M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree','M10_Antidepressant_PostHoc_AntidepressantUseVsHC']
model_all = ["M7_AOGroup_PostHoc_AdultAOVsHC"]
for i = 1:numel(LR_all)
    for j = 1:numel(model_all)
        LR = LR_all(i)
        model = model_all(j)
        MetaForR = load(strcat(result_path,sub,'_',LR,'/',model,'_ForR.mat'));

        TFiles = MetaForR.TFiles;
        L = numel(TFiles); 
        allLH = readcell(strcat(peak_path,model,'_',LR,'_result_processed.csv'), 'Delimiter', ',');

        % allLH = readmatrix([peak_path,sub,'_',LR,'result.csv'], 'Delimiter', ',');

        M = size(allLH, 1);
        N = size(allLH, 2);

        result = cell(L+1, M);
        % M = numel(csvread([peak_path,sub,'_',LR,'result.csv'], 0, 5)); 
        % result = cell(L, M);
        result(1, 2:M) = allLH(2:M,1)';


        for i = 2:L+1
            pathParts = strsplit(TFiles{i-1}, '/');
            data_path_enigma = strcat(site_path,pathParts{13},'/StatsFixModel7/',sub,'/AnatSurf',LR,'/Thickness/',erase(model,'_Meta'),'.gii');
            data_path_direct = strcat(site_path,pathParts{13},'/Stats/',sub,'/AnatSurf',LR,'/Thickness/','/fsaverage/',erase(model,'_Meta'),'.gii');
            if exist(data_path_enigma,'file')
                data_path = data_path_enigma;
            elseif exist(data_path_direct,'file')
                data_path = data_path_direct;
            end
            % filePath = fullfile(site_path, pathParts{12});

            fileData = y_ReadAll(char(data_path));

            %indices = csvread('All_LH.csv', 0, 2);
            result{i, 1} = pathParts{13};

            for j = 2:M
                a = allLH(j,5);
                result{i, j} = fileData(a{1},:);
            end
        end

        resultCell = cell2table(result);
        writetable(resultCell, strcat(save_path,model,'_',LR,'extract.csv'), 'Delimiter', ',');
    end
end
