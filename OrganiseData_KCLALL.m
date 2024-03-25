clear all
%% GET DATA STRUCTURE
workdir='/Users/e410377/Desktop/Ludo/AlexLudo/ReorganiseKdata';

scriptdir = workdir;
datadir = [workdir, '/Data/ProjectAI_ToTransfer/'];

%% MAKE OUTDIR 
outfolder = [workdir, '/Data/', '/reformatted/'];

% Check if the directory exists
if exist(outfolder, 'dir') ~= 7
    % If it doesn't exist, create the directory
    mkdir(outfolder);
    disp(['Directory ', outfolder, ' created successfully.']);
else
    % If it already exists, display a message
    disp(['Directory ', outfolder, ' already exists.']);
end

%%
patient_data=[outfolder, '/patient_data'];
% Call the function with the output directory
createDataStructure(patient_data);
% Get a list of patient_data
exportDirectoriesToWorkspace(patient_data);
% Get a list of directories inside datadir
exportDirectoriesToWorkspace(datadir)

%% LOAD PATIENT INFO
excelFilePath = fullfile(PatientsInfoDir, 'datasetPBR_included_PrjAI_new.xlsx');
patientsInfo= importdata(excelFilePath);
patientsInfoTable=readtable(excelFilePath);
pairedInfoTable = patientsInfo.textdata.RepeatedScan;
%%
% Extract variable names and values
variableNames = patientsInfoTable.Properties.VariableNames;

% Create variables in the workspace
for index = 1:length(variableNames)
    variableName = variableNames{index};
    variableValues = patientsInfoTable.(variableName);
    
    % Create a variable in the base workspace
    assignin('base', variableName, variableValues);
end
%%
affinity_bin=ones(size(Affinity)); var=(Affinity=="MAB"); affinity_bin(var==1)=2;
sex_bin=ones(size(Gender)); var=(Gender=="female"); sex_bin(var==1)=2;


%% NOW WE ORGANISE THE DATA ..
%Weight=str2double(Weight);
CorrectiveFactor = (Weight) ./ Dose;
PtStatusB = ones(size(PrjAI_ID));  hcIndices = contains(PrjAI_ID, "HC"); PtStatusB(hcIndices) = 0;
ScanNumber=zeros(size(PrjAI_ID)); matches = ismember(PrjAI_ID, pairedInfoTable(:, 2)); ScanNumber(matches) = 1;

FR_blood = {};
FR_blood_SUV = {};
LR_blood = [];
%%
% Loop through each ID in PrjAI_ID
for idIdx = 1:length(PrjAI_ID)
    if idIdx == 32 || idIdx == 49 || idIdx == 80 || idIdx == 85  %32 : 25 frames, 49, 27 frames; 80:25 frames;  85  has weird mid frames (shifted of 0.5)
        continue;
    else
        % Construct the file name
        matFileName = [PrjAI_ID{idIdx}, '.mat'];

        % Load the mat file
        matFilePath = fullfile(TimesPET, matFileName);
        load(matFilePath, 'Times');

        % Extract midFramesTimes
        midFrames = Times.midFrameTimes;
    
        filename = sprintf('%d', idIdx - 1);

        % Save midFrames in PETmidframe directory as a text file
        dlmwrite(fullfile(PETmidframe, [filename, '.txt']), midFrames);

        % Initialize matrix for start and end frames
        FrameStartEnd = zeros(size(midFrames, 1), 2);
        FrameStartEnd(1, 1) = 0;
        FrameStartEnd(1, 2) = 2 * midFrames(1);

        % Calculate start and end frames
        for k = 2:size(midFrames, 1)
            FrameStartEnd(k, 1) = FrameStartEnd(k - 1, 2);
            FrameStartEnd(k, 2) = FrameStartEnd(k - 1, 2) + 2 * (midFrames(k) - FrameStartEnd(k - 1, 2));
        end

        % Save start and end frames in PETframestartstop directory as a text file
        dlmwrite(fullfile(PETframestartstop, [filename, '.txt']), FrameStartEnd);

        % Calculate frame durations
        FrameDur = FrameStartEnd(:, 2) - FrameStartEnd(:, 1);

        % Save frame durations in PETframelength directory as a text file
        dlmwrite(fullfile(PETframelength, [filename, '.txt']), FrameDur);

        % Load blood data
        bloodvar = fullfile(BloodDir, matFileName);
        load(bloodvar, 'Blood');

        % ================== Now we get and save the FR blood ============
        timesBlood = Blood.values(:, strcmp(Blood.labels, 'Time'));
        PF = Blood.values(:, strcmp(Blood.labels, 'Parent Plasma Radioactivity'));
        WB = Blood.values(:, strcmp(Blood.labels, 'Whole Blood Radioactivity'));
        TP = Blood.values(:, strcmp(Blood.labels, 'Total Plasma Radioactivity'));

        % resolution
        dlmwrite(fullfile(FR_metabolite_corrected_signal_data, [filename, '.txt']), PF);
        dlmwrite(fullfile(FR_whole_blood_signal_data, [filename, '.txt']), WB);
        dlmwrite(fullfile(FR_signal_data, [filename, '.txt']), TP);

        % ================== Now we get and save the LR blood ============
        PF_short = downsampleBlood(FrameStartEnd, timesBlood, PF);
        WB_short = downsampleBlood(FrameStartEnd, timesBlood, WB);
        TP_short = downsampleBlood(FrameStartEnd, timesBlood, TP);
        dlmwrite(fullfile(metabolite_corrected_signal_data, [filename, '.txt']), PF_short);
        dlmwrite(fullfile(whole_blood_signal_data, [filename, '.txt']), WB_short);
        dlmwrite(fullfile(signal_data, [filename, '.txt']), TP_short);
        dlmwrite(fullfile(blood_time, [filename, '.txt']), timesBlood);

        % just storing blood for plot/check
        LR_blood = cat(2, LR_blood, PF_short); 
        FR_blood = [FR_blood, PF]; 
        FR_blood_SUV = [FR_blood_SUV, PF.*CorrectiveFactor(idIdx)]; 


        % Now split 4D PET in 3DxFrames .. 
        PET_4D_img = niftiread(fullfile(DynamicPET, [PrjAI_ID{idIdx}, '.nii.gz']));
        subject_folder = fullfile(image_data, [num2str(idIdx - 1)]);
        mkdir(subject_folder);

        for frame_number = 1:size(PET_4D_img, 4)
            frame_image = PET_4D_img(:, :, :, frame_number);

            % Generate the filename with the format idIdx_00framenumber.nii.gz
            filenameIMG = sprintf('%s_%04d.nii.gz', num2str(idIdx - 1), frame_number - 1);

            % Save the 3D image to the subject's folder
            niftiwrite(frame_image, fullfile(subject_folder, filenameIMG));
        end

        % ================== Now we save the clinical features ============ 
        file_path = fullfile(clinical_features_data, [filename, '.txt']);
        % Open the file for writing
        fid = fopen(file_path, 'w');
        % Check if the file is successfully opened
        if fid == -1
            error('Error opening the file for writing.');
        end
        % Write the lines to the file
        fprintf(fid, 'affinity(HAB=1):=%s\n', num2str(affinity_bin(idIdx)));
        fprintf(fid, 'sex(M=1):=%s\n', num2str(sex_bin(idIdx)));
        fprintf(fid, 'age(y):=%s\n', num2str(Age(idIdx)));
        fprintf(fid, 'SUBJUniqueIdentifier:=%d\n', SUBJUniqueIdentifier((idIdx)));
        fprintf(fid, 'Health Status(HC=0):=%d\n', (PtStatusB(idIdx)));
        fprintf(fid, 'Visit Number(baseline=0):=%d\n', (ScanNumber(idIdx)));
        fprintf(fid, 'Weight[kg]:=%.2f\n',(Weight(idIdx)));
        fprintf(fid, 'Dose[MBq]:=%.2f\n',(Dose(idIdx)));
        fprintf(fid, 'Tracer:=[11C]PBR28\n');
        fprintf(fid, 'Scanner:=KCL PET/MR\n');

        % Close the file
        fclose(fid);

        % ================== Now we save the TACs ============
        load(fullfile(TACsDir, matFileName), 'TACs');
        TAClabels=TACs.labels;
        % Initialize an empty array to store modified labels
        modifiedLabels = cell(size(TAClabels));
        
        % Loop through each file
        for TACindex = 1:numel(TAClabels)
            % Extract label from the file name (remove extension)
            [~, fname, ~] = fileparts(TAClabels(TACindex));
            
            % Check if the label contains '_R' or '_L'
            if contains(fname, '_R') || contains(fname, '_L')
                % Assign NaN if it contains '_R' or '_L'
                modifiedLabels{TACindex} = 'NO';
            else
                % Keep the label unchanged
                modifiedLabels{TACindex} = fname;
            end
        end
        subset=~strcmp(modifiedLabels, 'NO');
        TAC_subset = TACs.activity(:,subset); 
        subject_folderTACs = fullfile(time_activity_curves, [num2str(idIdx - 1)]);
        mkdir(subject_folderTACs);
        % Loop through each file
        for indexTACs = 1:size(TAC_subset,2)
            filenameTACs = sprintf('%s_%d', num2str(idIdx - 1), indexTACs);
            dlmwrite(fullfile(subject_folderTACs, [filenameTACs, '.txt']), TAC_subset(:,indexTACs));
        end
        
        % ================== Just some checks at the end ============
        if idIdx==1
            midFrames_template=midFrames;
            subset_template=subset;
         end

         % Check if midFrames is equal to midFrames_template
        if ~isequal(midFrames, midFrames_template)
            error('Error: midFrames is not equal to midFrames_template. Check the data.');
        else
            fprintf('difference in time is: %f\n',sum(midFrames_template-midFrames));
        end

         % Check if midFrames is equal to midFrames_template
        if ~isequal(subset, subset_template)
            error('Error: subset is not equal to subset_template. Check the data.');
        else
            fprintf('difference in subset is: %f\n',sum(midFrames_template-midFrames));
        end

        fprintf('All data file created successfully for subjID %d \n',idIdx-1);
    end
end


% Assuming you have already defined midFrames, LR_blood, timesBlood, and FR_blood

%%
figure(2);

% Subplot 1
subplot(2, 2, 1);
plot(midFrames(1:15), LR_blood(1:15, :));
title('LR Blood Data');
xlabel('Time');
ylabel('LR Blood');
set(gca,'FontSize',15)
grid on;
ylim([0,200])
% Subplot 2
subplot(2, 2, 2);
for plotindex = 1:70
    plot(timesBlood(1:800), FR_blood{plotindex}(1:800));
    hold on;
end
hold off;
title('FR Blood Data');
xlabel('Time');
ylabel('FR Blood');
grid on;
set(gca,'FontSize',15)
ylim([0,200])


% Subplot 2
subplot(2, 2, 3);
plot(midFrames(1:15), TAC_subset(1:15, :));
title('TACs (one subj)');
xlabel('Time');
ylabel('TACs');
set(gca,'FontSize',15)
grid on;

% Set white background for both subplots
set(gcf, 'Color', 'w');

% Subplot 2
subplot(2, 2, 4);
for plotindex = 1:70
    plot(timesBlood(1:800), FR_blood_SUV{plotindex}(1:800));
    hold on;
end
hold off;
title('FR Blood Data SUV');
xlabel('Time');
ylabel('FR Blood');
grid on;
set(gca,'FontSize',15)

%% ALL FUNCTIONS HERE

% createDataStructure
function createDataStructure(outfolder)
    % Create the directory if it doesn't exist
    if exist(outfolder, 'dir') ~= 7
        mkdir(outfolder);
        disp(['Directory ', outfolder, ' created successfully.']);
    else
        disp(['Directory ', outfolder, ' already exists.']);
    end

    % Define subdirectories
    subdirs = {'clinical_features_data', 'image_data', 'time_activity_curves', 'image_derived_input_function_signal_data', 'metabolite_corrected_signal_data', 'whole_blood_signal_data', 'signal_data','FR_metabolite_corrected_signal_data', 'FR_signal_data', 'FR_whole_blood_signal_data','PETmidframe','PETframelength','blood_time','PETframestartstop'};

    % Create subdirectories and assign to variables
    for index = 1:length(subdirs)
        subdir = fullfile(outfolder, subdirs{index});
        
        % Create the subdirectory if it doesn't exist
        if exist(subdir, 'dir') ~= 7
            mkdir(subdir);
            disp(['Directory ', subdir, ' created successfully.']);
        else
            disp(['Directory ', subdir, ' already exists.']);
        end
        
    end
end

% exportDirectoriesToWorkspace
function exportDirectoriesToWorkspace(directoryPath)
    % Get a list of directories inside the specified path
    dirList = dir(directoryPath);
    
    % Filter out non-directories and remove '.' and '..'
    dirList = dirList([dirList.isdir] & ~ismember({dirList.name}, {'.', '..'}));

    % Create variables for each directory in the workspace
    for indexlist = 1:length(dirList)
        dirName = dirList(indexlist).name;
        assignin('base', dirName, fullfile(directoryPath, dirName));
    end
end


function plas_short = downsampleBlood(time, plastime, blood)
    plastime=plastime';
    plasmet=cat(2,linspace(1,size(blood, 1), size(blood,1))',blood); 
    tstart=time(1,1)+1; pstart=tstart+1; 
    plas_short=zeros(size(time(:,1)));

    for m=1:length(time(:,1))
        tend=time(m,2);
        if tend >= plastime(end)
            tend=plastime(end)-1;
        end
    
        for x=1:length(plasmet(:,1))
            if eq(plasmet(x,1),tend)
                pend=plasmet(x,1)+1;
            end
        end
    
        plas_short(m,1)=mean(plasmet(pstart:pend,2));
        pstart=pend+1;
    end

end