function [] = Start_Dual_Color_Analysis(varargin)
% - - - - - - - - - - - - - - - - - - - - -

% Input:
% Start_Dual_Color_Analysis(), in this case the user navigates to the image 
%       stacks and also chooses the parent folder where the output analysis files will be saved
%   OR
% Start_Dual_Color_Analysis(DefaultPath), where DefaultPath is the directory to which 
%       the user will be automatically directed to find the image 
%       stacks. After choosing the stacks, the user then chooses the 
%       parent folder where the output analysis files will be saved.
%   OR
% Start_Dual_Color_Analysis(DefaultPath,SavePath), where DefaultPath is as above, and 
%       SavePath is the parent folder where the output analysis files will be saved

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. The information about the number of viruses bound, together 
% with the intensity in each color image, will be in the BindingDataToSave 
% structure, as defined in Find_And_Process_Virus.m.

% Note: This program has been designed to process many sets of images sequentially,
% but it has been tested with individual sets, so keep that 
% in mind if you choose to process many sets at once.

% Original version of script by Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048

% This is a slightly modified version of the original code to analyze Sendai
% virus binding data. The main differences are different parameters and options used in
% the analysis that are optimized for the Sendai virus data.

% Modifications by:
% Bob Rawle and Amy Lam, Williams College, 2021
% Published online in conjunction with:
% Lam et al, 2021, Single virus assay reveals membrane determinants and 
% mechanistic features of Sendai virus binding

% - - - - - - - - - - - - - - - - - - - - -
%Define which options will be used
[Options] = Setup_Options();
Threshold = Options.Threshold;
close all


    %First, we load the .tif files.  Should be an image stack.  We'll also
    %set up the save folder.
    if length(varargin) == 1
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed',...
            varargin{1},'Multiselect', 'on');
        SavePath = uigetdir(varargin{1},'Choose the directory where data folder will be saved');
    elseif length(varargin) == 2
        SavePath = varargin{1,2};
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed',...
            varargin{1,1},'Multiselect', 'on');
    else
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed', 'Multiselect', 'on');
        SavePath = uigetdir(DefaultPath,'Choose the directory where data folder will be saved');
    end

    % Create the save folder
        if  strcmp(Options.DualColor,'y')
             if  strcmp(Options.ReverseFind,'y')
                 DataFolderName = strcat('/Binding Analysis','/');
                 %DataFolderName = strcat('/Reverse Find Analysis','/');
             else 
                 DataFolderName = strcat('/Binding Analysis','/');
                 %DataFolderName = strcat('/Dual Color Analysis','/');
             end
        else
            DataFolderName = strcat('/Binding Analysis','/');
        end
        
        SaveDataPathname = strcat(SavePath,DataFolderName);
        mkdir(SaveDataPathname);

if iscell(StackFilenames) %This lets us know if there is more than one file
    NumberOfFiles = length(StackFilenames);
else
	NumberOfFiles = 1;
end

for i = 1:NumberOfFiles
    if iscell(StackFilenames) 
        CurrentFilename = StackFilenames{1,i};
    else
        CurrentFilename = StackFilenames;
    end
    CurrStackFilePath = strcat(DefaultPath,CurrentFilename);
    CharFileName = char(CurrentFilename);

    % Now we call the function to find the virus particles and extract
    % their fluorescence intensity
    [BindingDataToSave, OtherDataToSave] = ...
        Find_And_Process_Virus(CurrStackFilePath,Threshold,CurrentFilename,...
            i, DefaultPath,Options);

    save(strcat(SaveDataPathname,StackFilenames(1:end),Options.DataFileLabel,'-AnalysisFile','.mat'));
end

close all
disp('---------------------')
disp ('Thank you.  Come again.')

end