function Visualize_Dual_Color_Data_Zika2017(varargin)

    %Debugging
%     dbstop in Visualize_Find_Data at 66
    set(0, 'DefaultAxesFontSize',15)
    
    %First, we load the .mat data files.
    if length(varargin) == 1
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
            varargin{1},'Multiselect', 'on');
    elseif length(varargin) == 2
        DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
    else
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
    end
    
%     Calculate_Bleedthrough(DataFilenames,DefaultPathname)
    Calculate_Number_Single_Molecules(DataFilenames,DefaultPathname)
end

function Calculate_Number_Single_Molecules(DataFilenames,DefaultPathname)
    IntensityPerMolecule = 595;
        % 804 = 500 ms (from approximately 30 points)
    BleedthroughPercentage = 0;
        % TR into Alexa 488

    Color1List  = [];
    Color2RoughList =[];
    Color2GaussianList =[];
    NumberMoleculesList = [];
    
    if iscell(DataFilenames) 
        NumberFiles = length(DataFilenames);
    else
        NumberFiles = 1;
    end
    
    for CurrentFileNumber = 1:NumberFiles
        if iscell(DataFilenames) 
            CurrDataFileName = DataFilenames{1,CurrentFileNumber};
        else
            CurrDataFileName = DataFilenames;
        end

        CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);

        InputData = open(CurrDataFilePath);
        BindingDataToSave = InputData.BindingDataToSave;
        FileNumbers(CurrentFileNumber) = CurrentFileNumber;

        for b = 1:length(BindingDataToSave)
            VirusData = BindingDataToSave(b).VirusData;

            for j = 1:length(VirusData)
                CurrentVirusData = VirusData(j);
                if strcmp(CurrentVirusData.IsVirusGood,'y') &&...
                        CurrentVirusData.GaussianIntensity2 > 0
                    
                   Color1Intensity = CurrentVirusData.IntensityBackSub;
                   Color2IntensityRough = CurrentVirusData.RoughIntensity2;
                   Color2IntensityGauss = CurrentVirusData.GaussianIntensity2;
                   NumberMolecules = round((Color2IntensityGauss - ...
                       BleedthroughPercentage*Color1Intensity)/IntensityPerMolecule);

                   Color1List = [Color1List Color1Intensity];
                   Color2RoughList = [Color2RoughList Color2IntensityRough];
                   Color2GaussianList = [Color2GaussianList Color2IntensityGauss];
                   NumberMoleculesList = [NumberMoleculesList NumberMolecules];
                end
            end
        end 
    end
    
    FigureHandles.NumberMoleculesWindow = figure(1);
    set(FigureHandles.NumberMoleculesWindow,'Position',[1   479   451   300]);
    cla
    FigureHandles.Color2Window= figure(5);
    set(FigureHandles.Color2Window,'Position',[858   405   450   300]);
    cla
    FigureHandles.NumberNotZeroWindow = figure(2);
    set(FigureHandles.NumberNotZeroWindow ,'Position',[452   476   450   300]);
    cla
    FigureHandles.Color1Window = figure(3);
    set(FigureHandles.Color1Window,'Position',[452 50 450 300]);
    cla
    FigureHandles.DiagnosticWindow= figure(4);
    set(FigureHandles.DiagnosticWindow,'Position',[1 50 450 300]);
    cla
    
        
    set(0,'CurrentFigure',FigureHandles.Color1Window);
        NotZeroIndex = NumberMoleculesList > 0;
        XData = Color1List(NotZeroIndex);
        YData = Color2GaussianList(NotZeroIndex);
        CorrelationCoefficient = corrcoef(XData,YData);
        FigureTitle = strcat('TR VS 488 (DNA Only), r =', num2str(CorrelationCoefficient(2,1)));
        plot(XData,YData,'bo')
        title(FigureTitle)
        xlabel('Texas red Intensity');
        ylabel('488 Intensity');
        
%         DatatoShow = Color1List;
%         FigureTitle = strcat('Oregon Green Intensity, mean =',num2str(mean(DatatoShow)),...
%             '; n =',num2str(length(DatatoShow)));
%         hist(DatatoShow)
%         title(FigureTitle)
%         xlabel('Intensity');
%         ylabel('Number of Particles');
    
    set(0,'CurrentFigure',FigureHandles.Color2Window);
        NotZeroIndex = NumberMoleculesList > 0;
        XData = Color1List(NotZeroIndex);
        YData = Color2GaussianList(NotZeroIndex);
        DatatoShow = YData./XData;
        FigureTitle = strcat('488/TR (DNA Only), med =',num2str(median(DatatoShow)),...
            '; mean =',num2str(mean(DatatoShow)));
        hist(DatatoShow)
        title(FigureTitle)
        xlabel('Ratio');
        ylabel('Number of Particles');
%         
%         DatatoShow = Color2GaussianList;
%         FigureTitle = strcat('546 Gaussian Intensity, mean =',num2str(mean(DatatoShow)));
%         hist(DatatoShow)
%         title(FigureTitle)
%         xlabel('Intensity');
%         ylabel('Number of Particles');
        
    set(0,'CurrentFigure',FigureHandles.DiagnosticWindow);
        XData = Color1List;
        YData = Color2GaussianList;
        CorrelationCoefficient = corrcoef(XData,YData);
        FigureTitle = strcat('TR VS 488 (All), r =', num2str(CorrelationCoefficient(2,1)));
        plot(XData,YData,'bo')
        title(FigureTitle)
        xlabel('Texas red Intensity');
        ylabel('488 Intensity');
        
%         DatatoShow = Color2RoughList;
%         FigureTitle = strcat('546 Rough Intensity, mean =',num2str(mean(DatatoShow)));
%         hist(DatatoShow)
%         title(FigureTitle)
%         xlabel('Intensity');
%         ylabel('Number of Particles');

    set(0,'CurrentFigure',FigureHandles.NumberMoleculesWindow);
        DatatoShow = NumberMoleculesList;
        FractionWithDNA = length(NumberMoleculesList(NumberMoleculesList >0))/length(NumberMoleculesList);
        FigureTitle = strcat('All Viruses, med =',num2str(median(DatatoShow)),...
            '; Fract =', num2str(FractionWithDNA),...
            '; n =',num2str(length(DatatoShow)));
%         Edges = 0:10:max(DatatoShow)+10;
        Edges = 0:1:30;
        hist(DatatoShow,Edges)
        title(FigureTitle)
        xlabel('Number Of DNA-Lipids');
        ylabel('Number of Particles');
        xlim([0 max(Edges)]);
    
        
    set(0,'CurrentFigure',FigureHandles.NumberNotZeroWindow);
        NotZeroIndex = NumberMoleculesList > 0;
        DatatoShow = NumberMoleculesList(NotZeroIndex);
        FigureTitle = strcat('Viruses With DNA, med =',num2str(median(DatatoShow)),...
            '; mean =',num2str(mean(DatatoShow)), '; n =',num2str(length(DatatoShow)));
%         Edges = 0:5:max(DatatoShow)+5;
        Edges = 0:1:30;
        hist(DatatoShow,Edges)
        xlim([0 max(Edges)]);
        title(FigureTitle)
        xlabel('Number Of DNA-Lipids');
        ylabel('Number of Particles');
        
    disp('Thank you.  Come Again.')
    Stop_ThisWillCauseError
end

function Calculate_Bleedthrough(DataFilenames,DefaultPathname)
    Color1List  = [];
    Color2RoughList =[];
    Color2GaussianList =[];
    BleedthroughList = [];
    
    if iscell(DataFilenames) 
        NumberFiles = length(DataFilenames);
    else
        NumberFiles = 1;
    end
    
    for CurrentFileNumber = 1:NumberFiles
        if iscell(DataFilenames) 
            CurrDataFileName = DataFilenames{1,CurrentFileNumber};
        else
            CurrDataFileName = DataFilenames;
        end

        CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);

        InputData = open(CurrDataFilePath);
        BindingDataToSave = InputData.BindingDataToSave;
        FileNumbers(CurrentFileNumber) = CurrentFileNumber;

        for b = 1:length(BindingDataToSave)
            VirusData = BindingDataToSave(b).VirusData;

            for j = 1:length(VirusData)
                CurrentVirusData = VirusData(j);
                if strcmp(CurrentVirusData.IsVirusGood,'y') &&...
                        CurrentVirusData.GaussianIntensity2 > 0
                    
                   Color1Intensity = CurrentVirusData.IntensityBackSub;
                   Color2IntensityRough = CurrentVirusData.RoughIntensity2;
                   Color2IntensityGauss = CurrentVirusData.GaussianIntensity2;
                   BleedthroughRatio = Color2IntensityGauss/Color1Intensity;
%                    BleedthroughRatio = Color2IntensityRough/Color1Intensity;

                   Color1List = [Color1List Color1Intensity];
                   Color2RoughList = [Color2RoughList Color2IntensityRough];
                   Color2GaussianList = [Color2GaussianList Color2IntensityGauss];
                   BleedthroughList = [BleedthroughList BleedthroughRatio];
                end
            end
        end 
    end
    
    FigureHandles.BleedthroughWindow = figure(3);
    set(FigureHandles.BleedthroughWindow,'Position',[1 50 450 300]);
    cla
    FigureHandles.DiagnosticWindow= figure(4);
    set(FigureHandles.DiagnosticWindow,'Position',[452 50 450 300]);
    cla
    FigureHandles.Color1Window = figure(1);
    set(FigureHandles.Color1Window,'Position',[6   479   451   300]);
    cla
    FigureHandles.Color2Window = figure(2);
    set(FigureHandles.Color2Window,'Position',[472   476   450   300]);
    cla
        
    set(0,'CurrentFigure',FigureHandles.Color1Window);
        DatatoShow = Color1List;
        FigureTitle = strcat('Oregon Green Intensity, mean =',num2str(mean(DatatoShow)),...
            '; n =',num2str(length(DatatoShow)));
        hist(DatatoShow)
        title(FigureTitle)
        xlabel('Intensity');
        ylabel('Number of Particles');
    
    set(0,'CurrentFigure',FigureHandles.Color2Window);
        DatatoShow = Color2RoughList;
        FigureTitle = strcat('546 Rough Intensity, mean =',num2str(mean(DatatoShow)));
        hist(DatatoShow)
        title(FigureTitle)
        xlabel('Intensity');
        ylabel('Number of Particles');
        
    set(0,'CurrentFigure',FigureHandles.DiagnosticWindow);
        DatatoShow = Color2GaussianList;
        FigureTitle = strcat('546 Gaussian Intensity, mean =',num2str(mean(DatatoShow)));
        hist(DatatoShow)
        title(FigureTitle)
        xlabel('Intensity');
        ylabel('Number of Particles');
        
    set(0,'CurrentFigure',FigureHandles.BleedthroughWindow);
        DatatoShow = BleedthroughList;
        FigureTitle = strcat('Bleedthrough Ratio, mean =',num2str(mean(DatatoShow)));
        hist(DatatoShow)
        title(FigureTitle)
        xlabel('Ratio');
        ylabel('Number of Particles');
    
    disp('Thank you.  Come Again.')
    Stop_ThisWillCauseError
end