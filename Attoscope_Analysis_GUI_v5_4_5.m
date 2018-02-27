function Attoscope_Analysis_GUI
%///////////////ATTOSCOPE BEAM IMAGE ANALYSIS ////////////////////////
%/        by Maria Weikum (maria.weikum@desy.de)                    /
%/                  Version 5.4.5, 2018                               /
%/                Last edit: 09/02/2018                             /

%// Attoscope beam image analysis: converts attoscope beam image on screen
% into an image of the longitudinal bunch profile and outputs the RMS
% longitudinal beam size as well as potential microbunching//
%// for details, see .pdf documentation file

%// changes from version 5.4.4:
%   - use of a manually decided point for sine fit amplitude-starting-point
%  if "ideal amplitude" checkbox is not turned on
%   - deleted offset along x as a variable parameter for the sine fit; it
%   is now by default 0
%   - instead of flipping y coordinates (which causes reconstruction
%   errors), the image data array structure is flipped with flipud
%   - improved visual layout of GUI


%------------------GUI INITIATION--------------------------------
%----------------------------------------------------------------
%create and hide GUI
f=figure('Visible','off','Position',[160,300,1000,800]);
set(f,'Units','normalized');

%construct the figures within the GUI
ha1=axes('Units','Pixels','Position',[55,445,425,325]);
ha2=axes('Units','Pixels','Position',[530,445,425,325]);
set(ha1,'Units','normalized');
set(ha2,'Units','normalized');
set(f,'toolbar','figure');

%construct the text fields within the GUI
% nl,nr = no. of buttons per column 
% x0l,x0r,y0l,y0r = x- and y-positions of top buttons in left and right row
% dx0,dy0 = x- and y-dimensions of the buttons
% gapy = vertical gap between buttons
htext.name={'Drift length L1: undulator - TDS [m]:';'Drift length L2: TDS - screen [m]:';...
    'Electron energy [MeV]:';'Laser power [W]:';'Laser wavelength [m]:'; 'Laser spot size [m]:';...
    'Laser intensity RMS length [s]:';'Undulator strength K:';'No. of undulator periods:';...
    'Undulator wavelength [m]:';'Deflector strength A_rf = V_rf[V] / E[eV]:';...
    'Deflector RF wavelength [m]:';'Profile resolution [m]'};
nl=ceil(size(htext.name,1)/2); nr=size(htext.name,1)-nl; 
x0l=55; y0l=355; x0r=280; y0r=355; dx0=135; dy0=45; gapy=8;
for i=1:nl %left column of text fields
    htext.pos(i,1)={[x0l,(y0l-(i-1)*dy0-(i-1)*gapy),dx0,dy0]};
    htext.button(i,1)={['htext' num2str(i) 'l']};
end
for i=1:nr %right column of text fields
    htext.pos(nl+i,1)={[x0r,(y0r-(i-1)*dy0-(i-1)*gapy),dx0,dy0]};
    htext.button(nl+i,1)={['htext' num2str(i) 'r']};
end
for i=1:size(htext.button,1) %initialise button and set for automatic resizing
    htext.button{i,1}=uicontrol(f,'Style','text','String',htext.name{i,1},'Position',htext.pos{i,1});
    set(htext.button{i,1},'Units','normalized'); 
    align([htext.button{i,1}],'Center','None');
end

%construct the text input fields within the GUI and define variables and initial values
%for each field
%same position and size definitions as above in htext
hedit.initval={'0.5';'0.5';'44';'100e9';'10.6e-6';'1.1e-3';'3e-12';'2.5';'10';'4e-2';'0.2273';'0.026';'1e-7'};
hedit.var={'L1';'L2';'E';'PL';'lambda';'wL';'tauL';'K';'Nu';'lambda_u';'Arf';'lambda_rf';'dsin'};
for i=1:size(hedit.var,1)
    s.(hedit.var{i})=str2double(hedit.initval{i});
end
x0led=195; y0led=355; x0red=420; y0red=355; dx0ed=60; dy0ed=45; gapyed=8; %initial values
for i=1:nl %left column of input fields
    hedit.button(i,1)={['hedit' num2str(i) 'l']};
    hedit.pos(i,1)={[x0led,(y0led-(i-1)*dy0ed-(i-1)*gapyed),dx0ed,dy0ed]};
end
for i=1:nr %right column of input fields
    hedit.button(nl+i,1)={['hedit' num2str(i) 'l']};
    hedit.pos(nl+i,1)={[x0red,(y0red-(i-1)*dy0ed-(i-1)*gapyed),dx0ed,dy0ed]};
end
for i=1:size(hedit.button,1) %initialise button and set for automatic resizing
    hedit.button{i,1}=uicontrol(f,'Style','edit','String',hedit.initval{i,1},'Position',...
        hedit.pos{i,1},'Callback',{@hedit_Callback});
    set(hedit.button{i,1},'Units','normalized');
    align([hedit.button{i,1}],'Center','None');
    hedit.initval{i}=str2double(hedit.initval{i});
end

%construct and initialise remaining push buttons in GUI
%h_loadbutton: loads input image
h_loadbutton=uicontrol(f,'Style','pushbutton','String','Load data','Position',...
    [530,350,187.5,50],'Callback',{@h_loadbutton_Callback});
%h_bgbutton: loads input background image
h_bgbutton=uicontrol(f,'Style','pushbutton','String','Remove background','Position',...
    [530,280,187.5,50],'Callback',{@h_bgbutton_Callback});
%h_longbutton: reconstructs long. beam profile from input image (with
%background subtracted)
h_longbutton=uicontrol(f,'Style','pushbutton','String','Long. profile - single shot',...
    'Position',[530,210,187.5,50],'Callback',{@h_longbutton_Callback});
%h_longbuttonm: reconstructs long. beam profile from multiple single reconstructed
%profiles to be loaded in
h_longbuttonm=uicontrol(f,'Style','pushbutton','String','Long. profile - multi shot',...
    'Position',[530,140,187.5,50],'Callback',{@h_longbuttonm_Callback});
%h_savebutton: saves the reconstructed profile and some data about it
h_savebutton=uicontrol(f,'Style','pushbutton','String','Save data','Position',...
    [530,70,187.5,50],'Callback',{@h_savebutton_Callback});
%h_measbutton: measures the distance between two specific points along the reconstructed
%profile
h_measbutton=uicontrol(f,'Style','pushbutton','String','Measure distance','Position',...
    [767.5,210,187.5,50],'Callback',{@h_measbutton_Callback});
%h_microbunchbutton: detects microbunches within the reconstructed profile and calculates their
%width and distance
h_microbunchbutton=uicontrol(f,'Style','pushbutton','String','Measure microbunching',...
    'Position',[767.5,140,187.5,50],'Callback',{@h_microbunchbutton_Callback});
%h_helpbutton: opens up a help document with more information
h_helpbutton=uicontrol(f,'Style','pushbutton','String','HELP','Position',...
    [802.5,10,187.5,30],'Callback',{@h_helpbutton_Callback});

%normalise the units for positioning of the buttons
set(h_loadbutton,'Units','normalized');
set(h_bgbutton,'Units','normalized');
set(h_longbutton,'Units','normalized');
set(h_longbuttonm,'Units','normalized');
set(h_measbutton,'Units','normalized');
set(h_savebutton,'Units','normalized');
set(h_microbunchbutton,'Units','normalized');
set(h_helpbutton,'Units','normalized');

%construct and initialise check boxes in GUI
%h_deflcheck: check if deflector is turned on
h_deflcheck=uicontrol(f,'Style','checkbox','String','Deflector is turned on.',...
    'Position',[767.5,340,187.5,50],'Value',1,'Callback',{@h_deflcheck_Callback});
%h_undcheck: check if undulator is turned on
h_undcheck=uicontrol(f,'Style','checkbox','String','Undulator is turned on.',...
    'Position',[767.5,280,187.5,50],'Value',1,'Callback',{@h_undcheck_Callback});
%h_debugcheck: check for activating debugging features, e.g. additional plots
h_debugcheck=uicontrol(f,'Style','checkbox','String','Debug mode','Position',...
    [680,10,127.5,30],'Value',0,'Callback',{@h_debugcheck_Callback});
%h_idealparcheck: check if ideal value of laser wavelength should be used
h_idealparcheck=uicontrol(f,'Style','checkbox','String','wavelength','Position',...
    [450,35,95,30],'Callback',{@h_idealparcheck_Callback});
%h_idealparcheck2: check if ideal value of amplitude of sine-curve structure 
%in image should be used
h_idealparcheck2=uicontrol(f,'Style','checkbox','String','amplitude ','Position',...
    [545,35,90,30],'Callback',{@h_idealparcheck2_Callback});
%text to go with idealparcheck checkboxes
idealparbutton=uicontrol(f,'Style','text','String','Use ideal sine parameters:',...
    'Position',[270,28,180,30]);
%h_manualparcheck2: check if manually defined value of amplitude of sine-curve 
%structure in image should be used; this is either found from the mean
%position of high intensity data points in the signal or, for low-intensity
%signals, determined by the user through manually marking a position
h_manualparcheck2=uicontrol(f,'Style','checkbox','String','amplitude ','Position',...
    [545,5,90,30],'Callback',{@h_manualparcheck2_Callback});
%text to go with manualparcheck checkboxes
manualparbutton=uicontrol(f,'Style','text','String','Use manual sine parameters:',...
    'Position',[270,0,200,30]);

%normalise the units for positioning of the checkboxes
set(h_deflcheck,'Units','normalized');
set(h_undcheck,'Units','normalized');
set(h_debugcheck,'Units','normalized');
set(h_idealparcheck,'Units','normalized');
set(h_idealparcheck2,'Units','normalized');
set(idealparbutton,'Units','normalized'); 
set(h_manualparcheck2,'Units','normalized');
set(manualparbutton,'Units','normalized'); 

%define constants
import constants;
IA=4*pi*constants.eps0*constants.me*constants.c^3/constants.e; %Alfven current in Ampere
P0=IA*constants.me*constants.c^2/constants.e;

%define variables to be read in all nested functions
filename_wd='-'; filename_bg='-'; 
I=[]; X=[]; Y=[];
und_check=1; defl_check=1; msg_def=''; msg_und=''; 
debug_check=0; idealpar_check=0; idealpar_check2=0; manualpar_check2=0;
x_amp=0; refp=0;
 Arf=0; krf=0; A=0; k=0; L=0;
s0_final=[]; I_lin=[]; s0_final_lin=[]; I_lin_lin=[];
s0mean=0; s0RMS=0; meanFWHM=0; meanFWHMhole=0; ulimit=0; llimit=0;
medianFWHM=0; medianFWHMhole=0; stdFWHM=0; stdFWHMhole=0;
multi=[]; MultiFilename=[];

%calculates Arf,krf,A,k from initial input values
parameters_calc() 

%initialize GUI
set(f,'Name','Attoscope Analysis'); %title
movegui(f,'center') %recenter window
set(f,'Visible','on'); %turn visible


%-------------------CALLBACK FOR LOAD BUTTON:h_loadbutton -----------------
%--------------------------------------------------------------------------
%load a data file and plot it; this can be either .tif or .asc format
function h_loadbutton_Callback(source,eventdata)
    %...............reset microbunching and background values.............
    meanFWHM=0; meanFWHMhole=0;
    
    %...........read in input data.......................................
    filename_wd=uigetfile('.asc'); %opens window to choose input file
    if strfind(filename_wd,'.tif')~=0
        I=double(imread(filename_wd));
    elseif strfind(filename_wd,'.asc')~=0
        I=double(dlmread(filename_wd));
    end
    
    %...define image screen size and create x,y coordinates for the grid...
    screensize=inputdlg({'Enter the image screen size: x-direction [m]',...
          'Enter the image screen size: y-direction [m]'},'Screen Size',1,{'1.0','1.0'});
    x=linspace(-str2double(screensize{1,1})/2,str2double(screensize{1,1})/2,size(I,2));
    y=linspace(-str2double(screensize{2,1})/2,str2double(screensize{2,1})/2,size(I,1));
    %option to flip image data vertically; this may need to be applied,
    %depending on how the image is defined
    %y=fliplr(y);
    %I=flipud(I);
    
    [X,Y]=meshgrid(x,y);
    
    %................plot the input data..................................
    axes(ha1)
    surf(X,Y,I,'LineStyle','None')
    colorbar; view([0 90]); axis tight;
    xlabel('X [m]'); ylabel('Y [m]');
    title('Image at detector');
        
end


%------------------CALLBACK FOR BACKGROUND BUTTON:h_bgbutton-------------
%------------------------------------------------------------------------
%load a background data file to be subtracted from the input data
function h_bgbutton_Callback(source,eventdata)
    %..........read in input background data............................
    filename_bg=uigetfile('.asc'); %opens window to choose input file
     if strfind(filename_wd,'.tif')~=0
        I_bg=double(imread(filename_bg));
    elseif strfind(filename_wd,'.asc')~=0
        I_bg=double(dlmread(filename_bg));
    end
    
    %...............subtract background from signal image.............
    I=I-(I_bg); %subtract background
    I(I<0)=0; %set minimum signal level to zero
    %.......................replot updated input data..................
    axes(ha1)
    surf(X,Y,I,'LineStyle','None')
    colorbar; view([0 90]); axis tight;
    xlabel('X [m]'); ylabel('Y [m]');
    title('Image at detector');
end
    

%------------------CALLBACK FOR HEDIT FIELDS:hedit{1,i}------------------
%------------------------------------------------------------------------
%read out the input from the text input fields
function hedit_Callback(hObject,eventdata)
    ind=find([hedit.button{:,1}]==hObject);
    s.(hedit.var{ind})=str2double(get(hObject,'string'))
    parameters_calc() %recalculate Arf,krf,A,k
end


%----------------PARAMETERS_CALC()--------------------------------------
%-----------------------------------------------------------------------
%calculate theoretical parameters Arf, krf, A, k, L (see documentation for details)
function parameters_calc(source,eventdata)
    disp('- calculate theoretical parameters for undulator, laser and TDS')
    mev=constants.me*constants.c^2/(constants.e*1e6); %mc^2 in units of MeV
    Arf=s.Arf;
    %calculate krf = TDS RF wavenumber = 2pi/RF wavelength
    krf=2*pi/s.lambda_rf
    %calculate k = laser wavenumber = 2pi/laser wavelength
    k=2*pi/s.lambda
    %calculate A = laser modulator "interaction strength"=2K/gamma^2*sqrt(laser power/P0)*JJ*f
    JJ=besselj(0,(s.K^2/(4+2*s.K^2)))-besselj(1,(s.K^2/(4+2*s.K^2)))
    tRMSh=k*constants.c*s.tauL/(2*pi*s.Nu) ;
    gammar=sqrt(k*s.lambda_u/(4*pi)*(1+s.K^2/2));
    q=2*s.Nu*s.lambda_u/(k*s.wL^2) ;
    v=2*s.Nu*(s.E/mev-gammar)/gammar;
    ffun=@(z) q*cos(2*pi*v*z-2*atan(q*z))./(1+(q.*z).^2).*exp(-(z/(2*tRMSh)).^2); %assuming s=0
    fun=integral(ffun,-0.5,0.5)
    A=2*s.K/(s.E/mev)^2*sqrt(s.PL/P0)*JJ*fun
    %define total drift length L 
    L=s.L1+s.L2;
end


%--------------CALLBACK FOR H_DEFLCHECK, H_UNDCHECK-----------------------
%-------------------------------------------------------------------------
%set check boxes to appropriate value to designate undulator / deflector
%presence in experiment
function h_deflcheck_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        defl_check=1;
    else
        defl_check=0;
    end
end
function h_undcheck_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        und_check=1;
    else
        und_check=0;
    end
end

%--------------CALLBACK FOR H_DEBUGCHECK----------------------------------
%-------------------------------------------------------------------------
%update value of debug checkbox determining if additional information for
%debugging is provided for bunch profile reconstruction and microbunch 
%identification (on=additional plot, off=no extra data)
function h_debugcheck_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        debug_check=1;
    else
        debug_check=0;
    end
end

%--------------CALLBACK FOR H_IDEALPARCHECK----------------------------------
%-------------------------------------------------------------------------
%update value of ideal_par checkbox determining if the theoretical (ideal)
%parameter should be used for the wavelength when reconstructing
%the sine curve on-screen
function h_idealparcheck_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        idealpar_check=1;
    else
        idealpar_check=0;
    end
end

%--------------CALLBACK FOR H_IDEALPARCHECK2----------------------------------
%-------------------------------------------------------------------------
%update value of ideal_par2 checkbox determining if the theoretical (ideal)
%parameter should be used for the amplitude when reconstructing
%the sine curve on-screen
function h_idealparcheck2_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        idealpar_check2=1;
    else
        idealpar_check2=0;
    end
end

%--------------CALLBACK FOR H_MANUALPARCHECK2----------------------------------
%-------------------------------------------------------------------------
%update value of manual_par2 checkbox determining if the manual
%parameter should be used for the amplitude when reconstructing
%the sine curve on-screen
function h_manualparcheck2_Callback(hObject,eventdata)
    if (get(hObject,'Value') == get(hObject,'Max'))
        manualpar_check2=1;
    else
        manualpar_check2=0;
    end
end

%------------------CALLBACK FOR HELP BUTTON: h_helpbutton-----------------
%------------------------------------------------------------------------
%loads a pdf-file with further instructions and details on how the script
%works
function h_helpbutton_Callback(source,eventdata)
    helpfilename='Attoscope_Beam_Image_Analysis-v5_4_Documentation.pdf'
    if exist(helpfilename,'file')==2
        if isunix
            disp('unix')
            system(['evince ' helpfilename  ' &']);
        else
            disp('not unix')
            open(helpfilename);
        end
    else
        warndlg(['Warning! Could not find help file. Check for location of ' helpfilename ...
            ' in current directory.'],'Warning');
    end
end

%-----------CALLBACK FOR MULTIPLE LONG. PROFILE BUTTON: %h_longbuttonmbutton-------
%-------------------------------------------------------------------------
%loads in multiple reconstructed long. profiles (ideally with slightly
%different relative laser-electron beam phase) and combines them into an improved profile
function h_longbuttonm_Callback(source,eventdata)
    l=0; j=1; filename=[];
    disp('LONGIITUDINAL PROFILE RECONSTRUCTION: MULTI SHOT')
    
    %..........define which profile to load..........................
    while l<1
        % Construct a question dialog box to ask for how many profiles
        % should be loaded
        choice = questdlg('Would you like to load another profile?', 'Multi-shot reconstruction',...
    'Yes','No','Yes');
        %if the box is answered with yes, a profile file name is quested
        if strcmp(choice,'Yes')
            MultiFilename{j}=uigetfile('.txt');
            filename{j}=strrep(MultiFilename{j},'_','\_');
            j=j+1;
        else
            l=l+1;
        end
    end
    
    %....load and combine profiles to produce an improved reconstruction....
    if numel(MultiFilename)~=0
        %.....load profiles...............
        refp_m=cell([1,numel(MultiFilename)]);
        s0_final_m=cell([1,numel(MultiFilename)]);
        I_lin_m=cell([1,numel(MultiFilename)]);
        I_lin_tp=cell([1,numel(MultiFilename)]);
        for l=1:numel(MultiFilename)
            Aload=load(MultiFilename{l});
            %reference point to combine profiles, e.g. central intensity
            %peak
            refp_m{l}=Aload(1,1);
            %long. beam position
            s0_final_m{l}=Aload(1,2:size(Aload,2));
            %long. beam profiles with turning point sections set to NaN,
            %normalised
            I_lin_m{l}=Aload(4,2:size(Aload,2))/max(Aload(4,2:size(Aload,2)));
            %full long. beam profiles, normalised
            I_lin_tp{l}=Aload(2,2:size(Aload,2))/max(Aload(2,2:size(Aload,2)));
        end
        
        %if debug option is set: plot individual profiles overlapped
        if debug_check==1
            f1=figure();
            subplot(1,2,1)
            legendInfo=cell([1,numel(MultiFilename)]);
            for l=1:numel(MultiFilename)
                plot(s0_final_m{l},I_lin_tp{l});
                hold on;
                legendInfo{l}=filename{l};
            end
            xlabel('s0 [m]'); ylabel('I [arb. units]');
            title('Input profiles: before synchronisation');
            hold off;
        end

        %..move profiles along s0-direction to fit their reference points..
        dx=zeros([1,numel(MultiFilename)]);
        mins0=zeros([1,numel(MultiFilename)]);
        maxs0=zeros([1,numel(MultiFilename)]);
        for l=1:numel(MultiFilename)
            dx(l)=refp_m{1}-refp_m{l};
            s0_final_m{l}=s0_final_m{l}+dx(l);
            %find max and min values for each s0_final
            mins0(l)=min(s0_final_m{l});
            maxs0(l)=max(s0_final_m{l});
        end
        
        %if debug option is set: plot individual profiles overlapped after
        %synchronising ref. points
        if debug_check==1
            figure(f1);
            subplot(1,2,2)
            for l=1:numel(MultiFilename)
                plot(s0_final_m{l},I_lin_tp{l});
                hold on;
                legendInfo{l}=filename{l};
            end
            xlabel('s0 [m]'); ylabel('I [arb. units]');
            title('Input profiles: after synchronisation');      
            hold off;
        end

        %....cut off sides of profiles outside of s0left and s0right.......
        s0left=max(mins0);
        s0right=min(maxs0);
        s0_final_cut=cell([1,numel(MultiFilename)]);
        I_lin_cut=cell([1,numel(MultiFilename)]);
        I_lin_tp_cut=cell([1,numel(MultiFilename)]);
        for l=1:numel(MultiFilename)
            %find starting index at left side for each profile
            [~, ind_cutl]=min(abs(s0_final_m{l}-s0left))
            %find end index at right side for each profile
            if l==1
                [~, ind_cutr]=min(abs(s0_final_m{l}-s0right))
                diff_cut=ind_cutr-ind_cutl;
            else
                ind_cutr=ind_cutl+diff_cut;
            end
            %cut off profiles at left and right boundaries to have same
            %length along s0
            size(s0_final_m{l})
            s0_final_cut{l}=s0_final_m{l}(ind_cutl:ind_cutr);
            I_lin_cut{l}=I_lin_m{l}(ind_cutl:ind_cutr);
            I_lin_tp_cut{l}=I_lin_tp{l}(ind_cutl:ind_cutr);
        end
        
        %......define first loaded profile as baseline for final profile....
        s0_final_fin=s0_final_cut{1};
        I_lin_fin=I_lin_cut{1}; %profile without turning point values
        I_lin_fin_full=I_lin_tp_cut{1}; %full profile

        %.....add data points from other loaded profiles to turning point
        %regions of the baseline..................
        for j=1:numel(MultiFilename)
            %define turning point region indices
            ind_tp=find(isnan(I_lin_fin));
            %add data points to baseline turning point regions
            for l=1:numel(ind_tp)
                [~, ind_tp_n]=min(abs(s0_final_fin(ind_tp(l))-s0_final_cut{j}));
                if ~isnan(I_lin_cut{j}(ind_tp_n))
                    s0_final_fin(ind_tp(l))=s0_final_cut{j}(ind_tp_n);
                    I_lin_fin(ind_tp(l))=I_lin_cut{j}(ind_tp_n);
                    I_lin_fin_full(ind_tp(l))=I_lin_cut{j}(ind_tp_n);
                else
                    I_lin_fin_full(ind_tp(l))=I_lin_tp_cut{1}(ind_tp(l));
                end
            end  
        end 

        %....option to load in reference beam profile image to compare.....
        choice2 = questdlg('Would you like to load a reference profile to compare the reconstructed profile to?', ...
        'Multi-shot reconstruction','Yes','No','No');
        if strcmp(choice2,'Yes')
            %open a specified reference profile as a figure
            figname=uigetfile('.fig');
            h=openfig(figname,'invisible');
            %extract relevant x- and y-data from figure
            axesObjs = get(h, 'Children');
            dataObjs = get(axesObjs, 'Children');
            xdata = get(dataObjs, 'XData');
            ydata = get(dataObjs, 'YData');
            %cut off similarly to the other loaded profiles
            ind_cut_img=find(xdata<=s0left | xdata>=s0right);
            xdata(ind_cut_img)=[]; ydata(ind_cut_img)=[];
        end
        
        %if debug option is set: plot comparison between reconstr. profile
        %from single image, multiple images and possibly an original
        %reference profile
        if debug_check==1
            figure();
            plot(s0_final_cut{1},I_lin_tp_cut{1},'r');
            hold on;
            plot(s0_final_fin,I_lin_fin_full,'b',s0_final_fin,I_lin_fin,'k.')
            title('Attoscope reconstruction of a long. beam profile with multiple images at different phases')
            xlabel('s0 [m]');
            ylabel('I [arb. units]');
            if strcmp(choice2,'Yes')
                plot(xdata,ydata*max(I_lin_fin_full));
                legend('reconstr. from single image','reconstr. from multiple images',' ','original');
            else
                legend('reconstr. from single image','reconstr. from multiple images',' ');
            end
            hold off
        end

        %....plot recovered profile from multishot analysis as output....
        axes(ha2);
        plot(s0_final_fin,I_lin_fin_full);
        xlabel('S0 [m]'); ylabel('I [arb. units]');
        title('Reconstructed long. bunch profile - multishot');
        
        %........calculate the mean and standard deviation of the recovered
        %profile........................................................
        s0mean=sum(s0_final_fin.*I_lin_fin_full)./(sum(I_lin_fin_full));
        s0var=sum(I_lin_fin_full.*(s0_final_fin-s0mean).^2)./(sum(I_lin_fin_full)-1);
        s0RMS=sqrt(s0var);
        disp(['RMS length = ' num2str(s0RMS) 'm'])
        %add s0RMS value to plot
        axes(ha2);
        xcord=xlim(ha2);
        ycord=ylim(ha2);
        text(xcord(1)+0.1*abs(xcord(1)),ycord(2)-0.1*abs(ycord(2)),['RMS length = ' num2str(s0RMS) 'm'],'Parent',ha2);
        hold off

        %...............update values of global variables.............
        multi='y';
        refp=refp_m{1};
        s0_final=s0_final_fin;
        I_lin=I_lin_fin_full;  
    end
    
    clear Aload
    MultiFilename=[];

end
        
%-----------CALLBACK FOR LONG. PROFILE BUTTON:h_longbutton---------------
%------------------------------------------------------------------------
%run through algorithm for converting x-y beam image to an image of the 
%longitudinal profile of the beam 
function h_longbutton_Callback(source,eventdata)
    disp('LONGIITUDINAL PROFILE RECONSTRUCTION: SINGLE SHOT')
    tic %measure time for reconstruction algorithm to complete
    disp('- recenter & rescale Y to reverse deflector effect')
    %..........reset microbunching values............................
    meanFWHM=0; meanFWHMhole=0;
    
    %.........recenter image data in x and y direction................
    %find nonzero values of I & calculate mean of corresponding x values
    %****turned off in current version due to its unsuitability for ultrashort 
    %beam measurements, may need to be re-activated depending on sample
    %type
    %  meanX_Inonzero=sum(X(1,:).*sum(I,1))/sum(sum(I,1));
    %  meanY_Inonzero=sum(Y(:,1).*sum(I,2))/sum(sum(I,2));
    %center X and Y around zero
    %X=X-meanX_Inonzero;
    %Y=Y-meanY_Inonzero;
    
    %.........rescale Y to reverse deflector effect....................
    if defl_check==0 %if deflector was off, rescaling is skipped
        Y_new=Y;
        X_new=X;
        msg_def='Deflector is not used.\n';
             
    else %if deflector was on, Y is rescaled by Arf*krf*L
        parameters_calc() 
        
        %warning if streaking strength from input parameters is zero
        if (Arf*krf*s.L2)==0
            warndlg('Arf*krf*s.L2=0, recheck input parameters or uncheck deflector checkbox!'...
                ,'Warning');         
        end
     %option to include beam size in y without deflector which gives a
     %more accurate result as it considers effect of divergence in y
        %calculate beam RMS size in y direction of original screen image
        yinmean=sum(Y(:,1).*sum(I,2))./(sum(sum(I,2)));
        yinvar=sum(sum(I,2).*(Y(:,1)-yinmean).^2)./(sum(sum(I,2))-1);
        yinRMS=sqrt(yinvar);
        %input box for unstreaked beam size in y: can be delivered as
        %numerical input or through loading unstreaked beam image
        sigmyDcheck=questdlg('Do you have an estimate of the beam size in y without the deflector?', 'Undeflected beam size',...
            'Yes (numerical value in m)','Yes (screen image)','No','No');
        if strcmp(sigmyDcheck,'Yes (numerical value in m)')
            sigmyDinput=inputdlg('Please enter the beam size in y at the screen without deflection [m]',...
            'Beam Size in y',1);
            yunstrRMS=str2double(sigmyDinput{1});
            %rescale Arf to take into account effect of divergence on
            %streaking in y
            Arf=Arf/sqrt(1-(yunstrRMS/yinRMS)^2);
        elseif strcmp(sigmyDcheck,'Yes (screen image)')
            %choose input file and load unstreaked beam image
            filename_unstr=uigetfile('.tif'); 
            if strfind(filename_unstr,'.tif')~=0
                I_unstr=double(imread(filename_unstr));
            elseif strfind(filename_unstr,'.asc')~=0
                I_unstr=double(dlmread(filename_unstr));
            end
            screensize_unstr=inputdlg({'Enter the image screen size: x-direction [m]',...
          'Enter the image screen size: y-direction [m]'},'Screen Size',1,{num2str(2*max(X(:))),num2str(2*max(Y(:)))});
            y_unstr=linspace(-str2double(screensize_unstr{2,1})/2,str2double(screensize_unstr{2,1})/2,size(I_unstr,1));
            I_unstr=flipud(I_unstr);
            %calculate beam RMS size in y direction of unstreaked screen
            %image
            yunstrmean=sum(y_unstr'.*sum(I_unstr,2))./(sum(sum(I_unstr,2)));
            yunstrvar=sum(sum(I_unstr,2).*(y_unstr'-yunstrmean).^2)./(sum(sum(I_unstr,2))-1);
            yunstrRMS=sqrt(yunstrvar);
            %rescale Arf to take into account effect of divergence on
            %streaking in y
            Arf=Arf/sqrt(1-(yunstrRMS/yinRMS)^2);
        end

        Y_new=Y/(Arf*krf*s.L2);
        X_new=X;
        msg_def='Deflector is used.\n';
        
        %..........replot adjusted beam image..........................
        axes(ha1)
        surf(X_new,Y_new,I,'LineStyle','None')
        xlabel('X [m]'); ylabel('Y [m]');
        axis tight; view([0 90]); colorbar;
        title('Image at detector');
    end
    
    toc %measures elapsed time since starting point with command "tic"
    
    %........find reference point in y-outline of screen image...........
    disp('- find reference point in screen image')
    %screen image integrated along x
    Iliny=sum(I,2);
    %smoothed version of Iliny
    Ilinysm=smooth(Iliny,10);
    [maxIliny, imaxIliny]=max(Ilinysm);
    %reference point is automatically detected as maximum intensity
    %position along integrated screen image
    refp=Y_new(imaxIliny,1);
    
    %if debug option is set: plot integrated and smoothed screen images
    %with automatic reference point marked
    if debug_check==1
            figure;
            plot(Y_new(:,1),Iliny,'b',Y_new(:,1),Ilinysm,'r',refp,maxIliny,'ko');
            title('Position of reference point on integrated intensity profile');
            xlabel('y [m]'); ylabel('I integrated over x [a.u.]');
            legend('Integrated signal','Integrated smoothed signal','Reference point');
    end
    
    %mark reference point in screen image assuming x-position = 0
    axes(ha1)
    hold on; 
    refpx=0;
    plot3(refpx,refp,maxIliny,'or');
    %option for manual correction of reference point
    refcheck=questdlg('Is the reference point set to the maximum peak in the vertical direction?', ...
            'Single-shot reconstruction','Yes','No','No');
    if strcmp(refcheck,'No')
        dcm_obj1 = datacursormode(f); 
        datacursormode on;
        msg_x1=msgbox('Click the correct reference point, then press Return.');
        waitfor(msg_x1)
        cursor_info1=getCursorInfo(dcm_obj1);
        refp=cursor_info1.Position(2);
    end  
    hold off;
    axes(ha1)
    surf(X_new,Y_new,I,'LineStyle','None')
    xlabel('X [m]'); ylabel('Y [m]');
    axis tight; view([0 90]); colorbar;
    title('Image at detector');
    hold on; 
    plot3(refpx,refp,maxIliny,'ro');
    hold off;
    toc
    
    
    %.........calculate spread of curve along x as a measure of how.......%
    %..........smeared out the measured beam signal is.....................%
    %define suitable range in y-direction as region between +/- half a
    %laser wavelength from the central reference point
    [~, i1]=min(abs(Y_new(:,1)-(refp-s.lambda/2)));
    [~, i2]=min(abs(Y_new(:,1)-(refp+s.lambda/2)));
    i=i1:1:i2;
    Ylin=zeros([1,numel(i)]);
    Isum=zeros([1,numel(i)]);
    Icent=zeros([1,numel(i)]);
    %calculate the mean position and RMS width of the intensity signal
    %along the x-direction for each pixel row within the defined region in
    %y
    for l=1:numel(i)
        Xlin=X(i(l),:); Ylin(l)=Y(i(l),1);
        Ilin=I(i(l),:);
        s0mean(l)=sum(Xlin.*Ilin)./sum(Ilin);
        s0RMS(l)=sqrt(sum(Ilin.*(Xlin-s0mean(l)).^2)./(sum(Ilin)-1));
        [~, il1]=min(abs(Xlin-(s0mean(l)-s0RMS(l))));
        [~, il2]=min(abs(Xlin-(s0mean(l)+s0RMS(l))));
        [~, il3]=min(abs(Xlin-(s0mean(l))));
        Isum(l)=sum(Ilin(il1:il2));
        Icent(l)=Ilin(il3);
    end
    %calculate the mean and the maximum value of the signal width values
    %within the region in y
    SFhor=mean(s0RMS);
    SFhor2=max(s0RMS);
    disp('- find horizontal signal spread')
    disp(['Mean horizontal signal width (within the signal central region) = '...
        num2str(SFhor) 'm']);
    disp(['Max. horizontal signal width (within the signal central region) = '...
        num2str(SFhor2) 'm']);  
    
    
    %.....check if undulator was on to run full reconstruction...........
    if und_check==0 %if undulator was off, skip to line 450
        msg_und='Undulator and/or laser are not used. \n\n';
        %note: if undulator is unused X vs S0 is plotted, not Y vs S0
        I_lin=sum(I,2);
        I_lin=smooth(I_lin,10); I_lin=I_lin';
        s0_final=Y_new(:,1)';
        ind_tp=[];
        
        if debug_check==1
            f3=figure;
        end
    else %if undulator was on, S0 is reconstructed from X below
        msg_und='Undulator and/or laser are used. \n\n';
    
        %.....find maxima and minima of the sine distribution............
        disp('- find maxima and minima of the sine distribution')
        [I_max, I_imax]=max(I,[],2); %max I for each y-value
        x_max=zeros(1,size(Y_new,1)); y_max=zeros(1,size(Y_new,1));
        for l=1:size(Y_new,1) %x- and y-coordinate for each max I value
        x_max(l)=X_new(l,I_imax(l)); y_max(l)=Y_new(l,I_imax(l));
        end
        x_max(I_max==0)=NaN; y_max(I_max==0)=NaN; I_max(I_max==0)=NaN; 
        I_max=I_max';
        
        toc 
          
        %...........find amplitude of sine curve...................
        disp('- find amplitude of the sine distribution')
        j=1;
        %identify pixels with intensity equal to or larger than 75% of max
        %intensity
        for l=1:size(Y_new,1)
            mI80=find(I(l,:)>0.75*max(I(:)));
            if ~isempty(mI80)
                xmI80(j)=mean(X_new(l,mI80));
                ymI80(j)=mean(Y_new(l,mI80));
                j=j+1;
            end
        end
        %define amplitude of sine curve to fit as theoretically predicted
        %value; this can be changed by checking idealparcheck2 checkbox
        x_amp=abs(A)*L;
        disp(['theoretical value for sine amplitude: x_amp=' num2str(x_amp) 'm'])
        
        %if debug option is set: plot the screen image and the points for
        %sine amplitude finding onto an additional plot for debugging
        if debug_check==1
            f3=figure();
            subplot(1,2,1)
            z3=zeros(size(xmI80))+max(I(:));
            surf(X_new,Y_new,I,'LineStyle','None')
            xlabel('X [m]'); ylabel('Y [m]');
            axis tight; view([0 90]); colorbar;
            title('Measured screen image');
            hold on
            plot3(xmI80,ymI80,z3,'ro')
            hold off
        end
        
        toc 
        
        %...........fit a sine curve to the data .......................
        disp('- fit sine curve to image')
        %identify suitable points of I_max to fit a sine to
        indfit=find(I_max>0.1*max(I_max));
        xfit=x_max(indfit); yfit=y_max(indfit);

        if manualpar_check2==1
            %if sufficiently many pixels are found for xm80 and ym80 and they are well clustered,
            %choose their mean x-position as the sine signal amplitude;
            %otherwise choose sine amplitude manually by indicating a
            %position along x in the screen signal image in Figure 1
            stdmI80=std(abs(xmI80));
            if numel(xmI80)<3 || stdmI80>=mean(abs(xmI80))/2
                 dcm_obj1 = datacursormode(f); 
                 datacursormode on;
                 msg_x1=msgbox('Click the correct amplitude reference point, then press Return.');
                 waitfor(msg_x1)
                 cursor_info1=getCursorInfo(dcm_obj1);
                 x_amp=cursor_info1.Position(1);
                 disp(['manual value for sine amplitude: x_amp=' num2str(x_amp)])
            else
                x_amp=mean(abs(xmI80));
                disp(['use updated value from image for sine amplitude: x_amp=' num2str(x_amp)])
            end
        end
        
        %Various other options for defining amplitude and wavelength of the
        %sine curve to be fit to the signal
        %for each case, best parameters for a sine fit are determined
        if idealpar_check==1 && (idealpar_check2==1 || manualpar_check2==1)
            %both wavelength and amplitude pre-defined by ideal value from
            %input parameters and theory
            fitsin = @(b,x)  x_amp.*sin(k*x + 2*pi/b(1));
            fcn = @(b) sum((fitsin(b,yfit) - xfit).^2);
            smin=fminsearch(fcn,[-1]);
            dsin=s.dsin;
        elseif idealpar_check==1
            %only wavelength pre-defined by ideal value
            fitsin = @(b,x)  b(1).*sin(k*x + 2*pi/b(2));
            fcn = @(b) sum((fitsin(b,yfit) - xfit).^2);
            smin=fminsearch(fcn,[x_amp;-1]);
            dsin=s.dsin;
        elseif (idealpar_check2==1 || manualpar_check2==1)
            %only amplitude pre-defined by ideal value
            fitsin = @(b,x)  x_amp.*sin(2*pi/b(1)*x + 2*pi/b(2));
            fcn = @(b) sum((fitsin(b,yfit) - xfit).^2);
            smin=fminsearch(fcn,[(2*pi/k);-1]);
            dsin=s.dsin;
        else
            %nothing pre-defined by ideal value; an initial guess for the
            %sine curve amplitude is found from the manually estimated
            %value
            fitsin = @(b,x)  b(1).*(sin(2*pi/b(2)*x + 2*pi/b(3)));
            fcn = @(b) sum((fitsin(b,yfit) - xfit).^2);
            smin=fminsearch(fcn,[x_amp;(2*pi/k);-1]);
            dsin=s.dsin;
        end
        %define s0 as ysin and calculate the fitted sine for it
        ysin=min(Y_new(:)):dsin:max(Y_new(:));
        xsin=fitsin(smin,ysin);
        %plot the fitted sine on top of the screen image
        axes(ha1)
        hold on
        zsin=zeros(size(xsin))+max(I(:));
        plot3(xsin,ysin,zsin,'g')
        hold off
        legend('screen image','reference point for multi shot reconstr.','sine fit')
            
        
        % if debug option is set: add the fitted sine to the debugging plot
        if debug_check==1
            figure(f3)
            subplot(1,2,1)
            hold on
            zsin=zeros(size(xsin))+max(I(:));
            plot3(xsin,ysin,zsin,'g')
            hold off
        end
        
        toc    
        
        %.........find turning point regions in the fitted sine.........
        disp('- find turning points of fitted sine')
        %define turning point regions by their relative distance in the
        %x-direction
        dxsin=diff(xsin)/max(abs(diff(xsin)));
        ind_tp=find(abs(dxsin)<0.75);
        xsin_tp=xsin(ind_tp); ysin_tp=ysin(ind_tp);
        %organise chosen points into individual regions for each turning
        %point with a start and end point
        indendtp=find(diff(ind_tp)>1);
        indstarttp=indendtp+1;
        indendtp(indendtp==1)=[];
        indstarttp(indstarttp==numel(ind_tp))=[];
        %add potential start and end points for regions that were missed at
        %the top and bottom of the image
        if numel(indstarttp)>0 && numel(indendtp)>0
        if indstarttp(1)>indendtp(1)
            %disp('added extra start point for turning point regions')
            indstarttp=[1 indstarttp];
        end
        if indendtp(numel(indendtp))<indstarttp(numel(indstarttp))
            %disp('added extra end point for turning point regions')
            indendtp=[indendtp numel(ind_tp)];
        end
        indstarttp=unique(indstarttp); %delete any double entries
        indendtp=unique(indendtp);
        start_tp=ind_tp(indstarttp);
        end_tp=ind_tp(indendtp);
        else
            start_tp=[];
            end_tp=[];
            ind_tp=[0];
        end
        %if debug option is set: add turning point regions and their start
        %/ end points to the debugging plot
        if debug_check==1
            figure(f3)
            subplot(1,2,1)
            hold on
            zsin_tp=zeros(size(xsin_tp))+max(I(:))+10;
            plot3(xsin_tp,ysin_tp,zsin_tp,'yo')
            zsin=zeros(size(xsin))+max(I(:));
            plot3(xsin(start_tp),ysin(start_tp),zsin(start_tp),'r*')
            plot3(xsin(end_tp),ysin(end_tp),zsin(end_tp),'k*')
            hold off
            legend('screen image','I>75%max(I)','sine fit','turning point regions',...
                'start of each turning point region','end of each turning point region');
        end
        
        toc    
          
        %........find the intensity in the image by integrating along
        %the trace of the fitted sine curve......................
        disp('- integrate intensities between neighbouring points along the sine fit')    
        %sum intensities of pixels along the fitted sine trace between
        %neighbouring data points
        B=zeros(size(ysin));
        tracedis=zeros(size(ysin));
        sizeBtemp=zeros(size(ysin));
        
        %redefines I to consider sum over spread of sine curve along x for
        %each pixel along y
        I_new=zeros(size(I));
        for l=1:numel(Y_new(:,1))
            s0m=sum(X_new(l,:).*I(l,:))./sum(I(l,:));
            s0R=sqrt(sum(I(l,:).*(X_new(l,:)-s0m).^2)./(sum(I(l,:))-1));
            [~, il1]=min(abs(X_new(l,:)-(s0m-3*s0R)));
            [~, il2]=min(abs(X_new(l,:)-(s0m+3*s0R)));
            [~, il3]=min(abs(X_new(l,:)-(s0m)));
            I_new(l,il3)=sum(I(l,il1:il2));
        end    
        for l=1:numel(ysin)-1
            if any(l~=ind_tp)
                [~,~,Btemp]=improfile(X_new,Y_new,I,xsin(l:l+1),ysin(l:l+1),'bilinear');
                B(l)=sum(Btemp)-Btemp(numel(Btemp));
                sizeBtemp(l)=numel(Btemp);
                tracedis(l)=sqrt((xsin(l+1)-xsin(l)).^2+(ysin(l+1)-ysin(l)).^2);
            end
        end
        
        %definition of the final reconstructed bunch profile
        I_lin=B;
        s0_final=ysin;
        toc
        
        %........adapt intensity at turning point regions to follow a
        %straight line between their surrounding linear-region neighboring
        %points.................................
        %note: this is necessary as the intensity at the turning point
        %regions is strongly underestimated due to intrinsic spread of the e-beam 
        %in the horizontal direction
        disp('- correct for intensities around the turning point regions of the sine curve')
        grad=zeros([1,numel(start_tp)]);
        zcross=zeros([1,numel(start_tp)]);
        for l=1:numel(start_tp)
            %define relevant neighbouring points
            if start_tp(l)>3 && end_tp(l)<numel(s0_final)-3
                iIlin1=start_tp(l)-3:start_tp(l)-1; 
                iIlin2=end_tp(l)+1:end_tp(l)+3;
            elseif end_tp(l)<max(ind_tp)-3
                iIlin1=1:start_tp(l)-1; 
                iIlin2=end_tp(l)+1:end_tp(l)+3;
            elseif start_tp(l)>3
                iIlin1=start_tp(l)-3:start_tp(l)-1; 
                iIlin2=end_tp(l)+1:numel(s0_final);
            else
                iIlin1=1:start_tp(l)-1; 
                iIlin2=end_tp(l)+1:numel(s0_final);
            end
            if isempty(iIlin1)
                iIlin1=1;
            end
            if isempty(iIlin2)
                iIlin2=numel(s0_final);
            end
            Ilin1=mean(I_lin(iIlin1));
            Ilin2=mean(I_lin(iIlin2));
            xlin1=mean(s0_final(iIlin1));
            xlin2=mean(s0_final(iIlin2));
            %find fitted line between the two defined neighbour points
            grad(l)=(Ilin2-Ilin1)/(xlin2-xlin1); %gradient
            zcross(l)=Ilin2-grad(l)*xlin2; %zero crossing
            
             I_lin(start_tp(l):end_tp(l))=grad(l)*s0_final(start_tp(l):end_tp(l))+zcross(l);
        end
        
        I_lin(isnan(I_lin))=0;
        
        %if debug option is set: plot reconstructed bunch profile with turning 
        %point regions marked
        if debug_check==1   
            figure(f3);
            subplot(1,2,2)
            plot(s0_final,I_lin)
            xlabel('S0 [m]'); ylabel('I [arb. units]');
            title('Reconstructed bunch profile');
            hold on
            plot(s0_final(ind_tp),I_lin(ind_tp),'ro');
            hold off;
        end
    
        toc    

    end
     
    %subtract any remaining background in the bunch profile, assuming that
    %the very sides of the profile are equal to zero
    %****turned off in current version, may need to be re-activated depending on sample
    %type
%     disp('- subtract any remaining background noise level')
%     I_linmin=mean([I_lin(s0_final<=min(s0_final)+s.lambda) I_lin(s0_final>=max(s0_final)-s.lambda)]);
%     disp(['remaining background signal of amplitude ' num2str(I_linmin) ' subtracted'])
%     I_lin=I_lin-5/4*I_linmin;
%     I_lin(I_lin<0)=0;
  
    
    %define arrays for linear sections of image only
    s0_final_lin=s0_final;
    I_lin_lin=I_lin; 
    if any(ind_tp~=0)
    I_lin_lin(ind_tp)=NaN;
    end
    
    disp('- calculate s0RMS and plot profile')
    %if debug option is set: add peak data points to correct and their 
    %corrected values to the bunch profile plot
    if debug_check==1
        figure(f3)
        if und_check==1
            subplot(1,2,2)
        end
        hold on
        plot(s0_final,I_lin,'k')
        hold off
    end
    
    %plot the peak-corrected long. bunch profile
    axes(ha2);
    plot(s0_final,I_lin)
    xlabel('S0 [m]'); ylabel('I [arb. units]');
    title('Reconstructed long. bunch profile');

    %if debug option is set: add turning point data points to long. bunch profile
    if debug_check==1 && any(ind_tp~=0)
        axes(ha2);
        hold on
        plot(s0_final(ind_tp),I_lin(ind_tp),'ko')
        legend('long. bunch profile','corrected turning point regions')
        hold off
    end

    %........calculate the mean and standard deviation of the recovered
    %profile........................................................
    s0mean=sum(s0_final.*I_lin)./(sum(I_lin));
    s0var=sum(I_lin.*(s0_final-s0mean).^2)./(sum(I_lin)-1);
    s0RMS=sqrt(s0var);
    disp(['RMS length = ' num2str(s0RMS) 'm'])
    
    %add s0RMS value to plot
    axes(ha2);
    xcord=xlim(ha2);
    ycord=ylim(ha2);
    text(xcord(1)+0.1*abs(xcord(1)),ycord(2)-0.1*abs(ycord(2)),['RMS length = ' num2str(s0RMS) 'm'],'Parent',ha2);
    hold off
    
    %if debug option is set: add RMS width value to debug plot
    if debug_check==1
        figure(f3);
        if und_check==1
            subplot(1,2,2)
        end
        hold on;

        legend('reconstructed profile','corrected turning point regions','reconstructed profile after subtracting background');
        axis tight;
        xlabel('s_0 [m]'); ylabel('I [a.u.]');
        xcord=xlim(gca);
        ycord=ylim(gca);
        text(xcord(1)+0.1*abs(xcord(1)),ycord(2)-0.1*abs(ycord(2)),['RMS length = ' num2str(s0RMS) 'm'],'Parent',gca);
        hold off;
    end
    
    toc 
    disp(' ')
 
end

%---------------CALLBACK FOR MICROBUNCH BUTTON:h_microbunchbutton----------
%------------------------------------------------------------------------
%prompts user to choose two points on the graph as region of interest
%between which the length and distance of microbunching structures is
%measured
function h_microbunchbutton_Callback(source,eventdata)
    disp('MICROBUNCH MEASUREMENT')
    tic %measure time for reconstruction algorithm to complete
    disp('- quantify size and distance of microbunching structures in the profile')
    
    meanFWHM=0; meanFWHMhole=0;
    
    %......question dialog: choose whether to automatically detect
    %microbunch positions or manually define them.................
    choice_mb=questdlg('Please choose how to detect microbunches in the reconstructed profile.',...
        'Microbunching detection','automatic','manual','automatic');
    
    %......automatic microbunch detection..................
    if strcmp(choice_mb,'automatic')
        %......choose a region of interest in the profile inside of which to
        %calculate the width and distance of microbunches.................
        %input for the beginning of the region of interest
        dcm_obj1 = datacursormode(f);
        datacursormode on;
        msg_x1=msgbox('Please choose a region of interest for measuring microbunches. Click the beginning of the region of interest, then press Return.');
        waitfor(msg_x1)
        cursor_info1=getCursorInfo(dcm_obj1);
        llimit=cursor_info1.Position(1);
        %input for the end of the region of interest
        dcm_obj2 = datacursormode(f);
        datacursormode on;
        msg_x2=msgbox('Click the end of the region of interest, then press Return.');
        waitfor(msg_x2)
        cursor_info2=getCursorInfo(dcm_obj2);
        ulimit=cursor_info2.Position(1);
        %error if ulimit<llimit
        if ulimit<llimit
            warndlg('Upper limit is below lower limit of region of interest! Rechoose values.'...
                    ,'Warning');
        end
        %....choose all data points within the region of interest (ROI)...
        xregion=s0_final(s0_final>=llimit & s0_final<=ulimit);
        Iregion=I_lin(s0_final>=llimit & s0_final<=ulimit);
        %smooth the ROI profile to facilitate finding peaks and troughs
        xregion2=xregion; Iregion2=smooth(Iregion,10); Iregion2=Iregion2';
        
        %...find maxima and minima of profile to use for definition of individual
        %microbunches....
        [xregionmax, ixregionmax]=findpeaks(Iregion2,'MinPeakProminence',0.2*max(Iregion));
         if numel(ixregionmax)<=1
            warndlg('No microbunching structure was detected. Check again with manual mode.'...
                ,'Warning');
         end
        [xregionmin, ixregionmin]=findpeaks(-Iregion2,xregion2,'MinPeakDistance',...
            2/3*mean(abs(diff(xregion2(ixregionmax)))),'MinPeakProminence',0.2*max(Iregion));

        %....add profile end points as minima for peaks at each end of
        %profile....
        addhl=0; addh1=0;     
        if xregion2(ixregionmax(1))<ixregionmin(1)
            ixregionmin=[xregion2(1) ixregionmin];         
            xregionmin=[-Iregion2(1) xregionmin];
            addh1=1;
        end
        if xregion2(ixregionmax(numel(ixregionmax)))>ixregionmin(numel(ixregionmin))
            ixregionmin=[ixregionmin xregion2(numel(xregion2))]; 
            xregionmin=[xregionmin -Iregion2(numel(Iregion2))];
            addhl=1;
        end
        %...divide the bunch profile into regions of microbunches and holes
        %between them according to the extrema above....
        mpeaks=cell([1 numel(ixregionmin)-1]);
        mpeaksI=cell([1 numel(ixregionmin)-1]);
        for l=1:numel(ixregionmin)-1
            mpeaks{l}=xregion2(xregion2>=ixregionmin(l) & xregion2<ixregionmin(l+1));
            mpeaksI{l}=Iregion2(xregion2>=ixregionmin(l) & xregion2<ixregionmin(l+1));
        end
        mholes=cell([1 numel(ixregionmax)-1]);
        mholesI=cell([1 numel(ixregionmax)-1]);
        for l=1:numel(ixregionmax)-1
            mholes{l}=xregion2(xregion2>=xregion2(ixregionmax(l)) & xregion2<xregion2(ixregionmax(l+1)));
            mholesI{l}=Iregion2(xregion2>=xregion2(ixregionmax(l)) & xregion2<xregion2(ixregionmax(l+1)));
        end

    %..............manual microbunch detection............................    
    else
        %.....use full profile and smooth out for easier peak detection...
        xregion2=s0_final; Iregion2=smooth(I_lin,10); Iregion2=Iregion2';
        
        %.....manually define all peak positions...........
        dcm_obj1 = datacursormode(f);
        datacursormode on;
        msg_x1=msgbox('Mark all microbunch peaks, then press Return.');
        waitfor(msg_x1)
        cursor_info_peaks=getCursorInfo(dcm_obj1);
        
        %......manually define all trough positions......
        dcm_obj1 = datacursormode(f);
        datacursormode on;
        msg_x1=msgbox('Mark all microbunch troughs, then press Return.');
        waitfor(msg_x1)
        cursor_info_holes=getCursorInfo(dcm_obj1);
        
        %....extract peak and hole centres from defined positions.....
        xpeak=zeros([1 numel(cursor_info_peaks)]);
        Ipeak=zeros([1 numel(cursor_info_peaks)]);
        for l=1:numel(cursor_info_peaks)
            xpeak(l)=cursor_info_peaks(l).Position(1);
            Ipeak(l)=cursor_info_peaks(l).Position(2);
        end
        xhole=zeros([1 numel(cursor_info_holes)]);
        Ihole=zeros([1 numel(cursor_info_holes)]);
        for l=1:numel(cursor_info_holes)
            xhole(l)=cursor_info_holes(l).Position(1);
            Ihole(l)=cursor_info_holes(l).Position(2);
        end
        [xpeak, Ip]=sort(xpeak);
        [xhole, Ih]=sort(xhole);
        Ipeak=Ipeak(Ip);
        Ihole=-Ihole(Ih);
        
        %....add profile end points as minima for peaks at each end of
        %profile....
        addh1=0; addhl=0;
        if xpeak(1)<xhole(1)
            xhole=[xregion2(1) xhole]; 
            Ihole=[-Iregion2(1) Ihole];
            addh1=1;
        end
        if xpeak(numel(xpeak))>xhole(numel(xhole))
            xhole=[xhole xregion2(numel(xregion2))]; 
            Ihole=[Ihole -Iregion2(numel(Iregion2))];
            addhl=1;
        end
        Ihole = -Ihole;
                
        %...divide the bunch profile into regions of microbunches and holes
        %between them according to the extrema above....
        mholes=cell([1 numel(xpeak)-1]);
        mholesI=cell([1 numel(xpeak)-1]);
        for l=1:numel(xpeak)-1 
            mholes{l}=xregion2(xregion2>=xpeak(l) & xregion2<xpeak(l+1));
            mholesI{l}=Iregion2(xregion2>=xpeak(l) & xregion2<xpeak(l+1));
        end
        mpeaks=cell([1 numel(xhole)-1]);
        mpeaksI=cell([1 numel(xhole)-1]);
        for l=1:numel(xhole)-1
            mpeaks{l}=xregion2(xregion2>=xhole(l) & xregion2<xhole(l+1));
            mpeaksI{l}=Iregion2(xregion2>=xhole(l) &xregion2 <xhole(l+1));
        end
    end 

    %.................calculate the FWHM width of each microbunch.............
    %(note: it was tested to calculate the FWHM from the RMS value of
    %each peak assuming they are Gaussian, but this did not produce
    %good result; a fit of a Gaussian to each peak is also not sufficient)
    if exist('mpeaks','var')==1
        minpeaksI=cell([1 numel(mpeaks)]);
        maxpeaksI=cell([1 numel(mpeaks)]);
        centpeaks=zeros([1 numel(mpeaks)]);
        FWHM1=zeros([1 numel(mpeaks)]);
        FWHM1I=zeros([1 numel(mpeaks)]);
        FWHM2=zeros([1 numel(mpeaks)]);
        FWHM2I=zeros([1 numel(mpeaks)]);
        FWHM=zeros([1 numel(mpeaks)]);
        for l=1:numel(mpeaks)
            if addh1==1 && l==1
                minpeaksI{l}=mean(mpeaksI{l}(numel(mpeaksI{l})-5:numel(mpeaksI{l})));
            elseif addhl==1 && l==numel(mpeaks)
                minpeaksI{l}=mean(mpeaksI{l}(1:6));
            else
                minpeaksI{l}=mean([mpeaksI{l}(numel(mpeaksI{l})-5:numel(mpeaksI{l})) mpeaksI{l}(1:6)]); %minimum of peak
            end
            if strcmp(choice_mb,'automatic')
                maxpeaksI{l}=max(mpeaksI{l}); %maximum of peak
            else
                maxpeaksI{l}=Ipeak(l);
            end
            centpeaks(l)=sum(mpeaks{l}.*mpeaksI{l})./sum(mpeaksI{l});%mean(mpeaks{l}); %central position of peak
            lev50=(maxpeaksI{l}-minpeaksI{l})/2+minpeaksI{l}; %half maximum intensity
            signI=sign(mpeaksI{l}-lev50);
            FWHMind=find(diff(signI)~=0); %coordinate of half maximum level position
            interp1=(lev50-mpeaksI{l}(FWHMind(1)))/(mpeaksI{l}(FWHMind(1)+1)-mpeaksI{l}(FWHMind(1))); %interpolation term
            FWHM1(l)=mpeaks{l}(FWHMind(1))+interp1*(mpeaks{l}(FWHMind(1)+1)-mpeaks{l}(FWHMind(1))); %left half maximum position
            FWHM1I(l)=mpeaksI{l}(FWHMind(1))+interp1*(mpeaksI{l}(FWHMind(1)+1)-mpeaksI{l}(FWHMind(1))); %intensity at left half maximum position
            interp2=(lev50-mpeaksI{l}(FWHMind(numel(FWHMind))))/(mpeaksI{l}(FWHMind(numel(FWHMind))+1)-mpeaksI{l}(FWHMind(numel(FWHMind)))); %interpolation term
            FWHM2(l)=mpeaks{l}(FWHMind(numel(FWHMind)))+interp2*(mpeaks{l}(FWHMind(numel(FWHMind))+1)-mpeaks{l}(FWHMind(numel(FWHMind)))); %right half maximum position
            FWHM2I(l)=mpeaksI{l}(FWHMind(numel(FWHMind)))+interp2*(mpeaksI{l}(FWHMind(numel(FWHMind))+1)-mpeaksI{l}(FWHMind(numel(FWHMind)))); %intensity at right half maximum position
            FWHM(l)=FWHM2(l)-FWHM1(l); %full width at half maximum
        end
        FWHM=FWHM(isfinite(FWHM)); %delete all NaN and Inf values
        meanFWHM=mean(FWHM(FWHM>0));
        medianFWHM=median(FWHM(FWHM>0));
        stdFWHM=std(FWHM(FWHM>0));
    else
        warndlg('No microbunching structure was detected. If not done yet, check again with small microbunches feature on.'...
                ,'Warning');
    end
    disp(['mean of FWHM: ' num2str(meanFWHM) 'm, median of FWHM: ' num2str(medianFWHM) 'm, std. deviation of FWHM: ' num2str(stdFWHM) 'm'])
    
    %.............calculate the size of the holes as the distance between the FWHM
    %positions.......................................................
    if exist('mholes','var')==1 
        minholesI=cell([1 numel(mholes)]);
        maxholesI=cell([1 numel(mholes)]);
        centholes=zeros([1 numel(mholes)]);
        FWHM1holes=zeros([1 numel(mholes)]);
        FWHM1Iholes=zeros([1 numel(mholes)]);
        FWHM2holes=zeros([1 numel(mholes)]);
        FWHM2Iholes=zeros([1 numel(mholes)]);
        FWHMhole=zeros([1 numel(mholes)]);
        for l=1:numel(mholes)
            if strcmp(choice_mb,'automatic')
                minholesI{l}=min(mholesI{l}); %min of hole
            else
                minholesI{l}=Ihole(l);
            end
            maxholesI{l}=mean([mholesI{l}(numel(mholesI{l})-5:numel(mholesI{l})) mholesI{l}(1:6)]);%maximum of hole
            centholes(l)=mean(mholes{l}); %central position of hole
            lev50holes=(maxholesI{l}-minholesI{l})/2+minholesI{l}; %half maximum intensity
            signIholes=sign(mholesI{l}-lev50holes);
            FWHMindholes=find(diff(signIholes)~=0); %coordinate of half maximum level position
            
            interp1=(lev50holes-mholesI{l}(FWHMindholes(1)))/(mholesI{l}(FWHMindholes(1)+1)-mholesI{l}(FWHMindholes(1))); %interpolation term
            FWHM1holes(l)=mholes{l}(FWHMindholes(1))+interp1*(mholes{l}(FWHMindholes(1)+1)-mholes{l}(FWHMindholes(1))); %left half maximum position
            FWHM1Iholes(l)=mholesI{l}(FWHMindholes(1))+interp1*(mholesI{l}(FWHMindholes(1)+1)-mholesI{l}(FWHMindholes(1))); %intensity at left half maximum position
            interp2=(lev50holes-mholesI{l}(FWHMindholes(numel(FWHMindholes))))/(mholesI{l}(FWHMindholes(numel(FWHMindholes))+1)...
                -mholesI{l}(FWHMindholes(numel(FWHMindholes)))); %interpolation term
            FWHM2holes(l)=mholes{l}(FWHMindholes(numel(FWHMindholes)))+interp2*(mholes{l}(FWHMindholes(numel(FWHMindholes))+1)...
                -mholes{l}(FWHMindholes(numel(FWHMindholes)))); %right half maximum position
            FWHM2Iholes(l)=mholesI{l}(FWHMindholes(numel(FWHMindholes)))+interp2*(mholesI{l}(FWHMindholes(numel(FWHMindholes))+1)...
                -mholesI{l}(FWHMindholes(numel(FWHMindholes)))); %intensity at right half maximum position
            FWHMhole(l)=FWHM2holes(l)-FWHM1holes(l); %full width at half maximum
        end
        FWHMhole=FWHMhole(isfinite(FWHMhole));
        meanFWHMhole=mean(FWHMhole(FWHMhole>0));
        medianFWHMhole=median(FWHMhole(FWHMhole>0));
        stdFWHMhole=std(FWHMhole(FWHMhole>0));
    else
        meanFWHMhole=NaN;
        medianFWHMhole=NaN;
        FWHMhole=NaN;
        centholes=NaN;
    end
    disp(['mean of distance: ' num2str(meanFWHMhole) 'm, median of distance: ' num2str(medianFWHMhole) 'm, std. deviation of distance: ' num2str(stdFWHMhole) 'm'])    

    %..............add meanFWHMhole and meanFWHM values to plot.................
    str=sprintf(['Microbunching:\nWidth (FWHM) = ' num2str(meanFWHM) 'm\nDistance = ' num2str(meanFWHMhole) 'm']);
    xcord=xlim(ha2);
    ycord=ylim(ha2);
    t1=findobj(ha2,'Type','Text');
    if ~isempty(t1) %replace previous microbunching text box with a new one
        for l=1:size(t1,1)
            if size(t1(l).String,1)>1
            delete(t1(l));
            end
        end
    end
    text(xcord(1)+0.1*abs(xcord(1)),ycord(2)-0.2*abs(ycord(2)),str,'Parent',ha2);

    %if debug option is set: plot two extra plots for microbunching debugging details: one with
    %individual microbunch peaks, their smoothing spline, the max and
    %min points and the FWHM positions for each peak, another one showing
    %the calculated microbunch width and distance for each identified peak
    if debug_check==1
        figure();
        subplot(2,1,1)
        % Generate dummy info for plot handles "h" to create legend
        h = zeros(5,1);
        h(1) = plot(FWHM1(l),FWHM1I(l),'ro'); hold on;
        h(2) = plot(FWHM1(l),FWHM1I(l),'rx'); 
        h(3) = plot(FWHM1(l),FWHM1I(l),'yo');
        h(4) = plot(FWHM1(l),FWHM1I(l),'go');
        
        hold on
        if exist('mpeaks','var')==1
            for l=1:size(mpeaks,2)
                plot(mpeaks{l},mpeaksI{l});
                plot(FWHM1(l),FWHM1I(l),'ro',FWHM2(l),FWHM2I(l),'ro');
            end
        end
        plot(s0_final,I_lin,'y');
        if exist('mholes','var')==1
            for l=1:size(mholes,2)
                plot(FWHM1holes(l),FWHM1Iholes(l),'rx',FWHM2holes(l),FWHM2Iholes(l),'rx');
            end
        end
        if strcmp(choice_mb,'automatic')
           plot(xregion2(ixregionmax),xregionmax,'yo',ixregionmin,-xregionmin,'go');
        else
            plot(xpeak,Ipeak,'yo',xhole,Ihole,'go');
        end
        hold off;
        xlabel('s0 [m]'); ylabel('I [arb. units]');
        legend('FWHM points for definition of microbunch width','FWHM points for definition of microbunch distance',...
            'microbunch centre','centre point between neighbouring microbunches')
        axis tight;
        subplot(2,1,2)
        plot(centpeaks,FWHM,'bo-',centholes,FWHMhole,'ro-');
        xlabel('s0 [m]'); ylabel('FWHM width[m]');
        legend(['microbunch width'],['microbunch distance']);
    end     
    
    toc
    
end

%---------------CALLBACK FOR MEASURE BUTTON:h_measbutton-----------------
%------------------------------------------------------------------------
%prompts user to choose two points on the graphs and measures their distance in s0
function h_measbutton_Callback(source,eventdata)
     %..........input for first data point..........................
     dcm_obj1 = datacursormode(f);
     datacursormode on;
     msg_x1=msgbox('Click the first data point, then press Return.');
     waitfor(msg_x1)
     cursor_info1=getCursorInfo(dcm_obj1);
     x1=cursor_info1.Position(1);
     
     %.........input for second data point...........................
     dcm_obj2 = datacursormode(f);
     datacursormode on;
     msg_x2=msgbox('Click the second data point, then press Return.');
     waitfor(msg_x2)
     cursor_info2=getCursorInfo(dcm_obj2);
     x2=cursor_info2.Position(1);
     
     %......calculate and output distance between the two points..........
     distance=abs(x2-x1);
     msgbox(['The distance between the two points is ' num2str(distance) 'm.']);
end

%-----------------CALLBACK FOR SAVE DATA BUTTON:h_savebutton---------------
%--------------------------------------------------------------------------
%saves a log file with useful data and two output files with the reconstructed
%2D and lineout image data
function h_savebutton_Callback(source,eventdata)
    %.....open a log file to write important parameters in...........
    c=clock;
    %define file name depending on if data is based on single shot or multi
    %shot reconstruction
    if multi=='y'
        fileID=['Attoscope_BeamImage_MultiShot_' date '-' num2str(c(4)) '-' num2str(c(5))];
    else
        fileID=['Attoscope_BeamImage_' date '-' num2str(c(4)) '-' num2str(c(5))];
    end
    fid=fopen(['log-' fileID '.txt'],'at');
    %........writes input parameters into log file................
    %define data to write into log file depending on if data is based on
    %single shot or multi shot reconstruction
    if multi=='y'
        fprintf(fid,'Multishot-Reconstruction using files: \n');
        for j=1:numel(MultiFilename)
            fprintf(fid,[MultiFilename{j} '\n']);
        end
    else
        fprintf(fid, 'Setup input parameters: \n');
        for l=1:size(htext.name,1)
            fprintf(fid,[htext.name{l} num2str(s.(hedit.var{l})) '\n']);
        end
        fprintf(fid,['Arf=' num2str(Arf) ', krf=' num2str(krf) 'm, A=' num2str(A) ', k=' num2str(k) 'm \n']);
        fprintf(fid, ['Input Filenames: Signal: ' filename_wd ', Background: ' filename_bg '\n']);
        %.......writes pixel size of original image into log file..........
        dX=zeros([1 size(X,2)-1]);
        dY=zeros([1 size(Y,1)-1]);
        for l=1:size(X,2)-1
            dX(l)=X(1,l+1)-X(1,l);
        end
        for l=1:size(Y,1)-1
            dY(l)=Y(l+1,1)-Y(l,1);
        end
        fprintf(fid, ['Input pixel size: dX=' num2str(mean(dX)) 'm, dY=' num2str(mean(dY)) 'm \n \n']);
        %.........writes deflector and undulator status into log file..........
        fprintf(fid, [msg_def msg_und]);
        %........writes output data into log file...................
        fprintf(fid, 'From beam image: \n');
        fprintf(fid,['Beam features: mean_s0=' num2str(s0mean) 'm, s0_RMS=' num2str(s0RMS) 'm\n']);
        if meanFWHM~=0
            fprintf(fid,['Microbunching: FWHM width (mean) = ' num2str(meanFWHM) 'm, FWHM width (median) = ' num2str(medianFWHM) 'm, FWHM std. deviation = ' num2str(stdFWHM) 'm, Distance (mean) = ' num2str(meanFWHMhole)...
            'm, Distance (median) = ' num2str(medianFWHMhole) 'm, Distance std. deviation = ' num2str(stdFWHMhole) 'm (region of interest: ' num2str(llimit) 'm to ' num2str(ulimit) 'm)\n']);
        end
        fprintf(fid,['Output pixel size: ds0=' num2str(mean(diff(s0_final))) 'm\n']);
    end
    %.......open output files to write reconstructed lineout data and lineout data for linear parts 
    fid3=['output-' fileID '-lineout.txt'];
    %.......bring data into format to save; note that the reference point
    %for the image is saved as the first entry in s0_final...............
    s0_final_sav=[refp s0_final]; I_lin_sav=[0 I_lin]; s0_final_lin_sav=[0 s0_final_lin]; I_lin_lin_sav=[0 I_lin_lin];
    save(fid3,'s0_final_sav','I_lin_sav','s0_final_lin_sav','I_lin_lin_sav','-ascii')
    fprintf(fid, ['Lineout output data saved in ' fid3 '\n']);
    fclose('all');    
    disp(['- output data written into ' ['log-' fileID '.txt'] ' and ' fid3])
    
    %...............save plot on ha2.....................................
    fig=figure;
    plot(s0_final,I_lin)
    xlabel('S0 [m]'); ylabel('I [arb. units]');
    title('Reconstructed long. bunch profile');
    annotation('textbox',[0.2 0.8 0.1 0.1],'String',['s0RMS=' num2str(s0RMS) 'm'],'LineStyle','None');
    hold off
    print(fig,'-dpng','-r300',[fileID '.png']);
end


end
