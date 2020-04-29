function varargout = freqScan(varargin)
% freqScan MATLAB code for freqScan.fig
%      freqScan, by itself, creates a new freqScan or raises the existing
%      singleton*.
%
%      H = freqScan returns the handle to a new freqScan or the handle to
%      the existing singleton*.
%
%      freqScan('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in freqScan.M with the given input arguments.
%
%      freqScan('Property','Value',...) creates a new freqScan or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before freqScan_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to freqScan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 10-Apr-2019 11:37:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @freqScan_OpeningFcn, ...
                   'gui_OutputFcn',  @freqScan_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ImagerScript_fast
% testingAndor
% ImagerScript_fast_triggering
% ImagerScript_fast_bead_motion_recording
% freqScansTriggeringCommands
    global REGION_WIDTH;    % this is one half of one side of the integration area
    REGION_WIDTH = 8;       % for the current mag, with no binning 15 pixels  = 1 um
                            % so pixel size with 2x2 binning = 133 nm, -8:8 by -8:8 box ~ 1 um^2
    %% housekeeping / starting things up
    % clean up, also start a timer
    tic
    clc
    cla(handles.axes1)
    cla(handles.axes2)
    
    % TODO: MAKE THIS OBVIOUS set this flag to 2 (by default)
    handles.auto_double_quantum = 2;
    
    % if the hamamatsu preview is open, close it! camera can't do two
    % things at once
    h = findobj('Name', 'preview_fig');
    if ~isempty(h)
        close(h)
    end
    
    % check if RF is connected, if not, it'll error!
    if getV(handles.rf_connected) ~= 1
            error('CHECK YOUR RF CONNECTIONS')
    end
    
    %%%% here are various unique things you can do
    %% long frequency scan sweep
    long_sweep = 0;
    if long_sweep
        f_start = 2.82;
        f_end   = 2.92;
        f_step  = 0.0001;
        frequencies= f_start : f_step : f_end;
        num_freqs  = getN(handles.num_freqs);
        
        num_sweeps = round(length(frequencies) / num_freqs);
        for sweep = 1:num_sweeps
            start_index = (sweep-1)*num_freqs + 1;
            freq_start  = frequencies(start_index);
            end_index   = (sweep-1)*num_freqs - 1 + num_freqs;
            freq_end    = frequencies(end_index);

            setN(handles.startFreq, freq_start);
            setN(handles.endFreq, freq_end);

            file_name = ['podmrsweep' num2str(sweep) '_' num2str(freq_start) 'to' num2str(freq_end) 'GHz.mat']
            autofocus;
            freqScanSweepPulsed(handles, file_name);            
        end   
        disp('done');
        return;
    end
    %% auto double quantum: does one  sweep for the lower ODMR curve
    % and one for the upper and saves them in their own folder ('lower' and
    % 'upper')
    auto_double_quantum =0;
    if auto_double_quantum
        disp('auto double quantum on');
            if getV(handles.pulsed)
                autofocus;
                % pulsed and double quantum
%                 14N sample in fluid cell geometry:
%                start_and_end_freqs = [2.817 2.825 2.917 2.925];

%                 15N sample in face up geometry:
%                start_and_end_freqs = [2.851 2.853 2.883 2.885];

%                 1421 sample in flow cell geometry:
%                start_and_end_freqs = [2.815 2.825 2.914 2.924];

%                 1420 sample in flow cell geometry:
               start_and_end_freqs = [2.833 2.835 2.9065 2.9085];

%                1420 sample in open air geometry
%                 start_and_end_freqs = [2.859 2.865 2.877 2.883];
%                 get_new_array = input(['current range is ' newline num2str(start_and_end_freqs) ,...
%                    newline ', if this is okay,press enter, if it is not,' newline ,...
%                    'you can input your own now: (like this [a b c d]):'])
%                 if length(get_new_array) > 1
%                     start_and_end_freqs = get_new_array;
%                 end

                % get the new start and end frequencies 
                start1  = start_and_end_freqs(1);
                end1    = start_and_end_freqs(2);

                start2  = start_and_end_freqs(3);
                end2    = start_and_end_freqs(4);

                if start1 > end1 | start2 > end2
                    error('check your ranges!')
                end
                
                disp('running two single quantum freq scans')
                % run the first set (lower)
                power = getN(handles.power);
                handles.auto_double_quantum = 0; % for auto saving to a 'lower' folder
                setN(handles.startFreq, start1); % frequency in GHz
                setN(handles.endFreq, end1); % frequency in GHz
                freqScanSweepPulsed(handles); 
                
                % run the second set (upper) 
                setN(handles.power, power +1);
                handles.auto_double_quantum = 1; % for auto saving to a 'upper' folder
                setN(handles.startFreq, start2); % frequency in GHz
                setN(handles.endFreq, end2); % frequency in GHz
                freqScanSweepPulsed(handles);
               
                setN(handles.power, power);
                disp('DONE')
                return;
            end   
    else
        handles.auto_double_quantum = 2;
    end
    
    %% spin bath driving sweep
    bath_sweep = 0;
    if bath_sweep
        disp(['running spin bath sweep'])
        bath_power_index    = 0;    % can be 0, 1, 2, 3
        bath_freq_start     = 100;   % MHz
        bath_freq_end       = 500;   % MHz
        bath_freq_step      = 100;    % kHz
        bath_freq_step      = bath_freq_step * 1e-3;

        bath_freqs = bath_freq_start : bath_freq_step : bath_freq_end;
        synth = synthUSB;
        synth.init(bath_freqs(1), bath_power_index);
        
        for b = 1:length(bath_freqs)
            if ~mod(b, 10)
                autofocus;
            end
            bath_freq = bath_freqs(b);
            disp(['bath freq = ' num2str(bath_freq) ' MHz'])
            synth.change_rf(bath_freq);
            [~, ~, ~, ~, d] = freqScanSweepPulsed(handles);
            bath = d;
            bath.bath_freq = bath_freq;
            bath.bath_power_index = bath_power_index;

            bath_folder = 'spin_bath_sweep';
            bath_file   = ['freqscansbath_' num2str(bath_freq * 1e6) 'Hz.mat'];
            save([bath_folder '/' bath_file], 'bath')
        end
    end
    %% rf pulse duration sweep
    pulse_duration_sweep = 1;
    powers = (0:-1:-10);
    folder ='E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsed';
    if pulse_duration_sweep
        disp('DOING A PULSE DURATION SWEEP!');
                pulses = [300:100:2500]; 
                for p = 1:length(powers)
                    folderforthispower = [folder '/' num2str(powers(p)) 'mW'];
                    mkdir(folderforthispower);
                    power = powers(p);
                    for p2 = pulses
                        autofocus;
                        pulse = p2;
                        disp(['rf power set to ' num2str(power) ' dBm'])
                        setN(handles.power, power);
                        disp(['rf pulse duration set to ' num2str(pulse) ' ns'])
                        setN(handles.t_rf, pulse);
                        file_name = ['freq_scans_pulsed_' num2str(power) 'dBm_' num2str(pulse) 'pipulse.mat'];
                        freqScanSweepPulsed(handles, {file_name,folderforthispower});
                    end
                end
        disp('DONE')
        return;
    end
    
    %% optical pulse duration sweep
    opulse_duration_sweep = 0;
    dq=0;
    if opulse_duration_sweep
        disp('DOING AN OPTICAL PULSE DURATION SWEEP!');
                laser_power = input('enter laser power:');
                pulse_durations = [100:100:1200]; % 900 1200 1500 1800 2100 2400
                folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsed';
                for p = 1:length(pulse_durations)
                            
                    disp(['laser pulse duration set to ' num2str(pulse_durations(p)) ' ns'])   
                    setN(handles.t_laser, pulse_durations(p));
                    if ~dq                                        
                            file_name = ['freq_scans_pulsed_' num2str(laser_power) 'mW_' num2str(pulse_durations(p)) 'opticalpulse.mat'];
                            freqScanSweepPulsed(handles, {file_name, folder});
                    else
                            start_and_end_freqs = [2.819 2.8215 2.9195 2.922];
                            % get the new start and end frequencies 
                            start1  = start_and_end_freqs(1);
                            end1    = start_and_end_freqs(2);

                            start2  = start_and_end_freqs(3);
                            end2    = start_and_end_freqs(4);

                            if start1 > end1 | start2 > end2
                                error('check your ranges!')
                            end

                            disp('running two single quantum freq scans')
                            % run the first set (lower)
                            file_name1 = ['freq_scans_pulsed_lower' num2str(laser_power) 'mW_' num2str(pulse_durations(p)) 'opticalpulse.mat'];
                            power = getN(handles.power);
                            handles.auto_double_quantum = 0; % for auto saving to a 'lower' folder
                            setN(handles.startFreq, start1); % frequency in GHz
                            setN(handles.endFreq, end1); % frequency in GHz
                            autofocus;
                            freqScanSweepPulsed(handles,file_name1); 

                            % run the second set (upper) 
                            file_name2 = ['freq_scans_pulsed_upper' num2str(laser_power) 'mW_' num2str(pulse_durations(p)) 'opticalpulse.mat'];
                            setN(handles.power, power + 4.5);
                            handles.auto_double_quantum = 1; % for auto saving to a 'upper' folder
                            setN(handles.startFreq, start2); % frequency in GHz
                            setN(handles.endFreq, end2); % frequency in GHz
                            autofocus;
                            freqScanSweepPulsed(handles,file_name2);
                            setN(handles.power, power);
                    end
                end
        disp('DONE')
        toc
        return;
    end
    
    
    %% rf power sweep
    if getV(handles.rfPowerSweep)
        start_power = getN(handles.startRFPower);
        end_power = getN(handles.endRFPower);
        power_step = getN(handles.rfPowerStep);
        current_power = getN(handles.currentRFPower);

        if power_step < 0
            error('RF power step must be postive.')
        end
        power_array = start_power:power_step:end_power;
        
        if getV(handles.pulsed)
            disp('starting pulsed rf power sweep');
        else
            disp('starting CW rf power sweep');
        end
        
            for p = power_array
                setN(handles.currentRFPower, num2str(p));
                if get(handles.pulsed, 'Value')
                    disp(['rf power =', num2str(p),' dBm '])
                    freqScanSweepPulsed(handles);
                else
                    disp(strcat('rf power =',num2str(p),' dBm'))
                    freqScanSweep(handles);
                end
            end      
    end
    
    %% NORMAL FREQUENCY SCAN SWEEP
    if getV(handles.pulsed)
        disp('doing a pulsed sweep')
        autofocus;
        freqScanSweepPulsed(handles);
    else
        disp('doing a cw sweep')
        autofocus;
        freqScanSweep(handles);
    end
        
    set(handles.rf_connected, 'Value', 0);

disp('')
disp('DONE')
toc


% --- Executes just before freqScanTriggering is made visible.
function freqScan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to freqScanTriggering (see VARARGIN)

% Choose default command line output for freqScanTriggering
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes freqScanTriggering wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = freqScan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
% global feedbackCount;
% 
% feedbackCount = 0;
% % Use digital channel to turn switch on and off. 0 = Switch on, 1 = Switch
% % Off
%  
% data_sesh = daq.createSession('ni');
% data_sesh.Rate = 10000;
% data_sesh.addDigitalChannel('Dev1', 'Port0/Line6', 'OutputOnly'); 
% data_sesh.outputSingleScan(0);
% pause(0.1);
% data_sesh.release;



function startFreq_Callback(hObject, eventdata, handles)
% hObject    handle to startFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFreq as text
%        str2double(get(hObject,'String')) returns contents of startFreq as a double


% --- Executes during object creation, after setting all properties.
function startFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function power_Callback(hObject, eventdata, handles)
% hObject    handle to power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of power as text
%        str2double(get(hObject,'String')) returns contents of power as a double


% --- Executes during object creation, after setting all properties.
function power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exposureTime_Callback(hObject, eventdata, handles)
% hObject    handle to exposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exposureTime as text
%        str2double(get(hObject,'String')) returns contents of exposureTime as a double


% --- Executes during object creation, after setting all properties.
function exposureTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endFreq_Callback(hObject, eventdata, handles)
% hObject    handle to endFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFreq as text
%        str2double(get(hObject,'String')) returns contents of endFreq as a double


% --- Executes during object creation, after setting all properties.
function endFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numImages_Callback(hObject, eventdata, handles)
% hObject    handle to numImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numImages as text
%        str2double(get(hObject,'String')) returns contents of numImages as a double


% --- Executes during object creation, after setting all properties.
function numImages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Focus.
function Focus_Callback(hObject, eventdata, handles)
% hObject    handle to Focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiSettleAndFocus;
guiSettleAndFocusSmall;
% guiSettleAndFocus_partial_ccd
% guiSettleAndFocusSmall_partial_ccd


function x_center_Callback(hObject, eventdata, handles)
% hObject    handle to x_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_center as text
%        str2double(get(hObject,'String')) returns contents of x_center as a double


% --- Executes during object creation, after setting all properties.
function x_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_center_Callback(hObject, eventdata, handles)
% hObject    handle to y_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_center as text
%        str2double(get(hObject,'String')) returns contents of y_center as a double


% --- Executes during object creation, after setting all properties.
function y_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function size_Callback(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of size as text
%        str2double(get(hObject,'String')) returns contents of size as a double


% --- Executes during object creation, after setting all properties.
function size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ccd_size_Callback(hObject, eventdata, handles)
% hObject    handle to ccd_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ccd_size as text
%        str2double(get(hObject,'String')) returns contents of ccd_size as a double


% --- Executes during object creation, after setting all properties.
function ccd_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ccd_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_center_ccd_Callback(hObject, eventdata, handles)
% hObject    handle to y_center_ccd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_center_ccd as text
%        str2double(get(hObject,'String')) returns contents of y_center_ccd as a double


% --- Executes during object creation, after setting all properties.
function y_center_ccd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_center_ccd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binning_Callback(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binning as text
%        str2double(get(hObject,'String')) returns contents of binning as a double


% --- Executes during object creation, after setting all properties.
function binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function numSweeps_Callback(hObject, eventdata, handles)
% hObject    handle to numSweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numSweeps as text
%        str2double(get(hObject,'String')) returns contents of numSweeps as a double


% --- Executes during object creation, after setting all properties.
function numSweeps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numSweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentRFPower_Callback(hObject, eventdata, handles)
% hObject    handle to currentRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentRFPower as text
%        str2double(get(hObject,'String')) returns contents of currentRFPower as a double


% --- Executes during object creation, after setting all properties.
function currentRFPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startRFPower_Callback(hObject, eventdata, handles)
% hObject    handle to startRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startRFPower as text
%        str2double(get(hObject,'String')) returns contents of startRFPower as a double


% --- Executes during object creation, after setting all properties.
function startRFPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endRFPower_Callback(hObject, eventdata, handles)
% hObject    handle to endRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endRFPower as text
%        str2double(get(hObject,'String')) returns contents of endRFPower as a double


% --- Executes during object creation, after setting all properties.
function endRFPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endRFPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rfPowerStep_Callback(hObject, eventdata, handles)
% hObject    handle to rfPowerStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rfPowerStep as text
%        str2double(get(hObject,'String')) returns contents of rfPowerStep as a double


% --- Executes during object creation, after setting all properties.
function rfPowerStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rfPowerStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in rfPowerSweep.
function rfPowerSweep_Callback(hObject, eventdata, handles)
% hObject    handle to rfPowerSweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rfPowerSweep

 
% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)
% hObject    handle to abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abort


% --- Executes on key press with focus on abort and none of its controls.
function abort_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to abort (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function abort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in startPreview.
function startPreview_Callback(hObject, eventdata, handles)
% hObject    handle to startPreview (see GCBO)
    preview_hamamatsu(handles)
     
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stopPreview.
function stopPreview_Callback(hObject, eventdata, handles)
    preview_fig = findobj('type', 'figure', 'tag', 'preview_fig');
    close(preview_fig);
    ham = findall(0,'Tag','imaq.hamamatsu1_1');
    close(ham)
 
% hObject    handle to stopPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function text26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sweepTime_Callback(hObject, eventdata, handles)
% hObject    handle to sweepTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sweepTime as text
%        str2double(get(hObject,'String')) returns contents of sweepTime as a double


% --- Executes during object creation, after setting all properties.
function sweepTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sweepTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in rf_connected.
function rf_connected_Callback(hObject, eventdata, handles)
% hObject    handle to rf_connected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rf_connected


% --- Executes on button press in pulsed.
function pulsed_Callback(hObject, eventdata, handles)
% hObject    handle to pulsed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pulsed



function t_rf_Callback(hObject, eventdata, handles)
% hObject    handle to t_rf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_rf as text
%        str2double(get(hObject,'String')) returns contents of t_rf as a double


% --- Executes during object creation, after setting all properties.
function t_rf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_rf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_laser_Callback(hObject, eventdata, handles)
% hObject    handle to t_laser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_laser as text
%        str2double(get(hObject,'String')) returns contents of t_laser as a double


% --- Executes during object creation, after setting all properties.
function t_laser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_laser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freq_split_Callback(hObject, eventdata, handles)
% hObject    handle to freq_split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_split as text
%        str2double(get(hObject,'String')) returns contents of freq_split as a double


% --- Executes during object creation, after setting all properties.
function freq_split_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in double_quantum.
function double_quantum_Callback(hObject, eventdata, handles)
% hObject    handle to double_quantum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of double_quantum


% --- Executes on button press in bathSweep.
function bathSweep_Callback(hObject, eventdata, handles)
% hObject    handle to bathSweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bathSweep



function startBathFreq_Callback(hObject, eventdata, handles)
% hObject    handle to startBathFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startBathFreq as text
%        str2double(get(hObject,'String')) returns contents of startBathFreq as a double


% --- Executes during object creation, after setting all properties.
function startBathFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startBathFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepSize_Callback(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSize as text
%        str2double(get(hObject,'String')) returns contents of stepSize as a double


% --- Executes during object creation, after setting all properties.
function stepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endBathFreq_Callback(hObject, eventdata, handles)
% hObject    handle to endBathFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endBathFreq as text
%        str2double(get(hObject,'String')) returns contents of endBathFreq as a double


% --- Executes during object creation, after setting all properties.
function endBathFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endBathFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bath_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to bath_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bath_frequency as text
%        str2double(get(hObject,'String')) returns contents of bath_frequency as a double


% --- Executes during object creation, after setting all properties.
function bath_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bath_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function powerSetting_Callback(hObject, eventdata, handles)
% hObject    handle to powerSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of powerSetting as text
%        str2double(get(hObject,'String')) returns contents of powerSetting as a double


% --- Executes during object creation, after setting all properties.
function powerSetting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to powerSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_freqs_Callback(hObject, eventdata, handles)
% hObject    handle to num_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_freqs as text
%        str2double(get(hObject,'String')) returns contents of num_freqs as a double


% --- Executes during object creation, after setting all properties.
function num_freqs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
