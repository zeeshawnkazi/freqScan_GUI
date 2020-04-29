function [raw_rf_on, raw_rf_off, avg_scan, frequencies, d] = freqScanSweepPulsed(handles, varargin)
    % this function executes an ODMR frequency scan using pulsed
    cla(handles.axes1)
    cla(handles.axes2)
    %% define where things are saved
    % if doing auto dq, get the lower and upper folders ready
    % also define folder where data is saved and logs are saved
    if handles.auto_double_quantum == 0
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsed\lower';
        if exist(data_folder) ~= 2
            mkdir(data_folder);
        end
        
    elseif handles.auto_double_quantum == 1   
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsed\upper';
        if exist(data_folder) ~= 2
            mkdir(data_folder);
        end
    else
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsed';
    end
    log_folder  = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPulsedLog'; % averaged sweep figure is logged here
    power_sweep_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansPowerSweep'; 
    
    
    %% things you could change
    % type of lorentzian for auto fitting of the measured curve (number of dips) (usually 1, 2, 3, or 5)
    lorentzian_type = 1;
    
    if ~isempty(varargin)        
            if sum(size(varargin)) > 1
                argument = varargin{1};
                mat_file_name = argument{1};  
                data_folder = argument{2};
            else
                mat_file_name = varargin{1};
            end
    else
        mat_file_name = 'freq_scans_pulsed.mat';        
    end
    final_fig_name = 'avg_sweep_pulsed.fig';
    %% end of things you would probably want to change / housekeeping and setting up
       
    % get date and time for saving files
    dateandtime = get_dateandtime;
    
    % labels for ODMR curve plot
    fscan_xlabel = 'RF (GHz)';
    fscan_ylabel = 'normalized pl';
     
        
    %% pulse blaster constants
    LASER_STRETCH = 0; % we measured that the laser pulses are stretched by 50 ns
    LASER_DELAY = 600; % previous value was 840 ns. Not sure why this changed
    RF_STRETCH = 100; % TODO: measure this value

    % pulse blaster pins (where voltage signal controlling hardware is coming from)
    LASER_PIN = 0;
    ddq =0;
    if ddq
        RF_SWITCH_PIN = 2;
    else
        RF_SWITCH_PIN = 1;
    end
    CAMERA_PIN = 3;
    
    % get RF parameters from GUI
    freq_start                 = getN(handles.startFreq); % frequency in GHz
    freq_end                   = getN(handles.endFreq); % frequency in GHz
    num_sweeps                 = getN(handles.numSweeps);
    total_exposure_time_sec    = getN(handles.exposureTime) / 1000; % seconds. setting zero will result in minimum possible (but non-zero) exposure time
    num_freqs                  = getN(handles.num_freqs) + 1; % number of frequency steps
    display_RFstepsize(freq_start, freq_end, num_freqs)

    % if doing a RF power sweep, check that here (usually don't do this!)
    power       = getN(handles.power); %RF power in dBm (10*log10(Power in milliwatts))

    
    % get binning/imaging parameters
    binning       = getN(handles.binning);
    image_size    = getN(handles.ccd_size);

    checkbinningandimagesize(binning, image_size);
    
    % kinetic time and number of images
    exposure_time         = 0.004; % seconds! for pulsed, use the minimum exposure time and repeat that
    kin_time              = get_kin_time(exposure_time, image_size, binning); % total time to take an image (inc. readout and refresh)
    
    % get number of images per RF, and if not even, make even
    images_per_freq       = round(total_exposure_time_sec / kin_time);
    images_per_freq       = images_per_freq + mod(images_per_freq, 2); % ensure this value is even (rounds up)

    %% initialize data arrays
    avg_scan              = zeros(num_freqs, 1); % array of normalized PL - this data is plotted (ODMR curve)
    raw_rf_on             = zeros(image_size, image_size, images_per_freq / 2, num_freqs, num_sweeps); % all rf_on images
    raw_rf_off            = zeros(image_size, image_size, images_per_freq / 2, num_freqs, num_sweeps); % all rf_off images

    %% program Pulseblaster for pulsed ODMR
    t_laser               = getN(handles.t_laser);
    t_laser_actual        = t_laser - LASER_STRETCH; % 400 ns is how long the laser is stretched
    t_rf                  = getN(handles.t_rf);
    t_rf_actual           = t_rf - RF_STRETCH;
    res                   = 100;
    kin_time_ns           = kin_time * 10^9;
      
    disp('doing a single quantum measurement')
    % "single quantum": laser pulse, then rf pi pulse, then repeat.
    cycle_length          = t_laser_actual + t_rf;

    num_loops             = round(kin_time_ns / cycle_length);
    pb1 = PBInd([LASER_PIN, RF_SWITCH_PIN, CAMERA_PIN], cycle_length, res, 0); % auto stop is turned off, so we have to manually stop programming later
    PB.start_programming('PULSE_PROGRAM');
    
    
    for image = 1:(images_per_freq / 2)
        %% PICTURE 1 - RF ON
        pb1.on(CAMERA_PIN, 0, cycle_length);
        pb1.on(LASER_PIN, 0, t_laser_actual);
        pb1.on(RF_SWITCH_PIN, t_laser_actual, t_rf_actual);
        pb1.program([-LASER_DELAY, 0, 0], num_loops / 2); % program first half of duty cycle

        pb1.on(LASER_PIN, 0, t_laser_actual);
        pb1.on(RF_SWITCH_PIN, t_laser_actual, t_rf_actual);
        pb1.off(CAMERA_PIN, 0, cycle_length); % camera trigger is now off (second half of duty cycle)
        pb1.program([-LASER_DELAY, 0, 0], num_loops / 2); % program second half of duty cycle

        %% PICTURE 2 - RF OFF
        pb1.on(CAMERA_PIN, 0, cycle_length);
        pb1.on(LASER_PIN, 0, t_laser_actual);
        pb1.off(RF_SWITCH_PIN, t_laser_actual, t_rf_actual);
        pb1.program([-LASER_DELAY, 0, 0], num_loops / 2); % program first half of duty cycle

        pb1.off(CAMERA_PIN, 0, cycle_length); % camera trigger is now off (second half of duty cycle)
        pb1.on(LASER_PIN, 0, t_laser_actual);
        pb1.off(RF_SWITCH_PIN, t_laser_actual, t_rf_actual);
        pb1.program([-LASER_DELAY, 0, 0], num_loops / 2); % program second half of duty cycle

    end

    % print instructions to pulse blaster
    PB.inst_pbonly(0, 'STOP', 0, 2 * PB.MIN_INSTR_LENGTH);
    PB.stop_programming(); % since auto_stop is off, manually turn off
        
    if mod(cycle_length, res) ~= 0
        error(['please adjust t_laser (' num2str(t_laser) ') and t_rf (' num2str(t_rf) ') such that the total cycle length (' num2str(cycle_length) ') is divisble by ' num2str(res)])
    end
    
    %% initialize hardware setup 
    hardware = dq_freqScanPulsedHardware; % create instance of hardware class
    disp('initializing hardware')
    hardware.init(binning, image_size, exposure_time, images_per_freq, 2.87, power, num_freqs, num_sweeps,ddq);  % initialize instance
  
    %% loop through number of sweeps
    for i = 1:num_sweeps
        cla(handles.axes1)
        cla(handles.axes2)
        % check to see if user selected abort
        abort = getV(handles.abort);  
        if abort > 0
            disp('scan aborted')
            break;
        end
        disp(['running frequency sweep ' num2str(i) ' out of ' num2str(num_sweeps)])
        
        % initialize odmr sweep data arrays
        frequencies = linspace(freq_start, freq_end, num_freqs);
        pl_array           = zeros(num_freqs, 1);
        
        % step through RFs, take an image at each frequency
        for f = 1:length(frequencies)
            
            pause(0)
            abort = getV(handles.abort);  % check to see if user selected abort
            if abort > 0
                disp('scan aborted')
                break;
            end
            
            % change RF 
            hardware.set_freq(frequencies(f)); 
            % take data (including pulse sequence)
            images = hardware.capture();
                         
            % get RF on image and RF off image and separate them
            rawRfOnt    = images(:, :, 1:2:end); % every other image is RF on
            rawRfOfft   = images(:, :, 2:2:end); % every other image is RF off
            
            raw_rf_on(:, :, :, f, i)    = rawRfOnt; % rawRFUnclet
            raw_rf_off(:, :, :, f, i)   = rawRfOfft;
            
            % get center of integration region
            x0 = image_size / 2; 
            y0 = image_size / 2; 
            
            % total RF on images are summed up -> this summed image gives
            % us a photon count rate to compute sensitivity
            summed = sum(rawRfOnt, 3);
            
            % plot RF on
            imagesc(handles.axes2, summed);     
            colormap(handles.axes2, linspecer);
            colorbar(handles.axes2)
            title(handles.axes2,'RF on image')                   

            % display integration region
            rw = 2 * 8 / binning;
            rectangle(handles.axes2, 'Position',[x0 - (rw / 2), y0 - (rw / 2), rw, rw],...
            'LineWidth', 2, 'EdgeColor', 'red')
            hold off
            
            % for sensitivity, initialize number of camera counts
            photoelectrons = 0;
            % initialize normalized PL element (data to be plotted)
            norm = 0;
            
            % loop through each RFon/RFoff pair and get photoluminescence
            % divided by background 
            for k = 1:(images_per_freq / 2) 
                rf_on_image = rawRfOnt(:, :, k);
                rf_off_image = rawRfOfft(:, :, k);
                
                % get RF on and RF off counts in 1 um^2 centered around x0, y0
                pl = average_counts(rf_on_image, x0, y0, rw);
                bg = average_counts(rf_off_image, x0, y0, rw);
  
                % RF off counts
                photoelectrons = photoelectrons + bg;
                % normalized PL data point 
                norm = norm + pl / bg;
            end

            % normalize pl_array to 1
            norm = norm / (images_per_freq / 2); 
            % put this data in pl_array
            pl_array(f) = norm; % this is the final data point
            
            % plot freq scans
            plot(handles.axes1, frequencies(1:f), pl_array(1:f),'.-b');
            box(handles.axes1,'on')
            xlabel(handles.axes1,fscan_xlabel)
            ylabel(handles.axes1,fscan_ylabel)
            hold on
            axis(handles.axes1,'tight')
            

        end
        
        %% if user hit abort
        if abort ~= 1
            % average the data
            avg_scan = avg_scan + pl_array;
            norm_avg_scan = avg_scan / i;

            % plot freq scans
            cla(handles.axes1);
            plot(handles.axes1,frequencies, norm_avg_scan,'.-b');
            box(handles.axes1, 'on')
            xlabel(handles.axes1,fscan_xlabel)
            ylabel(handles.axes1,fscan_ylabel)
            hold(handles.axes1,'on')
            axis(handles.axes1, 'tight')
            
            % plot RF off
            imagesc(handles.axes2,summed);
            image_title = strcat('RF on Image');
            title(handles.axes2,image_title)
            hold(handles.axes2,'on')
            colorbar(handles.axes2)
            axis(handles.axes2,'tight')
            % display integration region
            rw = 2 * 8 / binning;
            rectangle(handles.axes2,'Position',[x0 - rw / 2, y0 - rw / 2, rw, rw],...
            'LineWidth', 2, 'EdgeColor', 'red')
            hold off
            
        else
            break
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % SAVE STUFF IN 'd' OBJECT
        % this script creates an object 'd' that all the data is stored in
        freq_scan_pulsed_stuff_to_save;
        disp('saving data for this sweep')
        % save data
%         save([data_folder '\' mat_file_name], 'd','-v7.3') 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    %% finishing / cleaning up    
    % kill hardware object
    hardware.kill()
            
    % fit and plot fit
    % initialize fit successful flag
    fit_successful = 0;
    
    %% try to fit!
    try   
        % get total exposure time for summed RF off image
        counts_exposure_time = exposure_time * images_per_freq / 2;
        
        [eta, dB, max_contrast, FWHM, counts, final_fit, fit_frequencies] = get_sensitivity(norm_avg_scan, counts_exposure_time, frequencies, photoelectrons, power, lorentzian_type, 0);
        plot(handles.axes1,fit_frequencies, final_fit, 'r', 'LineWidth', 2);
        hold on
        axis(handles.axes1,'tight')
        
        fit_successful = 1;
        
        % display sensitivity and FWHM 
        disp(['Sensitivity associated with this curve: ' num2str(eta * 10^9) ' nT/sqrt(Hz)'])
        disp(['FWHM of this curve: ' num2str(FWHM * 10^6) ' kHz'])
    catch me
        disp(me)
    end
    
    %% try to plot full sweep
    try 
        % plot average sweep figure (to be displayed and logged)
        figure(3); clf;
        plotaxes = gca;
        plot(plotaxes,frequencies, norm_avg_scan,'.-b');
        hold(plotaxes,'on')
        if fit_successful
            plot(plotaxes,fit_frequencies, final_fit, 'r', 'LineWidth', 2);
        end
        box(plotaxes,'on')
        xlabel(plotaxes,fscan_xlabel)
        ylabel(plotaxes,fscan_ylabel)
        hold(plotaxes,'on')
        figtitle = ['Pulsed ODMR curve: RF power = ' num2str(power) ' dBm, Num freqs = ' num2str(num_freqs) ' , Exposure time = ' num2str(exposure_time) ' s, \pi Pulse = ' num2str(t_rf) ' ns, '  num2str(num_sweeps) ' sweep(s)'];
        title(plotaxes, figtitle,'FontSize',8);
        axis(plotaxes,'tight')
    catch
        disp('a full sweep has not been run')
    end
    
    
    %% save data  
    if fit_successful
        disp('saving fit data')
        freq_scan_pulsed_FIT_stuff_to_save;  
    end
    
    if abort
        freq_scan_pulsed_stuff_to_save;
    end
    
    
    % save data for single RF power

    savefig(figure(3),[data_folder '\' final_fig_name]) 
    save([data_folder '\' mat_file_name], 'd', '-v7.3')

    figname = [dateandtime '_avgSweep.fig'];
    savefig(figure(3), [log_folder '\' figname])


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
end
