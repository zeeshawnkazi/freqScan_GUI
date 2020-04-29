function [raw_rf_on, raw_rf_off, avg_scan, frequencies] = freqScanSweep(handles)

    
    % some constants
    global REGION_WIDTH;
    lorentzian_type      = 2; % auto fit to this many dips (usually 2, 3 or 5)
    
    RF_STEPPED = 1; % set to 1 for individual frequency steps, 0 for continuous sweep
    % sweeps through frequencies, returns images, PL array, and
    % frequency array
    
    % get date and time for saving
    dateandtime = char(datetime);               % current date and time
    dateandtime(dateandtime==' ') = '_';        % set whitespace to underscore
    dateandtime(dateandtime=='-') = '';         % delete hyphens
    dateandtime(dateandtime==':') = '';         % delete colons
    
    if handles.auto_double_quantum == 0
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScans\lower';
        if exist(data_folder) ~= 2
            mkdir(data_folder);
        end
        
    elseif handles.auto_double_quantum == 1   
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScans\upper';
        if exist(data_folder) ~= 2
            mkdir(data_folder);
        end
    else
        data_folder = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScans';
    end
    
    log_folder  = 'E:\Dropbox\Projects\magPI\freqScan_GUI\freqScansLog';
    
    mat_file_name = 'freq_scans.mat';
    final_fig_name = 'avg_sweep.fig';
    % get RF parameters from GUI
    freq_start       = getN(handles.startFreq); % frequency in GHz
    freq_end         = getN(handles.endFreq); %frequency in GHz
    num_sweeps       = getN(handles.numSweeps);
    exposure_time    = getN(handles.exposureTime)/1000; % setting zero will result in minimum possible (but non-zero) exposure time
    num_freqs        = getN(handles.num_freqs) + 1; % number of frequencies to test
    
    if get(handles.rfPowerSweep,'Value')
        power        = getN(handles.currentRFPower);
    else
        power        = getN(handles.power); %RF power in dBm (10*log10(Power in milliwatts))
    end
    
    % get binning/imaging parameters
    binning = getN(handles.binning);
    
    binning_index = binning_pixels_to_index(binning);
    
    ccd_size_index =1;
    c_size = [64 128 256 512 1024 2048]; % possible image size options
    ccd_size = c_size(ccd_size_index);

    
    % kinetic time and number of images and number of frequencies to test (if stepping)
    kin_time = get_kin_time(exposure_time, ccd_size, binning)  ;
    sweep_time = num_freqs * 2 * exposure_time;
    images_per_freq       = 2; % takes one on, one off image per frequency
    
    % initialize data arrays
    summing_scan       = zeros(num_freqs, 1);
    indiv_scans        = zeros(num_freqs-1,num_sweeps);
    
    display_RFstepsize(freq_start, freq_end, num_freqs)
    
    if RF_STEPPED
        %% individual frequency steps
            raw_rf_on          = zeros(ccd_size,ccd_size, num_freqs, num_sweeps);
            raw_rf_off         = zeros(ccd_size,ccd_size, num_freqs, num_sweeps);
    else
        % analog sweep trigger output: 0 volts to 10 volts, with one image worth of zeros at the end for the RF off image
            max_voltage = 10; 
            output_data_sweeping = linspace(0,max_voltage,freqScanHardware.DAQ_RATE*sweep_time);
            output_data_off = linspace(0,0,freqScanHardware.DAQ_RATE*kin_time);
            output_data = [output_data_sweeping output_data_off];
            
            
            raw_rf_on          = zeros(ccd_size,ccd_size,num_images - 1, num_sweeps);
            raw_rf_off         = zeros(ccd_size,ccd_size, 1, num_sweeps);
    end
    
%     micro_autofocus
    laser_on
    for i = 1:num_sweeps
        
        if RF_STEPPED
            
            pl_array           = zeros(num_freqs, 1); 
            abort = get(handles.abort, 'Value');  % check to see if user selected abort
            if(abort == 1) 
                disp('scan aborted')
                break;
            end

            disp('initializing hardware')
            % initialize hardware setup
            hardware = dq_freqScanSteppedHardware;

            disp(['running frequency sweep ' num2str(i) ' out of ' num2str(num_sweeps)])

            hardware.init(binning_index, ccd_size_index, exposure_time, images_per_freq, 2.87, power, num_freqs, num_sweeps); % freq. is dummy variable for now
            frequencies = linspace(freq_start, freq_end, num_freqs);

            for f = 1:length(frequencies)
                abort = get(handles.abort, 'Value');  % check to see if user selected abort
                if(abort == 1) 
                    disp('scan aborted')
                     % kill hardware
                    hardware.kill()
                    break;
                end

                hardware.set_freq(frequencies(f));
                images = hardware.capture();
                
                if isempty(images) || (size(images, 3) < images_per_freq)
                    disp('error acquiring images... retrying')
                    images = hardware.capture();
                end
                
                % get RF on image and RF off image and separate them
                try
                    rawRfOnt    = images(:, :, 1); % every other image is RF on
                    rawRfOfft   = images(:, :, 2); % every other image is RF off
                catch
                    pause(1);
                end
                if num_sweeps == 1
                    raw_rf_on(:, :, f)    = rawRfOnt; % rawRFUncle_t
                    raw_rf_off(:, :, f)   = rawRfOfft;
                else
                    raw_rf_on(:, :, f, i)    = rawRfOnt; % rawRFUncle_t
                    raw_rf_off(:, :, f, i)   = rawRfOfft;
                end
                norm = 0;

                % establish center of integration region
                try
                    summed = sum(rawRfOfft, 3); % summed picture
                    % TODO: click and get spot
                    x0 = ccd_size / 2 ; 
                    y0 = ccd_size / 2 ;
%                     [x0, y0] = get_freqScan_center(summed); % calculate center based on summed image
                catch me
                    disp('invalid integration region. aborting.')
                     % kill hardware
                    hardware.kill()
                    return
                end

                % plot spot
                cla(handles.axes2);
                imagesc(handles.axes2,summed);
                colorbar(handles.axes2)
                image_title = strcat('RF Off Image - see integration region');
                
                % display integration region
                rw = 2 * REGION_WIDTH / binning;
                rectangle(handles.axes2,'Position',[x0 - rw / 2, y0 - rw / 2, rw, rw],...
                'LineWidth', 2, 'EdgeColor', 'red')
                hold off

                for k = 1:(images_per_freq / 2) % for each rf on/rf off pair
                    rf_on_image = rawRfOnt(:, :, k);
                    rf_off_image = rawRfOfft(:, :, k);

                    pl = average_counts(rf_on_image, x0, y0, rw);
                    bg = average_counts(rf_off_image, x0, y0, rw);
                    norm = norm + pl / bg;
                end
                
                photoelectrons = bg;
                norm = norm / (images_per_freq / 2); % this will normalize to 1
                pl_array(f) = pl_array(f) + norm; % this is the final data point
                
                % plot freq scans
                cla(handles.axes1);
                plot(handles.axes1,frequencies(1:f), pl_array(1:f),'.-b');
                box(handles.axes1,'on');
                xlabel(handles.axes1,'Frequency (GHz)')
                ylabel(handles.axes1,'Normalized PL')
                hold(handles.axes1,'on')
                axis(handles.axes1,'tight')
    %             xlim([frequencies(1) frequencies(end)])

            end

            if abort ~= 1
                summing_scan = summing_scan + pl_array;
                avg_scan = summing_scan / i;

                % plot freq scans
                cla(handles.axes1);
                plot(handles.axes1,frequencies, avg_scan,'.-b');
                box(handles.axes1,'on');
                xlabel(handles.axes1,'Frequency (GHz)')
                ylabel(handles.axes1,'Normalized PL')
                hold(handles.axes1,'on')
                axis(handles.axes1,'tight')
                
                % plot spot
                imagesc(handles.axes2,summed);
                image_title = strcat('RF Off Image - see integration region');
                title(handles.axes2,image_title)
                hold(handles.axes2,'on')
                colorbar(handles.axes2)

                % display integration region
                rw = 2 * REGION_WIDTH / binning;
                rectangle(handles.axes2,'Position',[x0 - rw / 2, y0 - rw / 2, rw, rw],...
                'LineWidth', 2, 'EdgeColor', 'red')
                hold off

                % kill hardware
                hardware.kill()
            else
                break
            end

        else
            %% continous sweep
            
            abort = get(handles.abort, 'Value');  % check to see if user selected abort
            if(abort == 1) 
                disp('scan aborted')
                break;
            end

            disp('initializing hardware')
            % initialize hardware setup
            hardware = freqScanHardware;
            hardware.init(binning_index, ccd_size_index, exposure_time, sweep_time, freq_start, freq_end, power);

            disp(['running frequency sweep ' num2str(i) ' of ' num2str(num_sweeps)])
            
            images = hardware.capture(output_data);

            % get RF on image and RF off image and separate them
            rawRfOnt    = squeeze(double(images(:,:,1:(num_images-1))));
            rawRfOfft   = squeeze(double(images(:,:,num_images)));

            raw_rf_on(:,:,:,i)    = rawRfOnt;
            raw_rf_off(:,:,1,i)   = rawRfOfft;
            
            % current scan is RF on image
            currentScan = raw_rf_on(:,:,:,i);
            
            
            
            % establish center of integration region
            try
                x0 = ccd_size / 2;
                y0 = ccd_size / 2;

                tempS = average_counts(currentScan, x0,y0, rw);
                photoelectrons = average_counts(raw_rf_off(:,:,1,i), x0, y0, rw);
                
            catch me
                disp('invalid integration region. aborting.')
                 % kill hardware
                hardware.kill()
                return
            end
            
            
            
            
            % average scans together
            indiv_scans(:,i) =  tempS(:,1) / photoelectrons;
            avg_scan = sum(indiv_scans,2)./i;
            frequencies = linspace(freq_start,freq_end,num_images-1);

            % plot freq scans
            cla(handles.axes1);
            axes(handles.axes1);
            plot(frequencies,avg_scan,'.-b');
            box on
            xlabel('frequency (GHz)')
            ylabel('pl')
            axis tight
            hold on


            % plot spot
            axes(handles.axes2);
            imagesc(raw_rf_on(:,:,2,i));
            image_title = strcat('RF Off Image - see integration region');
            title(image_title)
            hold on
            colorbar
            colormap gray;

            % display integration region
            rw = 2 * REGION_WIDTH;
            rectangle('Position',[x0 - REGION_WIDTH, y0 - REGION_WIDTH, rw, rw],...
            'LineWidth', 2, 'EdgeColor', 'red')
            hold off

            % kill hardware
            hardware.kill()
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STUFF IN 'd' OBJECT
        freq_scan_cw_to_save;
        
        if ~get(handles.rfPowerSweep,'Value')
            disp('saving data')
            % save data
        
            % save data for single RF power       
            delete([data_folder '\' mat_file_name], [data_folder '\' final_fig_name])
            
            save([data_folder '\' mat_file_name], 'd') 

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if mod(i,3)==0 %run autofocus on every third sweep
            micro_autofocus
        %end
    end
    
    % fit and plot fit
    fit_successful = 0;
    
    try   
        counts_exposure_time = kin_time;
        
        [eta, dB, max_contrast, FWHM, counts, final_fit, fit_frequencies] = get_sensitivity(avg_scan, counts_exposure_time, frequencies, photoelectrons, power, lorentzian_type, 0);
        axes(handles.axes1);
        plot(fit_frequencies, final_fit, 'r');
        hold on
        
        fit_successful = 1;
        
        eta
    catch me
        disp(me)
    end     
    
    % average sweep figure with fit

    figure(3); clf;
    plot(frequencies,avg_scan,'.-b')
    hold on
    if fit_successful
        plot(fit_frequencies, final_fit, 'r');
    end
    xlabel('frequency (GHz)');
    ylabel('pl');
    axis tight
    figtitle = ['RF power = ' num2str(power) ' dBm, Sweep time = ' num2str(sweep_time) ' sec, Exposure time = ' num2str(exposure_time) ' s, ' num2str(num_sweeps) ' sweep(s)'];
    title(figtitle,'FontSize',8);
    set(figure(3), 'Name', 'avg_scans');
        
    % plot freq scans
    axes(handles.axes1);
    hold on
    plot(frequencies,avg_scan,'.-b');
    box on
    xlabel('frequency (GHz)')
    ylabel('counts')
    hold off

    
    disp('saving data')
    % save data
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STUFF IN 'd' OBJECT
    freq_scan_cw_to_save;
    if fit_successful
        freq_scans_cw_FIT_to_save;  
    end

     if get(handles.rfPowerSweep,'Value')
        fileName = strcat('minus_', num2str(abs(power)), 'dbm_sweeps.mat');
        savefig(['rfPower_',num2str(power),'_avgSweep.fig'])
        save(fileName,'d','-v7.3') 
    else
        savefig(figure(3),[data_folder '\' 'avg_sweep.fig']) 
        save([data_folder '\' mat_file_name], 'd', '-v7.3')
        
        figname = ['avgSweep_' dateandtime '.fig'];
        savefig(figure(3), [log_folder '\' figname])
    end
        
    laser_off;    
end
