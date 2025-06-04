function [fOut, tOut, logSpec] = computeLogSpectrogram(sig, Fs_in, NFFT_val, NOverlap_val, fCut)
    % Helper function for EEG spectrogram calculation with improved error handling
    % Use hamming window for better spectral estimation
    
    % Input validation
    if isempty(sig) || length(sig) < NFFT_val
        error('Signal too short for specified NFFT window size');
    end
    
    % Ensure signal is finite
    if ~all(isfinite(sig))
        warning('Signal contains NaN or Inf values, attempting to interpolate');
        sig = filloutliers(sig, 'linear', 'DataVariables', 1);
        if ~all(isfinite(sig))
            sig(~isfinite(sig)) = 0; % Replace remaining non-finite values with zeros
        end
    end
    
    % Create window
    window_vec = hamming(NFFT_val);
    
    try
        % MATLAB's spectrogram with 'psd' gives Power Spectral Density
        [~, fAll_psd, tOut, P_psd] = spectrogram(double(sig), window_vec, NOverlap_val, NFFT_val, Fs_in, 'psd');
        
        % Check for valid outputs
        if isempty(P_psd) || ~all(isfinite(P_psd(:)))
            warning('Spectrogram computation produced invalid results');
        end
        
    catch ME
        error('Spectrogram computation failed: %s', ME.message);
    end
    
    % Select frequencies up to fCut
    keep_freqs = fAll_psd <= fCut;
    fOut = fAll_psd(keep_freqs);
    P_subset = P_psd(keep_freqs,:);
    
    % More robust handling of PSD values
    % Floor PSD values to prevent log(0) or log(negative) for ASD calculation
    % Smallest ASD ~1e-6 µV/√Hz => smallest PSD ~1e-12 (µV)^2/Hz
    min_psd_val = 1e-12;
    
    % Handle NaN and negative values more carefully
    P_subset(~isfinite(P_subset) | P_subset <= 0) = min_psd_val;
    P_subset(P_subset < min_psd_val) = min_psd_val;
    
    % Compute log10(ASD) with additional safety check
    % ASD = sqrt(PSD); units: µV/√Hz if input sig is in µV.
    asd_values = sqrt(P_subset);
    
    % Final check for finite values before log
    asd_values(~isfinite(asd_values) | asd_values <= 0) = sqrt(min_psd_val);
    
    logSpec = log10(asd_values);
    
    % Final validation
    if ~all(isfinite(logSpec(:)))
        warning('Log spectrogram contains non-finite values');
        logSpec(~isfinite(logSpec)) = log10(sqrt(min_psd_val));
    end
end 