function source_pac(meg_data, meg_trigger, time_window)

addpath('./software/eeglab/');
eeglab nogui;

addpath('./software/fieldtrip_latest/');
ft_defaults;

[val, index] = findpeaks(meg_trigger.trial{1});
index = index(val == 2);
val = val(val == 2);
disp([val; index]);

if (length(initial_trial:last_trial) ~= length(index))
    error('Number of trials do not match');
end

trial_info = [];
trial_info.index = index;
trial_info.practice_time = meg_data.time{1}(index);
trial_info.trial = initial_trial:last_trial;

phase_range = [4 8];
amp_range = [30 55];

EEG = fieldtrip2eeglab(meg_data);
EEG_p = pop_eegfiltnew(EEG, phase_range(1), phase_range(2));
EEG_a = pop_eegfiltnew(EEG, amp_range(1), amp_range(2));

phase_signal = rad2deg(angle(hilbert(EEG_p.data')))';
amp_signal = abs(hilbert(EEG_a.data'))';


n_bins = 18;
% Bining the phases
step_length = 360/n_bins;
phase_bins = -180:step_length:180;

n_trials = 36;
mi = nan(length(meg_data.label), n_trials);
p_dist = nan(length(meg_data.label), n_bins, n_trials);
fsample = meg_data.fsample;

duration = time_window(2) - time_window(1);
n_samples = fsample * duration + 1;
n_ch = length(meg_data.label);

phase_trials = nan(n_samples, n_ch, n_trials);
amp_trials = nan(n_samples, n_ch, n_trials);

for i_trial = 1:n_trials
    kp_index = find(trial_info.trial == i_trial);

    if ~isempty(kp_index)
        if (trial_info.practice_time(kp_index) + time_window(1) > 0 && trial_info.practice_time(kp_index) + time_window(2) < max(meg_data.time{1}))
            start_sample = dsearchn(meg_data.time{1}', trial_info.practice_time(kp_index) + time_window(1));
            end_sample = dsearchn(meg_data.time{1}', trial_info.practice_time(kp_index) + time_window(2));
            sample_window = start_sample:end_sample;

            phase_trials(:, :, i_trial) = phase_signal(:, sample_window)';
            amp_trials(:, :, i_trial) = amp_signal(:, sample_window)';
        end
    end
end

for i_trial = 1:n_trials
    for i_ch = 1:n_ch
        phase = squeeze(phase_trials(:, i_ch, i_trial));
        amp = squeeze(amp_trials(:, i_ch, i_trial));

        [~, phase_bins_ind] = histc(phase,phase_bins);

        amplitude_bins = nan(n_bins,1);
        for bin = 1:n_bins
            amplitude_bins(bin,1) = mean(amp(phase_bins_ind==bin), 1, 'omitnan');
        end

        P = amplitude_bins./repmat(sum(amplitude_bins),n_bins,1);
        p_dist(i_ch, :, i_trial) = P;
        mi(i_ch, i_trial) = 1+sum(P.*log(P))./log(n_bins);
    end
end

n_surrogates = 200;
rng(100);
if n_surrogates
    for i_trial = 1:n_trials
        for i_ch = 1:n_ch
            mi_surr = nan(n_surrogates,1);
            disp([num2str(i_trial) ' ' num2str(i_ch) ' Computing surrogate data...']);

            phase = squeeze(phase_trials(:, i_ch, i_trial));
            phase_signal = phase(:);
            [~, phase_bins_ind] = histc(phase_signal,phase_bins);

            surrogate = amp_trials(:, i_ch, i_trial);
            randind = randi(n_samples - 2*fsample, n_surrogates, 1) + fsample;

            for s=1:n_surrogates
                surrogate_signal = circshift(surrogate, randind(s));
                surrogate_signal = surrogate_signal(:);
                amplitude_bins_surr = zeros(n_bins,1);
                for bin = 1:n_bins
                    amplitude_bins_surr(bin,1) = mean(surrogate_signal(phase_bins_ind==bin), 1, 'omitnan');
                end
                P_surr = amplitude_bins_surr./repmat(sum(amplitude_bins_surr),n_bins,1);
                mi_surr(s) = 1+sum(P_surr.*log(P_surr))./log(n_bins);
            end
            [mn_mi(i_ch, i_trial), std_mi(i_ch, i_trial)] = normfit(mi_surr);
        end
    end
end

mi_z = (mi - mn_mi) ./ std_mi;


