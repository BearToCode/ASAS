function [mag, phase] = asymp_bode(G, w)
    % ASYMP_BODE computes the asymptotic Bode plot of a transfer function.
    %
    % [mag, phase] = asymp_bode(G, w)
    %
    % Input arguments:
    % G         [TransferFunction] transfer function to analyze
    % w         [1xN]           frequencies to evaluate the Bode plot at
    %
    % Output arguments:
    % mag       [1xN]           magnitude of the transfer function at w
    % phase     [1xN]           phase of the transfer function at w

    params = tf_params(G);

    db = @(x) 20 * log10(x);

    w_z = [1 ./ abs(params.t_i), params.a_n'];
    w_p = [1 ./ abs(params.T_i), params.w_n'];

    z = [params.zeros_real', params.zeros_cc'];
    p = [params.poles_real', params.poles_cc'];

    w_singularities = sort(unique([w_z, w_p]));

    if params.type == 0
        starting_line = @(omega) db(abs(params.gain));
    else
        starting_slope = -20 * params.type; % [dB/decade]
        omega_0 = abs(params.gain) ^ (1 / params.type);
        starting_line = @(omega) starting_slope * (log10(omega / omega_0));
    end

    mag_db = arrayfun(@(omega) mag_at(omega, params, w_singularities, w_z, w_p, starting_line), w);
    phase = arrayfun(@(omega) phase_at(omega, params, w_z, w_p, z, p), w);

    mag = 10 .^ (mag_db / 20);

    function slope = slope_at(omega, params, w_z, w_p)
        slope = -20 * params.type;
        slope = slope + length(find(w_z < omega)) * 20;
        slope = slope - length(find(w_p < omega)) * 20;
    end

    function mag_db = mag_at(omega, params, w_singularities, w_z, w_p, starting_line)
        overtaken_w = w_singularities(w_singularities < omega);

        if isempty(overtaken_w)
            mag_db = starting_line(omega);
            return;
        end

        overtaken_slopes = arrayfun(@(w) slope_at(w, params, w_z, w_p), overtaken_w(2:end));
        overtaken_w_log = log10(overtaken_w);

        w_log_diff = diff(overtaken_w_log);
        delta_mag_db = overtaken_slopes .* w_log_diff;

        starting_mag_db = starting_line(overtaken_w(1));

        current_slope = slope_at(omega, params, w_z, w_p);
        current_w_log = log10(omega);
        current_w_log_diff = current_w_log - overtaken_w_log(end);
        final_delta_mag_db = current_slope * current_w_log_diff;

        mag_db = starting_mag_db + sum(delta_mag_db) + final_delta_mag_db;
    end

    function phase = phase_at(omega, params, w_z, w_p, z, p)
        phase = -180 * (params.gain < 0);
        phase = phase - 90 * params.type;

        phase = phase + length(find(w_z < omega & real(z) < 0)) * 90;
        phase = phase + length(find(w_p < omega & real(p) > 0)) * 90;
        phase = phase - length(find(w_z < omega & real(z) > 0)) * 90;
        phase = phase - length(find(w_p < omega & real(p) < 0)) * 90;
    end

end

function params = tf_params(G)
    % TF_PARAMS find the parameters of a transfer function.
    %
    % params = tf_params(G)
    %
    % Input arguments:
    % G         [TransferFunction] transfer function to analyze
    %
    % Output arguments:
    % params    [Struct]        parameters of the transfer function
    % - zeros_real [1xN]        real zeros
    % - zeros_cc   [1xN]        complex conjugate zeros
    % - poles_real [1xN]        real poles
    % - poles_cc   [1xN]        complex conjugate poles
    % - gain       [1x1]          gain
    % - tc         [1x1]          transfer constant
    % - type       [1x1]          type of the system
    % - z_i        [1xN]          negated real cc zeros
    % - p_i        [1xN]          negated real cc poles
    % - n_i        [1xN]          damping ratios of the cc zeros
    % - a_n        [1xN]          natural frequencies of the cc zeros
    % - w_n        [1xN]          natural frequencies of the cc poles
    % - e_i        [1xN]          damping ratios of the cc poles
    % - t_i        [1xN]          time constants of the cc zeros
    % - T_i        [1xN]          time constants of the cc poles

    [tf_num, tf_den] = tfdata(G, 'v');

    tf_zeros = roots(tf_num);
    tf_poles = roots(tf_den);

    tf_zeros_nonzero = tf_zeros(tf_zeros ~= 0);
    tf_poles_nonzero = tf_poles(tf_poles ~= 0);

    tf_zeros_real = tf_zeros_nonzero(imag(tf_zeros_nonzero) == 0);
    tf_poles_real = tf_poles_nonzero(imag(tf_poles_nonzero) == 0);

    tf_zeros_cc = tf_zeros_nonzero(imag(tf_zeros_nonzero) ~= 0);
    tf_poles_cc = tf_poles_nonzero(imag(tf_poles_nonzero) ~= 0);

    tf_null_zeros = length(tf_zeros(tf_zeros == 0));
    tf_null_poles = length(tf_poles(tf_poles == 0));

    type = tf_null_poles - tf_null_zeros;

    % Remove singularities in zero
    tf_clean_num = polydiv(tf_num, poly(zeros(tf_null_zeros, 1)));
    tf_clean_den = polydiv(tf_den, poly(zeros(tf_null_poles, 1)));

    gain = tf_clean_num(end) / tf_clean_den(end);
    tc = tf_clean_num(find(tf_clean_num, 1)) / tf_clean_den(find(tf_clean_num, 1));

    z_i = -tf_zeros_real;
    p_i = -tf_poles_real;

    a_n = arrayfun(@(zero) sqrt(zero * (real(zero) - 1i * imag(zero))), tf_zeros_cc);
    w_n = arrayfun(@(pole) sqrt(pole * (real(pole) - 1i * imag(pole))), tf_poles_cc);

    n_i = arrayfun(@(idx) - (real(tf_zeros_cc(idx) / a_n(idx))), 1:length(tf_zeros_cc));
    e_i = arrayfun(@(idx) - (real(tf_poles_cc(idx) / w_n(idx))), 1:length(tf_poles_cc));

    t_i = 1 ./ z_i;
    T_i = 1 ./ p_i;

    params.zeros_real = tf_zeros_real;
    params.zeros_cc = tf_zeros_cc;
    params.poles_real = tf_poles_real;
    params.poles_cc = tf_poles_cc;
    params.gain = gain;
    params.tc = tc;
    params.type = type;
    params.z_i = z_i;
    params.p_i = p_i;
    params.a_n = a_n;
    params.w_n = w_n;
    params.n_i = n_i;
    params.e_i = e_i;
    params.t_i = t_i;
    params.T_i = T_i;
end
