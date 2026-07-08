function his_model(micro_speed)

t = 0:10000; % 0-10 s in milli-second
f_single = @(param, x) param(1) + param(2) ./ (1+exp((param(3)-x)*param(4))); % Single-sigmoidal function
f_double = @(param, x) param(1) + param(2) ./ (1+exp((param(3)-x)*param(4))) - param(5) ./ (1+exp((param(6)-x)*param(7))); % Double-sigmoidal function

opts = statset('nlinfit');
opts.RobustWgtFun = 'huber';

n_iter = 300;
rng('default');

syms y(x)
if ~isnan(micro_speed(1))
    n_samples = length(micro_speed);
    ex_min = min(micro_speed);
    ex_max = max(micro_speed);
    norm_speed = (micro_speed - ex_min) ./ (ex_max - ex_min);

    if (ex_min == ex_max)
        ex_param = [ex_min 0 0 0 0 0 0];
    else
        count = 0;
        ex_params = [];
        bic_value = [];
        for i = 1:n_iter
            tic;
            disp(['Count: ' num2str(count)]);
            count = count + 1;
            lb = [0 -1 0 0];
            ub = [1 1 10000 0.5];
            initial_p = [];
            k = length(lb);
            for p = 1:k
                initial_p(p) = unifrnd(lb(p), ub(p));
            end
            [temp, resnorm] = lsqcurvefit(f_single, initial_p, t, norm_speed, lb, ub, opts);
            if (temp(2) > 0)
                ex_params(count, :) = [temp nan nan nan];
            else
                ex_params(count, :) = [temp(1) nan nan nan -1*temp(2) temp(3:end)];
            end
            bic_value(count) = n_samples * log(resnorm / n_samples) + k * log(n_samples);
            toc;
        end

        for i = 1:n_iter
            tic;
            disp(['Count: ' num2str(count)]);
            count = count + 1;
            lb = [0 0 0 0 0 0 0];
            ub = [1 1 10000 0.5 1 10000 0.5];
            initial_p = [];
            k = length(lb);
            for p = 1:k
                initial_p(p) = unifrnd(lb(p), ub(p));
            end
            [ex_params(count, :), resnorm] = lsqcurvefit(f_double, initial_p, t, norm_speed, lb, ub, opts);
            k = length(initial_p);
            bic_value(count) = n_samples * log(resnorm / n_samples) + k * log(n_samples);
            toc;
        end

        [min_val, min_idx] = min(bic_value); 
        best_params = ex_params(min_idx, :);
        best_params(isnan(best_params)) = 0; 

        y(x) = - best_params(5) ./ (1+exp((best_params(6)-x)*best_params(7)));
        ydd = diff(y, x, 2);
        d2d = double(ydd(t)); % Used for determining HIS segment. 
    end
end
