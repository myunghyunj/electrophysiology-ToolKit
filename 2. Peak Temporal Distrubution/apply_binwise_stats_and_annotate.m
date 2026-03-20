function [pvals, test_used] = apply_binwise_stats_and_annotate( ...
    cent, edges, ...
    stimCountsPerFile, shamCountsPerFile, ...
    stimMeanCounts, shamMeanCounts, ...
    stimSEM, shamSEM, ...
    ax, alpha, isPaired)
% Bin-wise stats with plot stars (unpaired by default).
if nargin < 10 || isempty(alpha),   alpha   = 0.05; end
if nargin < 11 || isempty(isPaired), isPaired = false; end
assert(numel(cent) == numel(edges)-1, 'cent/edges size mismatch.');

Nb        = numel(cent);
pvals     = nan(1, Nb);
test_used = repmat({''}, 1, Nb);

fprintf('\n== Bin-wise STIM vs SHAM (alpha = %.3f) ==\n', alpha);

for b = 1:Nb
    x = stimCountsPerFile(:, b); x = x(~isnan(x));
    y = shamCountsPerFile(:, b); y = y(~isnan(y));
    [nx, px, nx_name] = i_isnormal(x, alpha);
    [ny, py, ny_name] = i_isnormal(y, alpha);

    if ~isPaired
        if nx && ny
            [~, p] = ttest2(x, y, 'Vartype','unequal', 'Alpha',alpha, 'Tail','both');
            used = 'Welch t-test';
        else
            p = ranksum(x, y);
            used = 'Mann-Whitney U';
        end
        fprintf('Bin %2d  %5.2f–%5.2f min : p = %.4g  %s  [norm: STIM %s p=%.3g, SHAM %s p=%.3g]\n', ...
            b, edges(b), edges(b+1), p, i_p2stars(p), nx_name, px, ny_name, py);
    else
        if numel(x) == numel(y) && numel(x) >= 2
            d = x - y; [nd, pd, nd_name] = i_isnormal(d, alpha);
            if nd
                [~, p] = ttest(x, y, 'Alpha',alpha, 'Tail','both'); used = 'paired t-test';
            else
                p = signrank(x, y); used = 'Wilcoxon signed-rank';
            end
            fprintf('Bin %2d  %5.2f–%5.2f min : p = %.4g  %s  [norm: DIFF %s p=%.3g]\n', ...
                b, edges(b), edges(b+1), p, i_p2stars(p), nd_name, pd);
        else
            if nx && ny, [~, p] = ttest2(x,y,'Vartype','unequal','Alpha',alpha,'Tail','both'); used='Welch t-test';
            else, p = ranksum(x,y); used='Mann-Whitney U'; end
            fprintf('Bin %2d  %5.2f–%5.2f min : p = %.4g  %s  [paired disabled; norm: STIM %s p=%.3g, SHAM %s p=%.3g]\n', ...
                b, edges(b), edges(b+1), p, i_p2stars(p), nx_name, px, ny_name, py);
        end
    end
    pvals(b)     = p;
    test_used{b} = used;
end

if ishghandle(ax)
    axes(ax); yyaxis left; hold on;
    yl = ylim; dy = 0.03 * diff(yl);
    for b = 1:Nb
        s = i_p2stars(pvals(b));
        if ~isempty(s)
            yTop = max([stimMeanCounts(b)+stimSEM(b), shamMeanCounts(b)+shamSEM(b), eps]);
            text(cent(b), yTop + dy, s, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom','FontWeight','bold','Clipping','off');
        end
    end
end
end

% ---- helpers ----
function [isnorm, p, name] = i_isnormal(v, alpha)
v = v(:); v = v(~isnan(v)); isnorm=false; p=NaN; name='N/A'; if numel(v) < 3, return; end
if exist('swtest','file')==2, try, [h,p]=swtest(v,alpha); isnorm=(h==0); name='Shapiro–Wilk'; return; end, end
if exist('lillietest','file')==2, [h,p]=lillietest(v,'Alpha',alpha); isnorm=(h==0); name='Lilliefors'; return; end
if exist('jbtest','file')==2, [h,p]=jbtest(v,alpha); isnorm=(h==0); name='Jarque–Bera'; return; end
if exist('adtest','file')==2, [h,p]=adtest(v,'Alpha',alpha); isnorm=(h==0); name='Anderson–Darling'; return; end
mu=mean(v); sg=std(v);
if sg>0 && exist('kstest','file')==2, [h,p]=kstest((v-mu)/sg,'Alpha',alpha); isnorm=(h==0); name='Kolmogorov–Smirnov'; return; end
end

function s = i_p2stars(p)
if isnan(p) || isinf(p) || p >= 0.05, s = '';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
else,             s = '*';
end
end
