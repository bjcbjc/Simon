function [locidx, filter, error] = readTumorSomaticCalls(fn, caller)
        
    parameters = {'mutect', [1, 2, 7], {}, '%s %f %s'; ...
        'varscan', [1,2,7,8], {}, '%s %f %s %s'; ...
        'strelka', [1,2, 7], {}, '%s %f %s'; ...
        'virmid', [1,2, 7], {}, '%s %f %s'};

    [~, cidx] = ismember(caller, parameters(:,1));
    if cidx == 0
        fprintf('unknown caller %s\n', caller);
        return
    end
    
    [output, status] = VCFFUNC.extract(fn, ...
        'returncol', parameters{cidx,2}, 'format', ...
        parameters{cidx,3}, 'parseformatstr', parameters{cidx,4});
    if status == 0
        if strcmp(caller, 'varscan')
            rmi = ~strcmp(output{3}, 'PASS');
            output{1}(rmi) = [];
            output{2}(rmi) = [];
            output{3}(rmi) = [];
            output{3} = regexp(output{4}, 'SS=(\d)', 'tokens');
            for i = 1:length(output{3})
                if isempty(output{3}{i})
                    output{3}{i} = '5';
                else
                    output{3}{i} = output{3}{i}{1}{1};
                end
            end
            output(4) = [];
            output{3}( strcmp(output{3}, '2') ) = {'PASS'};
        end
        locidx = gloc2index(numericchrm( output{1}), output{2} );
        valid = ~isnan(locidx);
        locidx = locidx(valid);
        filter = output{3}(valid);
        error = '';
    else
        locidx = [];
        filter = {};
        error = output;
    end
    
end