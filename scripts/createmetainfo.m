
metaTB.patienttb = {'patient5', '132540-10N', 'normal'; ...
    'patient5', '132540-1T', 'primary'; ...
    'patient5', '137064-1T', 'meta1'; ...
    'patient5', '139508-1T', 'meta2'; ...
    'patient10', '138381-4N', 'normal'; ...
    'patient10', '138381-2T', 'meta1'; ...
    'patient4', '130342-7N', 'normal'; ...
    'patient4', '130342-1T', 'primary'; ...
    'patient4', '130342-13T', 'meta1'; ...
    'patient4', '130342-18T', 'meta2'; ...
    'patient4', '130342-21T', 'meta3'; ...
    'patient4', '135098-7T', 'meta4'; ...
    'patient1', '104024-1N', 'normal'; ...
    'patient1', '104024-1T', 'primary'; ...
    'patient1', '137109-2T', 'meta1'; ...
    'patient1', '147137-2T', 'meta2'; ...
    'patient2', '114166-5T', 'meta1'; ...
    'patient2', '114166-7T', 'meta2'; ...
    'patient3', '126036-14N', 'normal';
    'patient3', '126036-1T', 'meta1'; ...
    'patient3', '121720-2T', 'primary'; ...
    'patient6', '134376-3N', 'normal'; ...
    'patient6', '134376-1T', 'primary'; ...
    'patient6', '142021-2T', 'meta1'; ...
    'patient6', '142021-6N', 'normal'; ...
    'patient7', '144424-20N', 'normal'; ...
    'patient7', '144424-10T', 'primary'; ...
    'patient8', '132054-6N', 'normal'; ...
    'patient8', '132054-2T', 'primary'; ...
    'patient9', '136661-7N', 'normal'; ...
    'patient9', '136661-3T', 'primary'; ...
    'patient11', '142646-22N', 'normal'; ...
    'patient11', '142646-12T', 'primary'; ...
    'patient11', '142646-2T', 'meta1'; ...
    'patient11', '147650-1T', 'meta2'; ...
    'patient12', '146038-4N', 'normal'; ...
    'patient12', '146038-2T', 'primary'; ...
    'patient13', '147805-14N', 'normal'; ...
    'patient13', '147805-2T', 'primary'; ...
    'patient14', '107740-2T', 'primary'};

metaTB.normalTumorPair = {'patient5', '132540-10N', '132540-1T', 'primary'; ...
    'patient5', '132540-10N', '137064-1T', 'meta1'; ...
    'patient5', '132540-10N', '139508-1T', 'meta2'; ...
    'patient10', '138381-4N', '138381-2T', 'meta1'; ...
    'patient4', '130342-7N', '130342-1T', 'primary'; ...
    'patient4', '130342-7N', '130342-13T', 'meta1'; ...
    'patient4', '130342-7N', '130342-18T', 'meta2'; ...
    'patient4', '130342-7N', '130342-21T', 'meta3'; ...
    'patient4', '130342-7N', '135098-7T', 'meta4'; ...
    
    'patient1', '104024-1N', '104024-1T', 'primary'; ...
    'patient1', '104024-1N', '137109-2T', 'meta1'; ...
    'patient1', '104024-1N', '147137-2T', 'meta2'; ...
        
    'patient3', '126036-14N', '126036-1T', 'meta1'; ...
    'patient3', '126036-14N', '121720-2T', 'primary'; ...
        
    'patient6', '134376-3N', '134376-1T', 'primary'; ...
    'patient6', '142021-6N', '142021-2T', 'meta1'; ...    
    
    'patient7', '144424-20N', '144424-10T', 'primary'; ...
    
    'patient8', '132054-6N', '132054-2T', 'primary'; ...
    'patient9', '136661-7N', '136661-3T', 'primary'; ...
    
    'patient11', '142646-22N', '142646-12T', 'primary'; ...
    'patient11', '142646-22N', '142646-2T', 'meta1'; ...
    'patient11', '142646-22N', '147650-1T', 'meta2'; ...
    
    'patient12', '146038-4N', '146038-2T', 'primary'; ...
    'patient13', '147805-14N', '147805-2T', 'primary'};
    

save data/metaTB.mat metaTB