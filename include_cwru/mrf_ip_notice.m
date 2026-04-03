function mrf_ip_notice()

% EULA UI output
lines_win = { ...
'MAGNETIC RESONANCE FINGERPRINTING (MRF) - IP NOTICE'
''
'Magnetic Resonance Fingerprinting (MRF) technology implemented in this software'
'is protected intellectual property owned by Case Western Reserve University'
'(CWRU). The underlying methods and related technology are subject to patent'
'protection.'
''
'MRF is protected by the following foundational patents'
'  - United States: 8,723,518'
'  - United States: 10,241,174'
'  - United States: 10,416,259'
'  - United States: 10,627,468'
'  - Europe: EP 2897523'
'  - Japan: 6557710'
'  - South Korea: 10-1674848'
'  - China: ZL2013800592667'
''
'The foundational scientific publication describing this technology is:'
'  Ma D, Gulani V, Seiberlich N, Liu K, Sunshine JL, Duerk JL, Griswold MA.'
'  Magnetic resonance fingerprinting. Nature. 2013 Mar 14;495(7440):187-192.'
'  doi: 10.1038/nature11971'
''
'Use of the MRF-specific source code and associated software components in this'
'repository is governed by a separate End User License Agreement (EULA).'
'The MRF technology remains the property of Case Western Reserve University'
'and is provided subject to the terms and conditions defined in the applicable EULA.'
''
'By accessing or using the MRF functionality, the user acknowledges and agrees'
'to comply with the terms described in the attached End User License Agreement'
'located in include_cwru/EULA.pdf.'
};

% check if EULA was acknowledged
eula_path = which('EULA.pdf');
eula_path = [eula_path(1:end-8) 'EULA_acknowledged.txt'];
if exist(eula_path, 'file') > 0
    mrf_ip_acknowledged = true;
else
    mrf_ip_acknowledged = false;
end

% ask user to acknowledge EULA.pdf
if ~mrf_ip_acknowledged
    fig = uifigure( ...
        'Name', 'MRF IP Notice', ...
        'Position', [100 100 800 700], ...
        'WindowStyle', 'modal');

    fig.UserData = false;
    fig.CloseRequestFcn = @(src,event) uiresume(src);

    uitextarea(fig, ...
        'Value', lines_win, ...
        'Editable', 'off', ...
        'Position', [20 70 760 610], ...
        'FontSize', 14);

    uibutton(fig, ...
        'Text', 'I Accept', ...
        'Position', [540 20 110 35], ...
        'ButtonPushedFcn', @(src,event) set(fig, 'UserData', true));

    uibutton(fig, ...
        'Text', 'Cancel', ...
        'Position', [670 20 110 35], ...
        'ButtonPushedFcn', @(src,event) uiresume(fig));

    addlistener(fig, 'UserData', 'PostSet', @(src,event) uiresume(fig));

    uiwait(fig);
    accepted = fig.UserData;
    delete(fig);

    if accepted
        fid = fopen(eula_path, 'w');               
        fprintf(fid, sprintf('EULA.pdf was acknowledged by user at %s', datetime()));        
        fclose(fid);
    else
        error(['MRF functionality cannot proceed because the IP notice was not accepted. ', ...
               'By accessing or using the MRF functionality, the user must acknowledge ', ...
               'and agree to comply with the terms described in the End User License ', ...
               'Agreement located in include_cwru/EULA.pdf.']);
    end
end

end