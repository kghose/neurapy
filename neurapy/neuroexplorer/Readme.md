Usage
-----
    from neurapy.neuroexplorer import nexio
    dt = nexio.read_nex('seqdms-jeff-07-21-2012.nex')

This will load all the events, spike time stamps, markers and so on from the nex file

    from neurapy.neuroexplorer import nexio
    dt = nexio.read_nex('seqdms-jeff-07-21-2012.nex', load=['neurons','markers'])

This will load only the neurons and markers from the file. For a list of allowed load strings see the inline documentation e.g. `nexio.read_nex?` in ipython

