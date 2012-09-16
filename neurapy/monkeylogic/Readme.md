`bhv_read`
---------
Reads in .bhv files produced by monkeylogic

    from neurapy.monkeylogic import bhv_read as brd

    bhv = brd.read_bhv(fname='my_bhv_file.bhv')


`bhv.keys()` will let you browse the keys, which are the same as the structs that
MonkeyLogic's matlab functions produce.


`moviemaker`
-----------
Generate movies of subject eye position

`python moviemaker.py -f='/my/bhvfile.bhv' -t=363 -x=1.0`

Do `python moviemaker.py --help` to get commandline options