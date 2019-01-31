reinitialize
set orthoscopic, on
run spectrumbar.py
load rmsf.pdb, pdb
bg_color white
spectrum b, blue_white_red
#as cartoon
set cartoon_putty_transform, 0
cartoon putty
unset cartoon_smooth_loops
unset cartoon_flat_sheets

#set_view (\
#    -0.991748571,    0.004799802,   -0.128077656,\
#    -0.124684490,   -0.267539561,    0.955443263,\
#    -0.029680848,    0.963532388,    0.265929431,\
#     0.000017576,    0.000065058, -209.263671875,\
#     1.196592808,   -2.955531120,  -17.951896667,\
#   157.330032349,  261.198211670,   20.000000000 )
png rmsf_putty.png, 1800, 1800, 300, 1

hide (all)
spectrumbar red, white, blue
#set_view (\
#    -0.999200165,    0.000329345,    0.039886702,\
#     0.008402883,   -0.975784063,    0.218567386,\
#     0.038992494,    0.218727157,    0.975003719,\
#     0.000000000,    0.000000000,  -28.356409073,\
#     5.000000000,    0.000000000,    0.000000000,\
#    22.356409073,   34.356407166,  -20.000000000 )
png bar.png, 600, 300, 300, 1

