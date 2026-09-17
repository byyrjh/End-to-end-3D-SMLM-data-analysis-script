Fitting protocol workflow for single molecule data
1. To review segmentation run MainUI.m and open one of PSFdata_**_m*.mat files
2. To prepare for fitting data and fitting info run experiment_seg_generate_m*.m
3. To create segment selection indice (This is specifically for single molecule data) run seg_vis_SM.m 
4. run data_refitting.m parm fitting_round = 1 
5. To generate 4D light sheet offset map "fitting_init_SM*.mat" run os_map_gen.m
6. run data_refitting.m again parm fitting_round = 2 


