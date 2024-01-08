"""
This script was written for CASA 5.1.1

Datasets calibrated (in order of date observed):
SB1: 2016.1.00484.L
     Observed 14 May 2017 and 17 May 2017 (2 execution blocks)
LB1: 2016.1.00484.L
     Observed 24 September 2017 and 04 November 2017 (2 execution blocks)

reducer: J. Huang
"""

""" Starting matter """
import os
execfile('reduction_utils.py')
skip_plots = True	# if True, can run script non-interactively

""" Input for loading data """
prefix  = 'GWLup' 
SB1_path = '/full_path/to_calibrated/msfile.ms'
LB1_path = '/full_path/to_calibrated/msfile.ms'
# Note that if you are downloading data from the archive, your SPW numbering 
# may differ from this script, depending on how you split your data out!
data_params = {'SB1': {'vis' : SB1_path,
                       'name' : 'SB1',
                       'field': 'GW_Lup',
                       'line_spws': np.array([0, 4]), # CO SPWs 
                       'line_freqs': np.array([2.30538e11, 2.30538e11]), 
                      }, 
               'LB1': {'vis' : LB1_path,
                       'name' : 'LB1',
                       'field' : 'GW_Lupi',
                       'line_spws': np.array([3,7]), # CO SPWs
                       'line_freqs': np.array([2.30538e11, 2.30538e11]), 
                      }
               }

""" Check data (various options here; an example) """
if not skip_plots:
    for i in data_params.keys():
        plotms(vis=data_params[i]['vis'], xaxis='channel', yaxis='amplitude', 
               field=data_params[i]['field'], ydatacolumn='data', 
               avgtime='1e8', avgscan=True, avgbaseline=True, iteraxis='spw')

""" Identify 50 km/s-wide region containing CO emission; then flag that and do 
    a spectral average to a pseudo-continuum MS """
for i in data_params.keys():      
    flagchannels_string = get_flagchannels(data_params[i], prefix, 
                                           velocity_range=np.array([-15, 25]))
    avg_cont(data_params[i], prefix, flagchannels=flagchannels_string)

""" Define simple masks and clean scales for imaging """
mask_pa  = 42  	# position angle of mask in degrees
mask_maj = 1.3 	# semimajor axis of mask in arcsec
mask_min = 1.1 	# semiminor axis of mask in arcsec
SB1_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
           ('15h46m44.71s', '-34.30.36.09', mask_maj, mask_min, mask_pa)
LB1_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
           ('15h46m44.71s', '-34.30.36.09', mask_maj, mask_min, mask_pa)
SB_scales = [0, 5, 10, 20]
LB_scales = [0, 5, 30, 100, 200]

if not skip_plots:
    """ Image each dataset individually """
    # images are saved in the format prefix+'_name_initcont_exec#.ms'
    image_each_obs(data_params['SB1'], prefix, mask=SB1_mask, scales=SB_scales, 
                   threshold='0.25mJy', interactive=False)
    image_each_obs(data_params['LB1'], prefix, mask=LB1_mask, scales=LB_scales, 
                   threshold='0.06mJy', interactive=False)

    """ Fit Gaussians to roughly estimate centers, inclinations, PAs """
    fit_gaussian(prefix+'_SB1_initcont_exec0.image', region=SB1_mask)
	#Peak : ICRS 15h46m44.710036s -34d30m36.08839s
    fit_gaussian(prefix+'_SB1_initcont_exec1.image', region=SB1_mask)
	#Peak : ICRS 15h46m44.708978s -34d30m36.12469s
    fit_gaussian(prefix+'_LB1_initcont_exec0.image', 
                 region='circle[[%s, %s], %.1farcsec]' % \
                        ('15h46m44.71s', '-34.30.36.09', 0.2))
		# region modified because of noisy image
	#Peak : ICRS 15h46m44.710940s -34d30m36.08832s
    fit_gaussian(prefix+'_LB1_initcont_exec1.image', region=LB1_mask)
	#Peak : ICRS 15h46m44.708871s -34d30m36.09063s
    	#PA of Gaussian component: 41.42 deg
    	#Inclination of Gaussian component: 38.55 deg

""" The emission centers are slightly misaligned.  So we split out the 
    individual executions, shift the peaks to the phase center, and reassign 
    the phase centers to a common direction. (Not in script: check that these 
    shifts do what you think they should by re-imaging!  They do for us.) """

""" Split out individual MSs for each execution """
split_all_obs(prefix+'_SB1_initcont.ms', prefix+'_SB1_initcont_exec')
split_all_obs(prefix+'_LB1_initcont.ms', prefix+'_LB1_initcont_exec')

""" Define a common direction (peak of 2nd LB EB) """
common_dir = 'J2000 15h46m44.70938s -034.30.36.075805'

""" Shift each MS so emission center is at same phase center """
SB0_shift = prefix+'_SB1_initcont_exec0_shift.ms'
os.system('rm -rf '+SB0_shift+'*')
fixvis(vis=prefix+'_SB1_initcont_exec0.ms', outputvis=SB0_shift, 
       field=data_params['SB1']['field'], 
       phasecenter='ICRS 15h46m44.710036s -34d30m36.08839s')
fixplanets(vis=SB0_shift, field=data_params['SB1']['field'], 
           direction=common_dir)
SB1_shift = prefix+'_SB1_initcont_exec1_shift.ms'
os.system('rm -rf '+SB1_shift+'*')
fixvis(vis=prefix+'_SB1_initcont_exec1.ms', outputvis=SB1_shift, 
       field=data_params['SB1']['field'], 
       phasecenter='ICRS 15h46m44.708978s -34d30m36.12469s')
fixplanets(vis=SB1_shift, field=data_params['SB1']['field'], 
           direction=common_dir)
LB0_shift = prefix+'_LB1_initcont_exec0_shift.ms'
os.system('rm -rf '+LB0_shift+'*')
fixvis(vis=prefix+'_LB1_initcont_exec0.ms', outputvis=LB0_shift, 
       field=data_params['LB1']['field'], 
       phasecenter='ICRS 15h46m44.710940s -34d30m36.08832s')
fixplanets(vis=LB0_shift, field=data_params['LB1']['field'],
           direction=common_dir)
LB1_shift = prefix+'_LB1_initcont_exec1_shift.ms'
os.system('rm -rf '+LB1_shift+'*')
fixvis(vis=prefix+'_LB1_initcont_exec1.ms', outputvis=LB1_shift, 
       field=data_params['LB1']['field'],
       phasecenter='ICRS 15h46m44.708871s -34d30m36.09063s')
fixplanets(vis=LB1_shift, field=data_params['LB1']['field'],
           direction=common_dir)


""" Now that everything is aligned, we inspect the flux calibration. """
if not skip_plots:
    """ Assign rough emission geometry parameters. """
    PA, incl = 41.5, 38.5

    """ Export MS contents into Numpy save files """
    for msfile in [prefix+'_SB1_initcont_exec0_shift.ms', 
                   prefix+'_SB1_initcont_exec1_shift.ms', 
                   prefix+'_LB1_initcont_exec0_shift.ms', 
                   prefix+'_LB1_initcont_exec1_shift.ms']:
        export_MS(msfile)

    """ Plot deprojected visibility profiles for all data together """
    plot_deprojected([prefix+'_SB1_initcont_exec0_shift.vis.npz', 
                      prefix+'_SB1_initcont_exec1_shift.vis.npz', 
                      prefix+'_LB1_initcont_exec0_shift.vis.npz', 
                      prefix+'_LB1_initcont_exec1_shift.vis.npz'],
                     fluxscale=[1.0, 1.0, 1.0, 1.0], PA=PA, incl=incl, 
                     show_err=False)
    # An obvious offset of LB1 and everything else; check how much

    """ Now inspect offsets by comparing against a reference """
    estimate_flux_scale(reference=prefix+'_SB1_initcont_exec0_shift.vis.npz', 
                        comparison=prefix+'_LB1_initcont_exec1_shift.vis.npz', 
                        incl=incl, PA=PA)
    #The ratio of comparison : reference is 0.82095
    #The scaling factor for gencal is 0.906 for your comparison measurement

    """ The 2nd LB execution offset is a problem; this calibrator is found to be
	an issue in other DSHARP datasets. """

""" Correct the flux scales where appropriate. """
rescale_flux(prefix+'_LB1_initcont_exec1_shift.ms', [0.906])



"""
SELF-CAL for short-baseline data
"""

""" Merge the SB executions back into a single MS """
SB_cont_p0 = prefix+'_SB_contp0'
os.system('rm -rf %s*' % SB_cont_p0)
concat(vis=[prefix+'_SB1_initcont_exec0_shift.ms', 
            prefix+'_SB1_initcont_exec1_shift.ms'], concatvis=SB_cont_p0+'.ms', 
            dirtol='0.1arcsec', copypointing=False)

""" Set up a clean mask """
mask_PA  = 42 
mask_maj = 1.3
mask_min = 1.1
mask_ra  = '15h46m44.709s'
mask_dec = '-34.30.36.076'
common_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
              (mask_ra, mask_dec, mask_maj, mask_min, mask_PA)

""" Initial clean """
tclean_wrapper(vis=SB_cont_p0+'.ms', imagename=SB_cont_p0, mask=common_mask, 
               scales=SB_scales, threshold='0.2mJy', savemodel='modelcolumn')

""" Define a noise annulus, measure the peak SNR in map """
noise_annulus = "annulus[[%s, %s],['%.2farcsec', '4.25arcsec']]" % \
                (mask_ra, mask_dec, 1.1*mask_maj) 
estimate_SNR(SB_cont_p0+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_SB_contp0.image
#Beam 0.274 arcsec x 0.229 arcsec (-83.62 deg)
#Flux inside disk mask: 86.60 mJy
#Peak intensity of source: 20.76 mJy/beam
#rms: 6.75e-02 mJy/beam
#Peak SNR: 307.49

""" Self-calibration parameters """
SB_contspws = '0~7' 
SB_refant   = 'DV18' 
SB1_obs0_timerange = '2017/05/13/00~2017/05/15/00'
SB1_obs1_timerange = '2017/05/16/00~2017/05/18/00'
 
""" First round of phase-only self-cal (short baselines only) """
SB_p1 = prefix+'_SB.p1'
os.system('rm -rf '+SB_p1)
gaincal(vis=SB_cont_p0+'.ms', caltable=SB_p1, gaintype='T', spw=SB_contspws, 
        refant=SB_refant, calmode='p', solint='30s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=SB_p1, xaxis='time', yaxis='phase',subplot=441, 
            iteration='antenna', timerange=SB1_obs0_timerange, 
            plotrange=[0,0,-180,180]) 
    plotcal(caltable=SB_p1, xaxis='time', yaxis='phase',subplot=441, 
            iteration='antenna', timerange=SB1_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=SB_cont_p0+'.ms', spw=SB_contspws, gaintable=[SB_p1], 
         interp='linearPD', calwt=True)

""" Split off a corrected MS """
SB_cont_p1 = prefix+'_SB_contp1'
os.system('rm -rf %s*' % SB_cont_p1)
split(vis=SB_cont_p0+'.ms', outputvis=SB_cont_p1+'.ms', datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=SB_cont_p1+'.ms', imagename=SB_cont_p1, mask=common_mask, 
               scales=SB_scales, threshold='0.1mJy', savemodel='modelcolumn')
estimate_SNR(SB_cont_p1+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_SB_contp1.image
#Beam 0.274 arcsec x 0.229 arcsec (-83.62 deg)
#Flux inside disk mask: 88.98 mJy
#Peak intensity of source: 22.10 mJy/beam
#rms: 3.26e-02 mJy/beam
#Peak SNR: 678.15

""" Second round of phase-only self-cal (short baselines only) """
SB_p2 = prefix+'_SB.p2'
os.system('rm -rf '+SB_p2)
gaincal(vis=SB_cont_p1+'.ms', caltable=SB_p2, gaintype='T', spw=SB_contspws, 
        refant=SB_refant, calmode='p', solint='18s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=SB_p2, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=SB1_obs0_timerange, 
            plotrange=[0,0,-180,180])
    plotcal(caltable=SB_p2, xaxis='time', yaxis='phase', subplot=441,
            iteration='antenna', timerange=SB1_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=SB_cont_p1+'.ms', spw=SB_contspws, gaintable=[SB_p2], 
         interp='linearPD', calwt=True)

""" Split off a corrected MS """
SB_cont_p2 = prefix+'_SB_contp2'
os.system('rm -rf %s*' % SB_cont_p2)
split(vis=SB_cont_p1+'.ms', outputvis=SB_cont_p2+'.ms', datacolumn='corrected')

""" Image the results; check the resulting map """
wflank_mask = [common_mask, 'circle[[15h46m44.796s, -34d30m33.59s],0.4arcsec]']
# After self-cal, there is a ~7-sigma source several arcsec to the NW; so we 
# add a flanking mask region to properly encompass it
tclean_wrapper(vis=SB_cont_p2+'.ms', imagename=SB_cont_p2, mask=wflank_mask, 
               scales=SB_scales, threshold='0.06mJy', savemodel='modelcolumn') 
estimate_SNR(SB_cont_p2+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_SB_contp2.image
#Beam 0.273 arcsec x 0.230 arcsec (-83.66 deg)
#Flux inside disk mask: 89.24 mJy
#Peak intensity of source: 22.29 mJy/beam
#rms: 3.23e-02 mJy/beam
#Peak SNR: 689.66

""" Amplitude self-cal (short baselines only) """
SB_ap = prefix+'_SB.ap'
os.system('rm -rf '+SB_ap)
gaincal(vis=SB_cont_p2+'.ms', caltable=SB_ap, gaintype='T', spw=SB_contspws, 
        refant=SB_refant, calmode='ap', solint='inf', minsnr=3.0, 
        minblperant=4, solnorm=False) 

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=SB_ap, xaxis='time', yaxis='amp', subplot=441, 
            iteration='antenna', timerange=SB1_obs0_timerange, 
            plotrange=[0,0,0,2])
    plotcal(caltable=SB_ap, xaxis='time', yaxis='amp', subplot=441, 
            iteration='antenna', timerange=SB1_obs1_timerange, 
            plotrange=[0,0,0,2])

""" Apply the solutions """
applycal(vis=SB_cont_p2+'.ms', spw=SB_contspws, gaintable=[SB_ap], 
         interp='linearPD', calwt=True)

""" Split off a corrected MS """
SB_cont_ap = prefix+'_SB_contap'
os.system('rm -rf %s*' % SB_cont_ap)
split(vis=SB_cont_p2+'.ms', outputvis=SB_cont_ap+'.ms', datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=SB_cont_ap+'.ms', imagename=SB_cont_ap, mask=wflank_mask, 
               scales=SB_scales, threshold='0.06mJy', savemodel='modelcolumn')
estimate_SNR(SB_cont_ap+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_SB_contap.image
#Beam 0.274 arcsec x 0.229 arcsec (-83.46 deg)
#Flux inside disk mask: 88.49 mJy
#Peak intensity of source: 22.27 mJy/beam
#rms: 3.04e-02 mJy/beam
#Peak SNR: 731.34



"""
SELF-CAL for the combined (short-baseline + long-baseline) data
"""

""" Merge the SB+LB executions into a single MS """
combined_cont_p0 = prefix+'_combined_contp0'
os.system('rm -rf %s*' % combined_cont_p0)
concat(vis=[SB_cont_ap+'.ms', prefix+'_LB1_initcont_exec0_shift.ms', 
            prefix+'_LB1_initcont_exec1_shift_rescaled.ms'], 
       concatvis=combined_cont_p0+'.ms', dirtol='0.1arcsec', copypointing=False)

""" Initial clean """
tclean_wrapper(vis=combined_cont_p0+'.ms', imagename=combined_cont_p0, 
               mask=wflank_mask, scales=LB_scales, threshold='0.05mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p0+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp0.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 88.86 mJy
#Peak intensity of source: 1.76 mJy/beam
#rms: 1.52e-02 mJy/beam
#Peak SNR: 116.30

""" Self-calibration parameters """
combined_contspws = '0~15'
combined_refant = 'DV09, DV18'
combined_spwmap = [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12] 
LB1_obs0_timerange = '2017/09/23/00~2017/09/25/00'
LB2_obs1_timerange = '2017/11/03/00~2017/11/05/00'

""" First round of phase-only self-cal (all data) """
combined_p1 = prefix+'_combined.p1'
os.system('rm -rf '+combined_p1)
gaincal(vis=combined_cont_p0+'.ms', caltable=combined_p1, gaintype='T', 
        combine='spw,scan', spw=combined_contspws, refant=combined_refant, 
        calmode='p', solint='900s', minsnr=1.5, minblperant=4) 

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=combined_p1, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB1_obs0_timerange, 
            plotrange=[0,0,-180,180]) 
    plotcal(caltable=combined_p1, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB2_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=combined_cont_p0+'.ms', spw=combined_contspws, 
         spwmap=combined_spwmap, gaintable=[combined_p1], interp='linearPD', 
         calwt=True, applymode='calonly')

""" Split off a corrected MS """
combined_cont_p1 = prefix+'_combined_contp1'
os.system('rm -rf %s*' % combined_cont_p1)
split(vis=combined_cont_p0+'.ms', outputvis=combined_cont_p1+'.ms', 
      datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=combined_cont_p1+'.ms', imagename=combined_cont_p1, 
               mask=wflank_mask, scales=LB_scales, threshold='0.045mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p1+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp1.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 89.12 mJy
#Peak intensity of source: 1.90 mJy/beam
#rms: 1.48e-02 mJy/beam
#Peak SNR: 128.42

""" Second round of phase-only self-cal (all data) """
combined_p2 = prefix+'_combined.p2'
os.system('rm -rf '+combined_p2)
gaincal(vis=combined_cont_p1+'.ms', caltable=combined_p2, gaintype='T', 
        combine='spw,scan', spw=combined_contspws, refant=combined_refant, 
        calmode='p', solint='360s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=combined_p2, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB1_obs0_timerange, 
            plotrange=[0,0,-180,180])
    plotcal(caltable=combined_p2, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB2_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=combined_cont_p1+'.ms', spw=combined_contspws, 
         spwmap=combined_spwmap, gaintable=[combined_p2], interp='linearPD', 
         calwt=True, applymode='calonly')

""" Split off a corrected MS """
combined_cont_p2 = prefix+'_combined_contp2'
os.system('rm -rf %s*' % combined_cont_p2)
split(vis=combined_cont_p1+'.ms', outputvis=combined_cont_p2+'.ms', 
      datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=combined_cont_p2+'.ms', imagename=combined_cont_p2, 
               mask=wflank_mask, scales=LB_scales, threshold='0.045mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p2+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp2.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 89.23 mJy
#Peak intensity of source: 1.95 mJy/beam
#rms: 1.48e-02 mJy/beam
#Peak SNR: 131.94
# ** map quality still improving **

""" Third round of phase-only self-cal (all data) """
combined_p3 = prefix+'_combined.p3'
os.system('rm -rf '+combined_p3)
gaincal(vis=combined_cont_p2+'.ms', caltable=combined_p3, gaintype='T', 
        combine='spw,scan', spw=combined_contspws, refant=combined_refant, 
        calmode='p', solint='180s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=combined_p3, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB1_obs0_timerange, 
            plotrange=[0,0,-180,180])
    plotcal(caltable=combined_p3, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB2_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=combined_cont_p2+'.ms', spw=combined_contspws, 
         spwmap=combined_spwmap, gaintable=[combined_p3], interp='linearPD', 
         calwt=True, applymode='calonly')

""" Split off a corrected MS """
combined_cont_p3 = prefix+'_combined_contp3'
os.system('rm -rf %s*' % combined_cont_p3)
split(vis=combined_cont_p2+'.ms', outputvis=combined_cont_p3+'.ms', 
      datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=combined_cont_p3+'.ms', imagename=combined_cont_p3, 
               mask=wflank_mask, scales=LB_scales, threshold='0.045mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p3+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp3.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 88.71 mJy
#Peak intensity of source: 2.06 mJy/beam
#rms: 1.47e-02 mJy/beam
#Peak SNR: 140.14
# ** map quality still improving **

""" Fourth round of phase-only self-cal (all data) """
combined_p4 = prefix+'_combined.p4'
os.system('rm -rf '+combined_p4)
gaincal(vis=combined_cont_p3+'.ms', caltable=combined_p4, gaintype='T', 
        combine='spw,scan', spw=combined_contspws, refant=combined_refant, 
        calmode='p', solint='60s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=combined_p4, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB1_obs0_timerange, 
            plotrange=[0,0,-180,180])
    plotcal(caltable=combined_p4, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB2_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=combined_cont_p3+'.ms', spw=combined_contspws, 
         spwmap=combined_spwmap, gaintable=[combined_p4], interp='linearPD', 
         calwt=True, applymode='calonly')

""" Split off a corrected MS """
combined_cont_p4 = prefix+'_combined_contp4'
os.system('rm -rf %s*' % combined_cont_p4)
split(vis=combined_cont_p3+'.ms', outputvis=combined_cont_p4+'.ms', 
      datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=combined_cont_p4+'.ms', imagename=combined_cont_p4, 
               mask=wflank_mask, scales=LB_scales, threshold='0.045mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p4+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp4.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 88.77 mJy
#Peak intensity of source: 2.25 mJy/beam
#rms: 1.46e-02 mJy/beam
#Peak SNR: 153.70

""" Fifth round of phase-only self-cal (all data) """
combined_p5 = prefix+'_combined.p5'
os.system('rm -rf '+combined_p5)
gaincal(vis=combined_cont_p4+'.ms', caltable=combined_p5, gaintype='T', 
        combine='spw,scan', spw=combined_contspws, refant=combined_refant, 
        calmode='p', solint='30s', minsnr=1.5, minblperant=4)

if not skip_plots:
    """ Inspect gain tables """
    plotcal(caltable=combined_p5, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB1_obs0_timerange, 
            plotrange=[0,0,-180,180])
    plotcal(caltable=combined_p5, xaxis='time', yaxis='phase', subplot=441, 
            iteration='antenna', timerange=LB2_obs1_timerange, 
            plotrange=[0,0,-180,180])

""" Apply the solutions """
applycal(vis=combined_cont_p4+'.ms', spw=combined_contspws, 
         spwmap=combined_spwmap, gaintable=[combined_p5], interp='linearPD', 
         calwt=True, applymode='calonly')

""" Split off a corrected MS """
combined_cont_p5 = prefix+'_combined_contp5'
os.system('rm -rf %s*' % combined_cont_p5)
split(vis=combined_cont_p4+'.ms', outputvis=combined_cont_p5+'.ms', 
      datacolumn='corrected')

""" Image the results; check the resulting map """
tclean_wrapper(vis=combined_cont_p5+'.ms', imagename=combined_cont_p5, 
               mask=wflank_mask, scales=LB_scales, threshold='0.045mJy', 
               savemodel='modelcolumn')
estimate_SNR(combined_cont_p5+'.image', disk_mask=common_mask, 
             noise_mask=noise_annulus)
#GWLup_combined_contp5.image
#Beam 0.035 arcsec x 0.026 arcsec (-89.67 deg)
#Flux inside disk mask: 88.66 mJy
#Peak intensity of source: 2.36 mJy/beam
#rms: 1.46e-02 mJy/beam
#Peak SNR: 161.19

# Additional phase-only cal on shorter intervals or amp self-cal attempts are 
# not helpful.



"""
Final outputs
"""

""" Save the final MS """
os.system('cp -r '+combined_cont_p5+'.ms '+prefix+'_continuum.ms')
os.system('tar cvzf '+prefix+'_continuum.ms.tgz '+prefix+'_continuum.ms')

""" Make a fiducial continuum image (based on experimentation) """
scales = [0, 20, 50, 100, 200]
tclean_wrapper(vis=combined_cont_p5+'.ms', imagename=prefix+'_continuum',
               mask=wflank_mask, scales=scales, threshold='0.05mJy', robust=0.5,
               uvtaper=['0.035arcsec', '0.015arcsec', '0deg'])
exportfits(prefix+'_continuum.image', prefix+'_continuum.fits')
