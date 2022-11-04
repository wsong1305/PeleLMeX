import numpy as np
from matplotlib.colors import LogNorm
import os 
import fnmatch

import yt
from yt.visualization.base_plot_types import get_multi_plot

def _hrr_norm(field, data):
    if(fuel =='tprf'):
        hrr_ref_laminar = 83023976018.02144

    else:
        hrr_ref_laminar = 61748905808.255035
    print('Using {} data to normalize HRR'.format(fuel))
    return data['HeatRelease']/hrr_ref_laminar

def _c(field, data):
    if(fuel =='tprf'):
        T_max = 2370.7240267599327
    else:
        T_max = 2342.648653573437
    return (data['temp']-343)/(T_max-343)


orient = "horizontal"
plot_size = 20. #size of region to plot in mm

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
# path = 'first_ignition/'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
# list_process=['plt00096']
print("To be processed: {} \n".format(list_process))

for i,name in enumerate(list_process):
        # ds = yt.load("{}/{}".format(path,name))


    ds1 = yt.load("{}/{}".format(path,name))  # load data

    time = float(ds1.current_time)

    # ds1.add_field(("gas", "c"),function=_c,sampling_type="cell")
    # ds1.add_field(("gas", "HRR_Normalized"),function=_hrr_norm,sampling_type="cell")

    # There's a lot in here:
    #   From this we get a containing figure, a list-of-lists of axes into which we
    #   can place plots, and some axes that we'll put colorbars.
    # We feed it:
    #   Number of plots on the x-axis, number of plots on the y-axis, and how we
    #   want our colorbars oriented.  (This governs where they will go, too.
    #   bw is the base-width in inches, but 4 is about right for most cases.
    fig, axes, colorbars = get_multi_plot(4, 1, colorbar=orient, bw=4)

    slc1 = yt.SlicePlot(
        ds1,
        "x",
        center=[0.0,0.0,plot_size/1.e+3/2.],
        fields=["temp",'Y(CH2O)','Y(OH)',"Y(POSF11498)","z_velocity","mixture_fraction"],
    )

    # slc1.annotate_contour('c',clim=(0.72,0.73))
    # slc2.annotate_contour('c',clim=(0.72,0.73))


    slc_frb1 = slc1.data_source.to_frb((plot_size/1.e+3), 512)
    # proj_frb = proj.data_source.to_frb((1.0, "Mpc"), 512)

    temp_axes    = [axes[0][0]]
    hrr_axes     = [axes[0][1]]
    species_axes1 = [axes[0][2]] # 
    species_axes2 = [axes[0][3]] # "z_velocity"


    # for dax, tax, vax1, vax2, vax3 in zip(temp_axes, hrr_axes, species_axes1,species_axes2,species_axes3):
#    for dax, tax, vax1, vax2 in zip(temp_axes, hrr_axes, species_axes1, species_axes2):
#
#        dax.xaxis.set_visible(False)
#        dax.yaxis.set_visible(False)
#
#        tax.xaxis.set_visible(False)
#        tax.yaxis.set_visible(False)
#
#        vax1.xaxis.set_visible(False)
#        vax1.yaxis.set_visible(False)
#
#        vax2.xaxis.set_visible(False)
#        vax2.yaxis.set_visible(False)


    # Converting our Fixed Resolution Buffers to numpy arrays so that matplotlib
    # can render them

    slc_temp1 = np.array(slc_frb1["temp"])
    slc_hrr1 = np.array(slc_frb1["Y(OH)"])
    slc_spec1_1  = np.array(slc_frb1["Y(CH2O)"])
    slc_spec1_2  = np.array(slc_frb1["Y(POSF11498)"])
    slc_vfrac    = np.array(slc_frb1["mixture_fraction"])
    
    extent = -plot_size/2.0, plot_size/2.0, 0.0, plot_size
    plots = [
        temp_axes[0]    .imshow(slc_temp1  , origin="lower",extent=extent),
        hrr_axes[0]     .imshow(slc_hrr1   , origin="lower",extent=extent),
        species_axes1[0].imshow(slc_spec1_1, origin="lower",extent=extent), #norm=LogNorm()),
        species_axes2[0].imshow(slc_spec1_2, origin="lower",extent=extent), #norm=LogNorm()),
        # species_axes3[0].imshow(slc_spec1_3, origin="lower"), #norm=LogNorm()),
        ]
    
    stoich_mixfrac = 0.045
    temp_axes    [0].contour(slc_vfrac, [stoich_mixfrac], colors='w', origin="lower",extent=extent)
    hrr_axes     [0].contour(slc_vfrac, [stoich_mixfrac], colors='w', origin="lower",extent=extent)
    species_axes1[0].contour(slc_vfrac, [stoich_mixfrac], colors='w', origin="lower",extent=extent)
    species_axes2[0].contour(slc_vfrac, [stoich_mixfrac], colors='w', origin="lower",extent=extent)

    # "temp" 
    plots[0].set_cmap("turbo")
    plots[0].set_clim((780, 2400.0))
    # plots[0].set_clim((1000.0, 1400.0))

    # "HRR_Normalized"
    plots[1].set_cmap("dusk")
    plots[1].set_clim((0.0, 2.e-3))

    # # "Y(CH4)"
    plots[2].set_cmap("turbo")
    plots[2].set_clim((0.0, 8.e-3))

    # # "Y(CH2O)"
    plots[3].set_cmap("turbo")
    # plots[3].set_clim((0.0, 3.5e-3))

    # # "z_velocity"
    # plots[4].set_cmap("jet")
    # plots[4].set_clim((0.0, 1.2e-4))



    titles = [
        r"$\mathrm{Temperature}$",
        r"$\mathrm{Y(OH)}$",
        r"$\mathrm{Y(CH2O)}$",
        r"$\mathrm{Y(POSF11498)}$",
        # r"$\mathrm{Y(CH3)}$",
    ]

    for p, cax, t in zip(plots, colorbars, titles):
        cbar = fig.colorbar(p, cax=cax, orientation=orient)
        cbar.set_label(t,fontsize=12)
    
    # axes[0][2].set_title('t = {:1.2f} [ms]'.format(time*1.e+3))
    # And now we're done!
    fig.savefig(f"{ds1}_scalars",dpi=300)
