import matplotlib.pyplot as plt
import yt
from yt.units.yt_array import YTQuantity

pUnit = YTQuantity(1, 'cm**2/s**2')


def CRPres(field,data):
            return 0.33333*(data['cray']*data['density'])

def GasPres(field,data):
                return data['pressure'] - 0.33333*(data['cray']*data['density']*pUnit)
yt.add_field(("gas","GasPres"), function=GasPres,units="g/cm/s**2")

yt.add_field(("gas","CRPres"), function=CRPres,units="g/cm**3")



# Load the dataset.
ts = yt.DatasetSeries("Plot Files/More Plot Files/parkerCRs_hdf5_plt_cnt_0*")

counter = 0
for ds in ts:

    print(ds.derived_field_list)
    # Create a sphere of radius 1 Mpc centered on the max density location.
    sp = ds.sphere("max", (1, "Mpc"))


    time = str(ds.current_time.in_units('Myr'))
    time = (time[:4]) if len(time) > 4 else time
    t = "{} Myrs".format(str(time))

    # Calculate and store the bulk velocity for the sphere.
    bulk_velocity = sp.quantities.bulk_velocity()
    sp.set_field_parameter('bulk_velocity', bulk_velocity)

    # Create a 1D profile object for profiles over radius
    # and add a velocity profile.
    prof = yt.create_profile(sp, 'y', ('flash', 'hrat'),
                             units = {'y': 'kpc'},
                             extrema = {'y': ((0.0, 'kpc'), (8.0, 'kpc'))},
                             logs = {'y': False},
                             weight_field='cell_mass')

    # Create arrays to plot.
    radius = prof.x
    mean = prof['flash', 'hrat']
    std = prof.standard_deviation['flash', 'hrat']

    # Plot the average velocity magnitude.
    plt.semilogy(radius, mean, label=t)
    # Plot the variance of the velocity magnitude.
   # plt.semilogy(radius, std, label='Standard Deviation')
    plt.xlabel('y (kpc)',fontsize=22)
    plt.ylim(1.E-32,5.E-28)
    plt.ylabel(r'CR Energy Loss Rate (erg $cm^{-3}$ $s^{-1}$',fontsize=18)
   # plt.title("t = 800 Myrs")
    plt.legend()
    plt.tight_layout()
    stri = str(counter).zfill(4)

    name = "hrat_profiles_{}.png".format(stri)

    plt.savefig(name)
    plt.close()

    counter = counter + 1
