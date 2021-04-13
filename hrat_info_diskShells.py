import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity
import matplotlib

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

# ts = yt.load("Plot Files/More Plot Files/parkerCRs_hdf5_plt_cnt_0*")
ts = yt.load("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0*")
rho_ex = []
times = []

avgCRPres0 = []
avgCRPres1 = []
avgCRPres2 = []
avgCRPres3 = []
avgCRPres4 = []
avgCRPres5 = []
avgGasPres0 = []
avgGasPres1 = []
avgGasPres2 = []
avgGasPres3 = []
avgGasPres4 = []
avgGasPres5 = []

#def numberDensity(field,data):
#        return data['density']

#yt.add_field(("gas","numberDensity"), function=grav,units="g/cm**2/s**2")

pUnit = YTQuantity(1, 'cm**2/s**2')

def Brho(field,data):
        return np.sqrt(data['magx']**2 + data['magy']**2 + data['magz']**2)/data['dens']

def CRPres(field,data):
            return 0.33333*(data['cray']*data['density'])

def GasPres(field,data):
                return data['pressure'] - 0.33333*(data['cray']*data['density']*pUnit)

yt.add_field(("gas","GasPres"), function=GasPres,units="g/cm/s**2")
yt.add_field(("gas","CRPres"), function=CRPres,units="g/cm**3")
yt.add_field(("gas","Brho"), function=Brho,units="G*cm**3/g")


massISM = []
for ds in ts:
        time = ds.current_time.in_units("Myr")
        times.append(time)
        dd = ds.all_data()
        print(ds.derived_field_list)
        print(dd.quantities.extrema("dx"))

        dsk1=ds.disk("center",[0,1,0],(100.0,"kpc"),(1.0,"kpc"))
        dsk2=ds.disk("center",[0,1,0],(100.0,"kpc"),(2.0,"kpc"))
        dsk3=ds.disk("center",[0,1,0],(100.0,"kpc"),(3.0,"kpc"))
        dsk4=ds.disk("center",[0,1,0],(100.0,"kpc"),(4.0,"kpc"))
        dsk5=ds.disk("center",[0,1,0],(100.0,"kpc"),(20.0,"kpc"))

        dskneg1 = dsk2-dsk1
        dskneg2 = dsk3-dsk2
        dskneg3 = dsk4-dsk3
        dskneg4 = dsk5-dsk4
       # dskneg5 = dsk6-dsk5
       # ad = ds.cut_region(sp,["obj['ism '] > 0.6"])
       # ad1 = ds.cut_region(dskneg1,["obj['ism '] > 0.01"])
       # ad2 = ds.cut_region(dskneg2,["obj['ism '] > 0.01"])
       # ad3 = ds.cut_region(dskneg3,["obj['ism '] > 0.01"])
       # ad4 = ds.cut_region(dskneg4,["obj['ism '] > 0.01"])
       # ad5 = ds.cut_region(dskneg5,["obj['ism '] > 0.01"])
        ad0 = dsk1
        ad1 = dskneg1
        ad2 = dskneg2
        ad3 = dskneg3
        ad4 = dskneg4
       # ad5 = dskneg5
        averageCRPres0 = ad0.quantities.weighted_average_quantity("cloo",weight="ones")
        avgCRPres0.append(abs(averageCRPres0))
        averageCRPres1 = ad1.quantities.weighted_average_quantity("cloo",weight="ones")
        avgCRPres1.append(abs(averageCRPres1))
        averageCRPres2 = ad2.quantities.weighted_average_quantity("cloo",weight="ones")
        avgCRPres2.append(abs(averageCRPres2))
        averageCRPres3 = ad3.quantities.weighted_average_quantity("cloo",weight="ones")
        avgCRPres3.append(abs(averageCRPres3))
        averageCRPres4 = ad4.quantities.weighted_average_quantity("cloo",weight="ones")
        avgCRPres4.append(abs(averageCRPres4))
       # averageCRPres5 = ad5.quantities.weighted_average_quantity("CRPres",weight="ones")
       # avgCRPres5.append(averageCRPres5)

        averageGasPres0 = ad0.quantities.weighted_average_quantity("hrat",weight="ones")
        avgGasPres0.append(averageGasPres0)
        averageGasPres1 = ad1.quantities.weighted_average_quantity("hrat",weight="ones")
        avgGasPres1.append(averageGasPres1)
        averageGasPres2 = ad2.quantities.weighted_average_quantity("hrat",weight="ones")
        avgGasPres2.append(averageGasPres2)
        averageGasPres3 = ad3.quantities.weighted_average_quantity("hrat",weight="ones")
        avgGasPres3.append(averageGasPres3)
        averageGasPres4 = ad4.quantities.weighted_average_quantity("hrat",weight="ones")
        avgGasPres4.append(averageGasPres4)
       # averageGasPres5 = ad5.quantities.weighted_average_quantity("GasPres",weight="ones")
       # avgGasPres5.append(averageGasPres5)
       # print("Total mass in sphere is %0.3e Msun" % (ism_mass.in_units('Msun')))
"""
print(massISM)
plt.semilogy(times,avgCRPres0,'k',label="Cooling + Photoionization Heating, 0 - 1 kpc")
#plt.plot(times,avgGasPres0,'k--',label="CR Heating, 0 - 1 kpc")
plt.semilogy(times,avgCRPres1,'r',label="Cooling + Photoionization Heating, 1 - 2 kpc")
#plt.plot(times,avgGasPres1,'r--',label="CR Heating, 1 - 2 kpc")
plt.semilogy(times,avgCRPres2,'g',label="Cooling + Photoionization Heating, 2 - 5 kpc")
#plt.plot(times,avgGasPres2,'g--',label="CR Heating, 2 - 5 kpc")
plt.semilogy(times,avgCRPres3,'m',label="Cooling + Photoionization Heating, 5 - 10 kpc")
#plt.plot(times,avgGasPres3,'m--',label="CR Heating, 5 - 10 kpc")
plt.semilogy(times,avgCRPres4,'b',label="Cooling + Photoionization Heating, 10 - 20 kpc")
#plt.plot(times,avgGasPres4,'b--',label="CR Heating, 10 - 20 kpc")
plt.xlabel("Time (Myrs)",fontsize=20)
plt.ylabel(r"dE/dt (erg cm$^{-3}$ s$^{-1}$ )",fontsize=20)
plt.ylim(1e-32,1e-27)
plt.legend(loc='lower center',ncol=2,fontsize='x-small')
plt.tight_layout()
plt.savefig("cooling_heating_diskShells.pdf")
plt.close()

"""

print(massISM)
#plt.semilogy(times,avgCRPres0,'k',label="Cooling, height = 0 - 1 kpc")
plt.semilogy(times,avgGasPres0,'k--',label="CR Heating, height = 0 - 1 kpc")
#plt.semilogy(times,avgCRPres1,'r',label="Cooling + Photo Heating, 1 - 2 kpc")
#plt.semilogy(times,avgGasPres1,'r--',label="CR Heating, 1 - 2 kpc")
#plt.semilogy(times,avgCRPres2,'g',label="Cooling, height = 1 - 2 kpc")
plt.semilogy(times,avgGasPres2,'g--',label="CR Heating, height = 1 - 2 kpc")
#plt.semilogy(times,avgCRPres3,'m',label="Cooling, height = 2 - 3 kpc")
plt.semilogy(times,avgGasPres3,'m--',label="CR Heating, height = 2 - 3 kpc")
#plt.semilogy(times,avgCRPres4,'b',label="Cooling + Photo Heating, 10 - 20 kpc")
#plt.semilogy(times,avgGasPres4,'b--',label="CR Heating, 10 - 20 kpc")
plt.xlabel("Time (Myrs)",fontsize=20)
plt.ylabel(r"dE/dt (erg cm$^{-3}$ s$^{-1}$ )",fontsize=20)
plt.ylim(1e-34,5e-26)
plt.legend(loc='upper center',ncol=2,fontsize='small')
plt.tight_layout()
plt.savefig("Plots/CRHeating_heights.pdf")
plt.close()


"""

print(massISM)
#plt.plot(times,avgCRPres0,'k',label="Cooling + Photoionization Heating, 0 - 1 kpc")
plt.semilogy(times,avgGasPres0,'k--',label="CR Heating, 0 - 1 kpc")
#plt.plot(times,avgCRPres1,'r',label="Cooling + Photoionization Heating, 1 - 2 kpc")
#plt.semilogy(times,avgGasPres1,'r--',label="CR Heating, 1 - 2 kpc")
#plt.plot(times,avgCRPres2,'g',label="Cooling + Photoionization Heating, 2 - 5 kpc")
plt.semilogy(times,avgGasPres2,'g--',label="CR Heating, 2 - 5 kpc")
#plt.plot(times,avgCRPres3,'m',label="Cooling + Photoionization Heating, 5 - 10 kpc")
plt.semilogy(times,avgGasPres3,'m--',label="CR Heating, 5 - 10 kpc")
#plt.plot(times,avgCRPres4,'b',label="Cooling + Photoionization Heating, 10 - 20 kpc")
#plt.semilogy(times,avgGasPres4,'b--',label="CR Heating, 10 - 20 kpc")
plt.xlabel("Time (Myrs)",fontsize=20)
plt.ylabel(r"dE/dt (erg cm$^{-3}$ s$^{-1}$ )",fontsize=20)
#plt.ylim(1e-17,5e-13)
plt.ylim(1e-32,1e-27)
plt.legend(loc='lower left',ncol=2,fontsize='x-small')
plt.tight_layout()
plt.savefig("CRHeating_diskShells.pdf")
plt.close()

"""
