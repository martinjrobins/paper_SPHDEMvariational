from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt
mlab.options.offscreen = True

if __name__ == '__main__':
    n = 23
    fdragz = np.zeros(n)
    f0z = np.zeros(n)
    fdragz2 = np.zeros(n)
    f0z2 = np.zeros(n)
    factors = np.arange(0.1,2.4,0.1)
    for i,factor in zip(range(n),factors):
        f_drag = np.zeros(3)
        f0 = np.zeros(3)
        f_drag2 = np.zeros(3)
        f02 = np.zeros(3)
        for timestep in range(50,100):
            filename = 'results_diameter_%0.1f_hmode_%d_viscmode_%d/dem%05d.vtu'%(factor,0,0,timestep)
            src = mlab.pipeline.open(filename)
            f_drag = f_drag + src.outputs[0].point_data.get_array(5).to_array()[0,:]
            f0 = f0 + src.outputs[0].point_data.get_array(4).to_array()[0,:]
            filename = 'results_diameter_%0.1f_hmode_%d_viscmode_%d/dem%05d.vtu'%(factor,0,1,timestep)
            src = mlab.pipeline.open(filename)
            f_drag2 = f_drag2 + src.outputs[0].point_data.get_array(5).to_array()[0,:]
            f02 = f02 + src.outputs[0].point_data.get_array(4).to_array()[0,:]
        f_drag = f_drag / 50;
        f0 = f0 / 50;
        f_drag2 = f_drag2 / 50;
        f02 = f02 / 50;
        
        fdragz[i] = f_drag[2]
        f0z[i] = f0[2]
        fdragz2[i] = f_drag2[2]
        f0z2[i] = f02[2]
        
    plt.figure()
    plt.hold(True)
    dem_vol = (4.0/3.0)*np.pi*(factors*0.0006/2.0)**3
    dem_mass = dem_vol*2500.0
    sph_visc = 8.9e-07
    sph_dens = 1000.0
    stokes_drag = np.ones(n)*3.0*np.pi*sph_visc*sph_dens*sph_visc
    plt.plot(factors*0.0006/0.0003,dem_mass*fdragz,label='drag (macro)')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z,label='pressure (macro)')
    plt.plot(factors*0.0006/0.0003,-dem_mass*fdragz2,label='drag (micro)')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z2,label='pressure (micro)')
    plt.plot(factors*0.0006/0.0003,stokes_drag,label='stokes')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z+dem_mass*fdragz,label='drag+pressure (macro)')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z2-dem_mass*fdragz2,label='drag+pressure (micro)')


    plt.legend()
    plt.show()
    
        
        
    