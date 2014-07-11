from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt
mlab.options.offscreen = True

if __name__ == '__main__':
    n = 18
    fdragz = np.zeros(n)
    f0z = np.zeros(n)
    factors = np.arange(0.1,1.9,0.1)
    for i,factor in zip(range(n),factors):
        f_drag = np.zeros(3)
        f0 = np.zeros(3)
        for timestep in range(50,100):
            filename = 'results_diameter_%0.1f/dem%05d.vtu'%(factor,timestep)
            src = mlab.pipeline.open(filename)
            f_drag = f_drag + src.outputs[0].point_data.get_array(5).to_array()[0,:]
            f0 = f0 + src.outputs[0].point_data.get_array(4).to_array()[0,:]
        f_drag = f_drag / 50;
        f0 = f0 / 50;
        
        fdragz[i] = f_drag[2]
        f0z[i] = f0[2]
        
    plt.figure()
    plt.hold(True)
    dem_vol = (4.0/3.0)*np.pi*(factors*0.0006/2)**3
    dem_mass = dem_vol*2500.0
    sph_visc = 8.9e-07
    sph_dens = 1000.0
    stokes_drag = np.ones(n)*3.0*np.pi*sph_visc*sph_dens*sph_visc
    plt.plot(factors*0.0006/0.0003,dem_mass*fdragz,label='drag')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z,label='pressure')
    plt.plot(factors*0.0006/0.0003, stokes_drag,label='stokes')
    plt.plot(factors*0.0006/0.0003,dem_mass*f0z+dem_mass*fdragz,label='drag+pressure')

    plt.legend()
    plt.show()
    
        
        
    