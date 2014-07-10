import mlab
import numpy as np

if __name__ == '__main__':
    for factor in np.arange(0.1,1.1,0.1):
        for timestep in range(100):
            filename = "results_diameter_%0.1/dem%05d"%(factor,timestep)
            src = mlab.pipeline.open(filename)
            print src.outputs[0].point_data.scalars.to_array().shape
    