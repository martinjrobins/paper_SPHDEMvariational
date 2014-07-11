from multiprocessing import Pool
import os

def f(x):
    factor = x/10.0
    for viscmode in range(2):
        d = "results_diameter_%0.1f_hmode_%d_viscmode_%d"%(factor,0,viscmode)
        print d
        if not os.path.exists(d):
            os.makedirs(d)
        os.system("./single_particle %0.1f %d %d"%(factor,0,viscmode))
    

if __name__ == '__main__':
    pool = Pool(processes=4) 
    pool.map(f, range(1,24))