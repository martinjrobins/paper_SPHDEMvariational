from multiprocessing import Pool
import os

def f(x):
    factor = x/10.0
    os.system("./single_particle %0.1f"%factor)

if __name__ == '__main__':
    pool = Pool(processes=4) 
    pool.map(f, range(1,11))