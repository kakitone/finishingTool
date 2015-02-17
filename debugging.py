from multiprocessing import Pool
import time 


def f(x, kk):
    print "in"
    t =0 
    for i in range(110):
        t = ((t + i) %10 * 10 )  /10 + 10
    time.sleep(1)
    
    print x + kk
    
    

    
def cube(x):
    return x**3



def test1():
    t= time.time()
    p = Pool(4)
    kk = 10
    results = []
    for x in range(1,7):
        results.append(p.apply_async(f, args=(x,kk)))
    output = [p.get() for p in results]
    print output
    print time.time() - t

def test2():
    pool = Pool(processes=4)
    results = [pool.apply_async(cube, args=(x,)) for x in range(1,100)]
    print results
    output = [p.get() for p in results]
    print(output)
    
test1()
