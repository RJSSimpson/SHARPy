'''@package PyBeam.Utils.Misc
@brief      Miscellaneous functions.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       16/01/2013
@pre        None
@warning    None
'''

import time
import math


class Timer:
    """timer class to be used with 'with' statements to time
    function calls."""    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

if __name__ == '__main__':
    "example: test timer"
    def test():
        for i in range(500):
            for j in range(500):
                print(math.sqrt(float(i*j)))
    
    
    with Timer() as t:
        test()
    
    print('Request took %.03f sec.' % t.interval)