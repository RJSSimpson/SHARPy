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
import datetime as dt


class Timer:
    """@brief timer class to be used with 'with' statements to time
    function calls."""    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

class milliTimer:
    """@brief microsecond precision timer class."""
    def __enter__(self):
        self.start = dt.datetime.now()
        return self
    
    def __exit__(self, *args):
        self.end = dt.datetime.now()
        self.interval = (self.end - self.start).microseconds/1.0e6
     
     
        
def isodd(num):
    """@brief returns True if Num is odd"""
    return bool(num & 1)

def iNode2iElem(iNode, NumNodes, NumNodesElem):
        """@brief Work out what element and sub-element node we are in.
        @param iNode Node number.
        @param NumNodes Total number of nodes.
        @param NumNodesElem Number of nodes per element.
        @return iElem Element number.
        @return iiElem Sub-element node number.
        """
        
        # Work out what element we are in.
        if iNode == 0:
            iElem = 0
        elif iNode < NumNodes-1:
            iElem = int(iNode/(NumNodesElem-1))
        elif iNode == NumNodes-1:
            iElem = int((iNode-1)/(NumNodesElem-1))
            
        # Work out what sub-element node we are in.
        # TODO: Is the iiElem numbering from LtR 0,2,1!?
        if NumNodesElem == 2:
            if iNode < NumNodes-1:
                iiElem = 0
            elif iNode == NumNodes-1:
                iiElem = 1
        elif NumNodesElem == 3:
            iiElem = 0
            if iNode == NumNodes-1:
                iiElem = 1 
            elif isodd(iNode):
                iiElem = 2
        
        return iElem, iiElem

if __name__ == '__main__':
    "example: test timer"
    def test():
        for i in range(500):
            for j in range(500):
                print(math.sqrt(float(i*j)))
        
    with Timer() as t:
        test()
    
    print('Request took %.03f sec.' % t.interval)