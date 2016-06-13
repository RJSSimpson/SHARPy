python
import sys
sys.path.insert(0, '/home/rob/Software/gdb_printers/Eigen')
from printers import register_eigen_printers
register_eigen_printers (None)
end
dir lib:main
