import numpy as np
import pandas as pd
import collections as cl
import glob
import re
from warnings import warn

from scipy.interpolate import LinearNDInterpolator, griddata

class PolarInterpolator(object):
    """
    A class for interpolation of airfoil drag
    """

    def __init__(self, path:str):
        """initializes PolarInterpolator with Polar from path
        
        :param path: search path for polar files
        :type path: str
        """

        # create emtpy data dictionary
        self.data = cl.defaultdict(dict)
        
        # import polars from path
        self._importpath(path)

    def _importpath(self, path:str):
        """import polar files in path
        
        :param path: search path
        :type path: str
        """

        
        self.interpolators = {}
 
        for path in glob.glob(path+'/*_Re*_M*.txt'):
            
            self._loadfile(path)
            
        self._create_interpolation()

    def _loadfile(self, filename:str):
        """load specific polar file 
        
        :param filename: polar file's filename
        :type filename: str
        :raises Exception: Execption(no data found)
        """

        with open(filename, 'r') as afile:
            
            tmpdata = []
            
            for ii, line in enumerate(afile):
                #print((ii,line))
                if ii>10:
                    linedata = np.fromstring(line.strip(), dtype=float, sep=' ')
                    if len(linedata)>0:
                        tmpdata.append(linedata)
                else:
                    if ii == 2:
                        tmpairfoilname = line[22:].strip()
                    #elif ii == 4 and 'Reynolds number fixed' not in line:
                    #    print('wrong polar type in file {}'.format(filename))
                    #    return False
                    elif ii == 7:
                        reg = re.compile(r'\s*Mach\s*=\s*(\d+.\d+)\s*Re\s*=\s*(\d+.\d+\s?e\s+-?\d+)')
                        try:
                            restr = reg.match(line).groups()[1]
                            tmpreynolds = float(restr.replace(' ',''))
                        except:
                            raise Exception('no reynolds data found')
                                    
            columns = 'alpha CL CD CDp Cm Top_Xtr Bot_Xtr Cpmin Chinge XCp'.split()
            
            tmpdata = pd.DataFrame(np.array(tmpdata), columns=columns)
            
            self.data[tmpairfoilname][tmpreynolds] = tmpdata
            
    def _create_interpolation(self):
        """create interpolation for self.data
        
        Uses LinarNDInterpolator as primary Interpolator and 
        stores points, and values for backup interpolator .
        """

        self.backup_data = {}
        for airfoilname, af_data in self.data.items():
            
            values = []
            points = []
            
            for reynolds, rey_data in af_data.items():
                for ii, row in rey_data.iterrows():
                    values.append(row.CD)
                    points.append((reynolds, row.CL))
                    
            self.interpolators[airfoilname] = \
                LinearNDInterpolator(np.array(points), np.array(values), rescale=True)
            self.backup_data[airfoilname] = (points, values)
                 
    def interpolate(self, airfoilname:str, reynolds:float, c_l:float):
        """Interpolate drag for given input
        
        :param airfoilname: name of airfoil
        :type airfoilname: str
        :param reynolds: reynoldsnumber
        :type reynolds: float
        :param c_l: lift coefficient
        :type c_l: float
        :raises Exception: Airfoil not found
        :return: interpolated drag
        :rtype: float
        """
        
        # check if airfoil exists in data
        if not airfoilname in self.interpolators.keys():
                raise Exception('No Polar Data for {}'.format(airfoilname))

        # use default interpolation to find value
        value = float(self.interpolators[airfoilname](reynolds, c_l))

        # if no value could be found use nearest data point
        if np.isnan(value):
            
            warn('Point outside interpolation range, using nearest data point for {} at Re={} and c_l={}'.format(airfoilname,
                    reynolds, c_l))
            value = float(griddata(*self.backup_data[airfoilname],(reynolds, c_l), method='nearest'))

        
        return value