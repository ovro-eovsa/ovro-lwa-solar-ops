import numpy as np
import os,logging
from casatools import table


class fast_vis_ant_selector():
    '''
    This goal of this class to find a selection of "good" antennas which
    will be used for getting the fast visibilities. Here "good" antennas 
    just mean that the antenna health is ok to the best of our knowledge.
    The antenna selection does not take into account the performance of the
    PSF. It only finds the nearest available unflagged antenna with which
    it can replace a flagged antenna, which exist in the fast vis antenna 
    list.
    :param bcal: Bandpass table using which we determine the flagged antennas
    :type bcal: str
    :param antfile: This is the configuration file which is used in OVRO-LWA
                    to set the configurations
    :type antfile : str
    :param badantsfile: This is an optional parameter. This file is written
                        by the flag_bad_ants function in the flagging module
                        of the OVRO-LWA solar package. This is essentially a 
                        comma separated list of bad antenna numbers in zero 
                        indexing
    :type badantsfile: str
    :param num_bins: The number of bins into which the X and Y axis will be
                     separated. This is primarily for computational efficiency.
                     The code will search over all antennas if needed, to find
                     a good antenna.
    :param antenna_coords: This is the EW and NS coordinates of the antennas. If
                            not provided, these will be loaded using the 
                            lwa_ewns_coords module
    :type antenna_coords: Optional, a numpy array of shape (352,2/3/4..).
    :param antenna_names: This is the list containing the antenna names. Should be
                            in same format as that provided in lwa_ewns_coords 
                            module. If not provided, will be loaded from the default
                            module
    :type antenna_names: str
    :param figfile: This is the filename to which the antenna plot will be saved
                    by default, if plotting is requested. Default is "ant_coords.png"
    :type figfile: str
    '''
    def __init__(self, bcal, antfile, badantsfile=None, num_bins=5,\
                    antenna_coords=None,antenna_names=None, figfile="ant_coords.png"):
        self.bcal=bcal
        self.antfile=antfile
        self.num_fast_ants=48
        self.antenna_coords=antenna_coords
        self.antenna_names=antenna_names
        self.num_bins=num_bins
        self.badantsfile=badantsfile
        self.figfile=figfile
    
    @property
    def bcal(self):
        return self._bcal
    
    @bcal.setter
    def bcal(self,value):
        if os.path.isdir(value):
            self._bcal=value
        else:
            logging.error("bandpass table does not exist")
            raise IOError("bandpass table does not exist")
            
    @property
    def antfile(self):
        return self._antfile
    
    @antfile.setter
    def antfile(self,value):
        if os.path.isfile(value):
            self._antfile=value
        else:
            logging.error("Antenna file does not exist")
            raise IOError("Antenna file does not exist")
    
    @property
    def antenna_coords(self):
        return self._antenna_coords
    @antenna_coords.setter
    def antenna_coords(self,value):
        if value is None :
            logging.info("Antenna coordinates not provided. "+\
                            "Will load from default")
            import lwa_ewns_coords as coords
            self._antenna_coords=coords.antenna_coords
            
        else:
            self._antenna_coords=value
            
    @property
    def antenna_names(self):
        return self._antenna_names
    @antenna_names.setter
    def antenna_names(self,value):
        if value is None :
            logging.info("Antenna names not provided. "+\
                            "Will load from default")
            import lwa_ewns_coords as coords
            self._antenna_names=coords.antenna_names
            
            
        else:
            self._antenna_names=value
                

    def get_fast_ant_list(self):
        with open(self.antfile,'r') as f1:
            lines=f1.readlines()
        
        for j,line1 in enumerate(lines):
            if 'fast_vis_ants' in line1:
                break
        self.fast_vis_ants=[None]*self.num_fast_ants
        j+=1
        for i in range(self.num_fast_ants):
            self.fast_vis_ants[i]=lines[j].strip()[2:].replace('-','')
            j+=1
            
    def make_antenna_grid(self):
        '''
        This function makes a grid of size (num_bins+1) x (num_bins+1)
        Each cell of the grid contains the antennas which are located
        inside it and their properties. This is helpful for finding the
        nearest antennas, which now we only need to search over a small
        number of antennas.
        '''
        self.xcoords=self.antenna_coords[:,0]
        self.ycoords=self.antenna_coords[:,1]
        
        xcen=np.mean(self.xcoords)
        ycen=np.mean(self.ycoords)
        
        self.xcoords-=xcen
        self.ycoords-=ycen
        
        self.minx=np.min(self.xcoords)
        self.miny=np.min(self.ycoords)
        
        self.grid_sepx=(np.max(self.xcoords)-np.min(self.xcoords))/self.num_bins
        self.grid_sepy=(np.max(self.ycoords)-np.min(self.ycoords))/self.num_bins
        
        self.grid=[]
        for i in range(self.num_bins+1):
            self.grid.append([None]*(self.num_bins+1))
        
        self.num_slow_ants=len(self.antenna_coords)
        
        self.antenna_props={}
        
        for i in range(self.num_slow_ants):
            index2=int((self.xcoords[i]-self.minx)/self.grid_sepx)
            index1=int((self.ycoords[i]-self.miny)/self.grid_sepy)
            if self.grid[index1][index2] is None:
                self.grid[index1][index2]={'num_ant':0,self.antenna_names[i]:{}}
            
            self.grid[index1][index2]['num_ant']+=1
            self.grid[index1][index2][self.antenna_names[i]]={'coords':np.array([self.xcoords[i],self.ycoords[i]])}
            self.grid[index1][index2][self.antenna_names[i]]['flagged']= self.antenna_flags[i]
            self.antenna_props[self.antenna_names[i]]=[index1,index2]
        

    def get_flagged_slow_ants(self):
        '''
        This function determines the antennas which are flagged. For the 
        bandpass table, I assume that an antenna is flagged, if 50% or more 
        channels is flagged. Also I put an antenna as flagged, if it is flagged
        in both polarisations.
        '''
        tb=table()
        tb.open(self.bcal)
        try:
            flag=tb.getcol('FLAG')
        finally:
            tb.close()
        
        flag=np.array(flag,dtype=int)
        flag_chan_avged=np.mean(flag,axis=1)
        flag_chan_avged[flag_chan_avged>0.5] = 1
        flag_chan_avged[flag_chan_avged<=0.5]= 0
        flag_chan_avged=np.array(flag_chan_avged,dtype=int)
        flag_pol_avged=np.array(np.mean(flag_chan_avged,axis=0),dtype=int)
        self.antenna_flags=np.array(flag_pol_avged, dtype=bool)
        if self.badantsfile:
            if os.path.isfile(self.badantsfile):
                badants=np.array(np.genfromtxt(self.badantsfile,delimiter=','),dtype=int)
 #               print (self.antenna_flags[badants])
                self.antenna_flags[badants]=True
 #               print (self.antenna_flags[badants])
            else:
                print ("Badants file provided does not exist. Ignoring")

        
    def get_flagged_fast_ants(self):
        self.flagged_fast_ants=[]
        for ant_name in self.fast_vis_ants:
            index1,index2=self.antenna_props[ant_name]
            if self.grid[index1][index2][ant_name]['flagged']:
                self.flagged_fast_ants.append(ant_name)
            else:
                del self.grid[index1][index2][ant_name]
                
        
    def find_replacement(self,ant_name):
        '''
        This function tries to find the best suited antenna using
        which we can replace the flagged antenna. First the code
        searches only inside the grid cell in which the flagged antenna
        lies. If it does not find any, then it moves to the neighbouring 
        cells and repeats the search, till it finds a suitable one. If
        the search covers the entire grid, and still fails, the code
        prints a message, that the antenna cannot be replaced and exits.
        '''
        index1,index2=self.antenna_props[ant_name]
        flagged_ant_coords=self.grid[index1][index2][ant_name]['coords']
        
        sep=0
        chosen_ant=None
        while chosen_ant is None:
            min_ind1=max(0,index1-sep)
            max_ind1=min(index1+sep,self.num_bins-1)
            min_ind2=max(0,index2-sep)
            max_ind2=min(index2+sep,self.num_bins-1)
            
        
            max_dist=1e8

            for ind1 in range(min_ind1,max_ind1+1):
                for ind2 in range(min_ind2,max_ind2+1):
                    if self.grid[ind1][ind2] is None:
                        continue
                        
                    keys=self.grid[ind1][ind2].keys()
                    
                    for key in keys:
                        if key=='num_ant':
                            continue
                        if key==ant_name or self.grid[ind1][ind2][key]['flagged']:
                            continue
                        current_ant_coords=self.grid[ind1][ind2][key]['coords']
                        dist=np.sqrt((current_ant_coords[0]-flagged_ant_coords[0])**2+\
                                        (current_ant_coords[1]-flagged_ant_coords[1])**2)
                        if dist<max_dist:
                            max_dist=dist
                            chosen_ant=key
                            chosen_ant_ind1=ind1
                            chosen_ant_ind2=ind2
                            
            sep+=1
            if ((min_ind1==0 and max_ind1==self.num_bins-1) and \
                        (min_ind2==0 and max_ind2==self.num_bins-1)) or chosen_ant:
                break
        if chosen_ant:
            del self.grid[chosen_ant_ind1][chosen_ant_ind2][chosen_ant]
        return chosen_ant
            
        
    def replace_flagged_ants(self):
        for ant_name in self.flagged_fast_ants:
            new_ant=self.find_replacement(ant_name)
            if not new_ant:
                print ("Could not replace ",ant_name)
            else:
                print ("Replace "+ant_name+" by "+new_ant)
                index=self.fast_vis_ants.index(ant_name)
                self.fast_vis_ants[index]=new_ant
            
    def make_fast_ant_list(self,do_plot=False):
        self.get_flagged_slow_ants()
#        print("Flagged slow ants found")
        self.get_fast_ant_list()
#        print ("Obtained the fast ant list")
        self.make_antenna_grid()
#        print ("Antenna grid made")
        self.get_flagged_fast_ants()
#        print ("got the flagged fast vis ants")
        if do_plot:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            mpl.use('Agg')
            fig=plt.figure()
            ax=fig.add_subplot(111)
            self.plot_antennas(ax,size=5,color='r')
            
        self.replace_flagged_ants()
#        print ("replaced ants")
        
        
        if do_plot:
            self.plot_antennas(ax,size=3,color='b')
            ax.set_aspect('equal')
            plt.savefig(self.figfile)
            plt.close()

    def plot_antennas(self,ax,size=5,color='r'):
        for ant_name in self.fast_vis_ants:
            index=self.antenna_names.index(ant_name)
            ax.plot(self.xcoords[index],self.ycoords[index],'o',color=color,markersize=size)
        
     
    

if __name__== '__main__':
    bcal='/home/surajit/Downloads/random/20240505_100407_73MHz.bcal'
    antenna_file='lwa_config_calim_std.yaml'
    fv=fast_vis_ant_selector(bcal=bcal,antfile=antenna_file,badantsfile='20240505_100407_73MHz.badants')
    fv.make_fast_ant_list(do_plot=True)
    print (fv.fast_vis_ants)
        
                
        
        
        
        
    
    
            
    
        
    
