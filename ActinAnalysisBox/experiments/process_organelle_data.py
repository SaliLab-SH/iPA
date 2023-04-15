# coding = utf-8


'''
output all distance and angle analysis results

By angdi  and bing

'''


from math import nan
import os, sys
import re
from tabnanny import check
from importlib_metadata import method_cache
import seaborn as sns

import numpy as np
from common.parser import arg
from common.preprocess import *
import glob 
from multiprocessing.dummy import Pool as ThreadPool
from common import preprocess, distances_generator
from scipy.spatial.distance import cdist
import torch 
import copy

filedir = os.path.dirname(os.path.abspath(__file__))
os.chdir(filedir)
sys.path.append(os.path.pardir)

class process_organelles():
    def __init__(self, dataid) -> None:
        self.dataid = dataid

        self.vesicle_type = 0
        self.actin_type = 1
        self.microtube_type = 2
        self.pm_up_type = 7
        self.pm_down_type = 8

        self.vesiclefilename = 'vesicle.mrc'
        self.actinfilename ='actin_filled_points.json'
        self.microtubefilename = 'microtube_filled_points.json'
        self.pm_up_filename = 'UpPM.mrc'
        self.pm_down_filename = 'DownPM.mrc'

        self.actinanglefilename = 'actin_to_pm_angles.json'
        self.plotdata_overlap = False

        self.get_datanewpath()
        self.load_plotdataname()

    def load_datatype(self, file_idx):
        'load data type to see what type of organelle is available for calculation'
        df = pd.read_excel(os.path.join(arg.root_dir, arg.parameter_dir, 'Organelle_Sum.xlsx'))
        organelles = df.values
        for line in organelles:
            if file_idx in line:
                organelle_label = line[1:]
        assert len(organelle_label) > 0
        ## ISG	Actin	MT
        return organelle_label
    
    def check_type(self, organelle1, organelle2):
        'check if organelle type is able to make further calculation'
        organelle_type = self.load_datatype(self.dataid)
        if organelle_type[organelle1] == 1 and organelle_type[organelle2] == 1:
            return True
        else: return False

    def check_organelle_type(self,organelletype):
        # ISG	Actin	MT
        return self.check_type(organelletype[0], organelletype[1])

    def get_filelist(self, dir, Filelist):
        newDir = dir
        if os.path.isfile(dir):
            Filelist.append(dir)
            # Filelist.append(os.path.basename(dir))
        elif os.path.isdir(dir):
            for s in os.listdir(dir):
                #if s == "xxx":
                    #continue
                newDir=os.path.join(dir,s)
                self.get_filelist(newDir, Filelist)
        return Filelist


    def check_data_updated(self, organelletypelist):
        '''
        check if data in updated files
        '''

        self.organelletypelist = organelletypelist
        
        filelist = self.get_filelist(os.path.join(arg.root_dir, arg.processed_data_dir),[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 

        self.vesiclefilepath, self.actinfilepath, self.microtubefilepath= str(),str(),str()
        self.actinanglefilepath = str()
        self.pmupfilepath, self.pmdownfilepath = str(), str()

        if 0 in self.organelletypelist: # vesicle 
            for name in filelist:
                if self.vesiclefilename in name:            self.vesiclefilepath = name
        if 1 in self.organelletypelist: # actin 
            for name in filelist:
                if self.actinfilename in name:              self.actinfilepath = name
                elif self.actinanglefilename in name:       self.actinanglefilepath = name
        if 2 in self.organelletypelist: # microtube 
            for name in filelist:
                if self.microtubefilename in name:          self.microtubefilepath = name
        if 7 in self.organelletypelist: # pm_up
            for name in filelist:
                if self.pm_up_filename in name:             self.pmupfilepath = name
        if 8 in self.organelletypelist: # pm_up
            for name in filelist:
                if self.pm_down_filename in name:           self.pmdownfilepath = name


        # check data
        if 0 in self.organelletypelist: # vesicle 
            if not os.path.exists(self.vesiclefilepath):        return False
        if 1 in self.organelletypelist: # actin 
            if not os.path.exists(self.actinfilepath):          return False
            if not os.path.exists(self.actinanglefilepath):     return False
        if 2 in self.organelletypelist: # microtube 
            if not os.path.exists(self.microtubefilepath):      return False

        return True


    def check_data_available(self, organelleidxlist):
        'check whether raw data exists'

        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        filelist = self.get_filelist(rawdata_dir,[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        datapath = '//'.join(filelist[0].split('\\')[:-1])
        filenamelist = preprocess.load_files(datapath)
        
        self.get_datanewpath()

        if self.vesicle_type in organelleidxlist:
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
        if self.actin_type in organelleidxlist:
            if not os.path.exists(filenamelist[2]): sys.exit('Actin data lacked.')
            else: self.actin_rawdata_name = filenamelist[2]
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
        if self.microtube_type in organelleidxlist:
            if not os.path.exists(filenamelist[4]): sys.exit('Microtube data lacked.')
            else: self.microtube_rawdata_name = filenamelist[4]
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
   


    def get_datanewpath(self):
        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        filelist = self.get_filelist(rawdata_dir,[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        datapath = '//'.join(filelist[0].split('\\')[:-1])
        filenamelist = preprocess.load_files(datapath)
        self.condition = filenamelist[0].split('//')[-2]
        self.location = filenamelist[0].split('//')[-3]
        self.datanewpath = f'{arg.root_dir}/{arg.processed_data_dir}/{self.location}/{self.condition}/{self.dataid}'

    def generate_processed_data(self):
        self.check_data_available(self.organelletypelist)

        self.datanewpath = f'{arg.root_dir}/{arg.processed_data_dir}/{self.location}/{self.condition}/{self.dataid}'
        if not os.path.exists(self.datanewpath): os.makedirs(self.datanewpath)
        
        if self.actin_type in self.organelletypelist:           self.preprocess_actin_data()
        if self.vesicle_type in self.organelletypelist:         self.preprocess_vesicle_data()
        if self.microtube_type in self.organelletypelist:       self.preprocess_microtube_data()
        
        self.check_data_updated(self.organelletypelist)

    def preprocess_vesicle_data(self):
        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        vesicleoutputname = f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}'
        vesicle_data = mrcfile.open(self.vesicle_rawdata_name, permissive=True).data
        with mrcfile.new(vesicleoutputname, overwrite=True) as mrc:
            mrc.set_data(vesicle_data )

    def preprocess_actin_data(self):
        self.preprocess_vesicle_data()
        actindict = preprocess.fill_actin_point_gap(self.actin_rawdata_name, )
        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        
        with open(actinoutputpath, 'w') as f:
            json.dump(actindict, f)

        self.preprocess_actin_angle_data()

    def preprocess_actin_angle_data(self):

        # angles are coming from actin rawdata file

        actinangledict = preprocess.get_actin_to_pm_angle(self.actin_rawdata_name, )

        actinangleoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinanglefilename}'

        with open(actinangleoutputpath, 'w') as f:
            json.dump(actinangledict, f)
        print('saved', actinangleoutputpath)

    def preprocess_microtube_data(self):
        self.preprocess_vesicle_data()
        microtubedict = preprocess.fill_microtube_point_gap(self.microtube_rawdata_name,)
        microtubeoutputpath = f'{self.datanewpath}/{self.dataid}_{self.microtubefilename}'
        with open(microtubeoutputpath, 'w') as f:
            json.dump(microtubedict, f)
        print('saved', microtubeoutputpath)



    def load_plotdataname(self):

        self.singledataoutplotpath = f'{arg.root_dir}/{arg.output_dir}/singledata/'
        self.conditiondataoutplotpath = f'{arg.root_dir}/{arg.output_dir}/conditions/'
        self.singledatavis_outplotpath = f'{arg.root_dir}/{arg.output_dir}/singlevisdata/'

        self.actin_surround_actin_vect_map_filename = f'actin_surround_actin_vect_map.csv'
        self.actin_surround_actin_plotname = f'actin_surrounding_actin_probability_density_ditribution.png'

        self.actin_surround_actin_vect_dist_map_filename = f'actin_surround_actin_vect_dist_map.csv'
        self.actin_surround_actin_dist_plotname = f'actin_surrounding_actin_dist_probability_density_ditribution.png'


        self.actin_to_vesicle_filename = f'actin_to_vesicle_distances.csv' 
        self.actin_to_vesicle_plotname = f'actin_distance_to_vesicle_distance_distribution.png' 

        self.actin_to_vesicle_dist_angle_pair_filename = f'actin_to_vesicle_distances_angle_pair.csv' 
        self.actin_to_vesicle_dist_angle_pair_plotname = f'actin_distance_to_vesicle_distance_angle_pair_distribution.png' 


        self.actin_vesicle_vis_filename = f'space_vis_actin_no_vesicle.mrc'
        self.microtube_vesicle_vis_filename = f'space_vis_microtube_no_vesicle.mrc'

        self.microtube_surround_actin_vect_map_filename = f'microtube_surround_actin_vect_map.csv'
        self.microtube_surround_actin_plotname = f'microtube_surrounding_actin_probability_density_ditribution.png'

        self.microtube_surround_actin_vect_dist_map_filename = f'microtube_surround_actin_vect_dist_map.csv'
        self.microtube_surround_actin_dist_plotname = f'microtube_surrounding_actin_dist_probability_density_ditribution.png'

        self.microtube_to_vesicle_filename = f'microtube_to_vesicle_distances.csv'
        self.microtube_to_vesicle_plotname = f'microtube_to_vesicle_distance_distribution.png'

        self.actin_angle_distance_pair_filename = f'actin_angle_distance_pair.csv'
        self.actin_angle_distance_pair_plotname = f'actin_angle_distance_pair_distribution.png'

        self.pm_aggrated_actin_angles_pairedangle_filename = f'pm_to_actin_angles_pairedangle_60.csv'


class check_data(process_organelles):

    def __init__(self, dataid) -> None:
        super().__init__(dataid)
    
    def check_organelle(self):
        # ISG	Actin	MT
        return self.check_type(self.actin_type, self.vesicle_type)


    def step(self):
        # ISG	Actin	MT

        if self.check_type(self.vesicle_type, self.actin_type):
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.match_vesicle_actin()


    def match_vesicle_actin(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)  
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        if len(actin_data) == 0:
            print('no filament')
        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data


            if np.sum(vesicle_mem_mask) == 0:
                print('no vesicle')
            else:
                actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

                filament_coords, mask, voxel_size_xyz = actin_data_new, vesicle_mem_mask, voxel_size_xyz

                filament_coords_modified ={}    
                filamentnames = list(filament_coords.keys())
                for idx in range(len(filamentnames)):
                    curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
                    filament_coords_modified[f'{filamentnames[idx]}'] =[]
                    for i, coord_new in enumerate(curfilamentcoordslst):
                        coord_new = [coord_new[0]/voxel_size_xyz[0], coord_new[1]/voxel_size_xyz[1], coord_new[2]/voxel_size_xyz[2]]
                        filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)


                coords_lst = []
                for idx in range(len(filamentnames)):
                    curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
                    coords_lst.extend(curfilamentcoordslst)

                coords = np.array(coords_lst).reshape(-1,3)

                voxel_coords = coords

                coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
                coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
                coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]



                if not (np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0) or not (np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0) or not (np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0):
                    print('x',np.max(coords_x),np.min(coords_x), mask.shape[0])
                    print('y',np.max(coords_y),np.min(coords_y), mask.shape[1])
                    print('z',np.max(coords_z),np.min(coords_z), mask.shape[2])

                assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
                assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
                assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0



class actin_to_actin(process_organelles):
    def __init__(self,dataid):
        super(actin_to_actin, self).__init__(dataid)


    def check_actintype(self):
        # ISG	Actin	MT

        return self.check_type(self.actin_type, self.actin_type)

    def step(self):

        if not self.check_actintype():
            print(f'{self.dataid} actin_to_actin data not exist.')
        else:
            organelletypelst = [self.actin_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_surround_actin_map_plot(overlap = self.plotdata_overlap)



    def actin_surround_actin_map_plot(self, overlap):
        self.get_datanewpath()
        self.actin_surround_actin_vect_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_map_filename}'
        self.actin_surround_actin_vect_dist_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_dist_map_filename}'
        if not os.path.exists(self.actin_surround_actin_vect_map_filepath) or not os.path.exists(self.actin_surround_actin_vect_dist_map_filepath) or overlap == True:
            self.actin_surround_actin_map()
        

        with open(f'{self.actin_surround_actin_vect_map_filepath}', 'r', encoding='utf-8') as f:
            filament_vect_map_lst = json.load(f)

        filament_vect_map_lst_new = []
        vect_count_n = 0
        for filament in filament_vect_map_lst:
            curcoord = [ np.array(vect) for vect in filament ]
            filament_vect_map_lst_new.append(curcoord)
            vect_count_n += len(curcoord)
        filament_vect_map_lst = filament_vect_map_lst_new

         
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_surround_actin_plotname}')
        plt.close()

        filament_vect_dist_map_lst = pd.read_csv(self.actin_surround_actin_vect_dist_map_filepath)

        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_surround_actin_dist_plotname}')
        plt.show()
        plt.close()



    def actin_surround_actin_map(self):
        '''
        1. shift bias
        2. generate map
        '''
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        with open(f'{actinoutputpath}', 'r') as f:
            actin_data = json.load(f)

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        vesicle_mem_mask = mrcfile.open(f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}', permissive=True).data

        imgsize = vesicle_mem_mask.shape

        actin_vect_map_lst = distances_generator.actin_to_actin_distance(actin_data_new, imgsize, voxel_size_xyz)
        actin_vect_map_filename = self.actin_surround_actin_vect_map_filepath
        actin_vect_map_lst_save = []


        for actin in actin_vect_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in actin ]
            actin_vect_map_lst_save.append(cur_vect)


        with open(self.actin_surround_actin_vect_map_filepath, 'w', encoding='utf-8') as f:
                json.dump(actin_vect_map_lst_save, f)

        self.actin_surround_actin_dist_map()


    def actin_surround_actin_dist_map(self):
        
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        with open(f'{actinoutputpath}', 'r') as f:
            actin_data = json.load(f)

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        vesicle_mem_mask = mrcfile.open(f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}', permissive=True).data

        imgsize = vesicle_mem_mask.shape

        dist = distances_generator.actin_to_actin_distance(actin_data_new, imgsize, voxel_size_xyz) 

        
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_dist_map_filename}', index=False)




class actin_to_vesicle(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	
        return self.check_type(self.actin_type, self.vesicle_type)

    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_vesicle_plot(overlap = self.plotdata_overlap)


    def actin_to_vesicle_plot(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}'
        if not os.path.exists(self.actin_to_vesicle_filepath) or overlap == True:
            self.actin_to_vesicle_distance()


        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        dist = pd.read_csv(self.actin_to_vesicle_filepath)

        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_to_vesicle_plotname}')
        # plt.show()
        plt.close()
        

    def actin_to_vesicle_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}', index=False)


class actin_to_vesicle_dist_angle_pair(process_organelles): 
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	
        return self.check_type(self.actin_type, self.vesicle_type)

    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_vesicle_dist_angle_pair_plot(overlap = self.plotdata_overlap)

    def actin_to_vesicle_dist_angle_pair_plot(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_dist_angle_pair_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_dist_angle_pair_filename}'
        if not os.path.exists(self.actin_to_vesicle_dist_angle_pair_filepath) or overlap == True:
            self.actin_to_vesicle_dist_angle_pair_map()


        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        dist_angle_pairs = pd.read_csv(self.actin_to_vesicle_dist_angle_pair_filepath)
        


    def actin_to_vesicle_dist_angle_pair_map(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        with open(self.actinanglefilepath, 'r') as f:
            actin_angle_data = json.load(f)    

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        dist, angle = distances_generator.coords_to_mem_distance_angle_pair_generator(actin_data_new, actin_angle_data, vesicle_mem_mask, voxel_size_xyz)


        dist_angle_savepd = pd.DataFrame({'Distance':dist,'Angle':angle})
    
        dist_angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_dist_angle_pair_filename}', index=False)


class actin_to_vesicle_vis(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT
        return self.check_type(self.actin_type, self.vesicle_type)

    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle_vis data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.actin_to_vesicle_visimg(overlap = self.plotdata_overlap)

    def actin_to_vesicle_visimg(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}'

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        if len(actin_data) == 0:
            print('no actin')
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            vis_img = np.zeros_like(vesicle_mem_mask)
        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data

            actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

            # dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

            vis_img = distances_generator.filadistance2mem_vis(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

        vis_img3 = vis_img.astype(np.int16)
        outputname = f'{self.singledatavis_outplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_vesicle_vis_filename}'
        with mrcfile.new(outputname, overwrite=True) as mrc:
            mrc.set_data(vis_img3)








class microtube_to_actin(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_microtube_actin_type(self):
        # ISG	Actin	MT	\
        return self.check_type(self.microtube_type, self.actin_type)


    def step(self):

        if not self.check_microtube_actin_type():
            print(f'{self.dataid} microtube_to_actin data not exist.')

        else:
            organelletypelst = [self.microtube_type, self.actin_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.microtube_to_actin_map_plot(overlap = self.plotdata_overlap)


    def microtube_to_actin_map_plot(self, overlap):
        self.get_datanewpath()
        # self.microtube_to_actin_map()

        self.microtube_surround_actin_vect_map_filepath= f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}'
        self.microtube_surround_actin_vect_dist_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_dist_map_filename}'
        # Check the data if need to skip
        if not os.path.exists(self.microtube_surround_actin_vect_map_filepath) or not os.path.exists(self.microtube_surround_actin_vect_dist_map_filepath) or overlap == True:
            self.microtube_to_actin_map()

 

    def microtube_to_actin_map(self):
        # shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(arg, self.dataid)

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)


        filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data) # unit = A
        filament_vect_2_MT_map_lst_save =[]
        
        for filament in filament_vect_2_MT_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in filament ]
            filament_vect_2_MT_map_lst_save.append(cur_vect)

        with open(f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}', 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)

        self.microtube_to_actin_dist_map()


    def microtube_to_actin_dist_map(self):

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)


        filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data)
        filament_vect_2_MT_map_lst_save =[]
        
        for filament in filament_vect_2_MT_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in filament ]
            filament_vect_2_MT_map_lst_save.append(cur_vect)

        with open(f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}', 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)



        dist = distances_generator.actin_to_microtube_shortest_distance(actin_data, MT_data)

        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        filelist = get_filelist(rawdata_dir,[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        datapath = '//'.join(filelist[0].split('\\')[:-1])
        filenamelist = preprocess.load_files(datapath)
        actin_rawdata_name = filenamelist[2]
        actinangledict = preprocess.get_actin_to_pm_angle(actin_rawdata_name, )
        actinangle_lst = list(actinangledict.values())

        dist_savepd = pd.DataFrame({'shortest distance': dist, 'actin_angle': actinangle_lst})
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_microtube_surround_actin_test1128.csv', index=False) 



class microtube_to_vesicle(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    

    def step(self):
        organelletypelst = [self.microtube_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} microtube_to_vesicle data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_vesicle_plot(overlap = self.plotdata_overlap)
        print('Done with microtube_to_vesicle.')


    def microtube_to_vesicle_plot(self, overlap):
        self.microtube_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_to_vesicle_filename}'
        if not os.path.exists(self.microtube_to_vesicle_filepath) or overlap == True:
            self.microtube_to_vesicle_distance()

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        dist2 = pd.read_csv(self.microtube_to_vesicle_filepath)
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_to_vesicle_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()


    def microtube_to_vesicle_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(MT_data_new, vesicle_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(self.microtube_to_vesicle_filepath, index=False)


class microtube_to_vesicle_vis(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT
        return self.check_type(self.microtube_type, self.vesicle_type)

    def step(self):

        organelletypelst = [self.microtube_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} microtube_vesicle_vis data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_vesicle_visimg(overlap = self.plotdata_overlap)


    def microtube_to_vesicle_visimg(self, overlap):
        self.get_datanewpath()
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        if len(MT_data) == 0:
            print('no microtube data')
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            vis_img = np.zeros_like(vesicle_mem_mask)

        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)

            vis_img = distances_generator.filadistance2mem_vis(MT_data_new, vesicle_mem_mask, voxel_size_xyz)

        vis_img3 = vis_img.astype(np.int16)
        outputname = f'{self.singledatavis_outplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_vesicle_vis_filename}'
        with mrcfile.new(outputname, overwrite=True) as mrc:
            mrc.set_data(vis_img3)



class actin_angle_distance_pair(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_angle_distance_pair'

    def step(self):
        organelletypelst = [self.actin_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_angle_distance_pair_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def step2(self):
        organelletypelst = [self.actin_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_angle_distance_pair_plot2(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def actin_angle_distance_pair_plot(self, overlap):
        
        self.actin_angle_distance_pair_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_angle_distance_pair_filename}'
        self.curdatapath = self.actin_angle_distance_pair_filepath
        self.actin2actin_angle_distance_pair_filepath = f'{self.datanewpath}/{self.dataid}_actin2actin_angle_distance_pair_filename.csv'
        self.curdatapath2 = self.actin2actin_angle_distance_pair_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.actin_angle_distance_pair_map()


        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)

        dist3 = pd.read_csv(self.curdatapath)
        filaments_min_d_lst, filaments_angle_lst =  dist3['Distances'], dist3['Angles']
        if dist3.shape[0] == 0:
            print(f'{self.dataid} only have one filament.')
        else:
            plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_angle_distance_pair_plotname}'
            plt.savefig(plotname)
            # plt.show()
            plt.close()

    def actin_angle_distance_pair_plot2(self, overlap):
        
        self.actin2actin_angle_distance_pair_filepath = f'{self.datanewpath}/{self.dataid}_actin2actin_angle_distance_pair_filename.csv'
        self.curdatapath2 = self.actin2actin_angle_distance_pair_filepath
        if not os.path.exists(self.curdatapath2) or overlap == True:
            self.actin_angle_distance_pair_map2()



    def actin_angle_distance_pair_map2(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)
        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        imgsize = vesicle_mem_mask.shape
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)     
        

        filaments_min_d_lst, filaments_angle_lst, filament_angle_byfilament = self.actinangle2dist2(filament_data_new, imgsize, voxel_size_xyz)

        actin2actin_distangle_savepd = pd.DataFrame(data= {'Distances':filaments_min_d_lst,
                                          'Angles': filament_angle_byfilament
                                          })
        actin2actin_distangle_savepd.to_csv(self.curdatapath2, index=False)

    def actin_angle_distance_pair_map(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)
        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        imgsize = vesicle_mem_mask.shape
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)

        with open(self.actinanglefilepath, 'r') as f: # 20230209
            actin_angle_data = json.load(f)        
        
        filament_angle_byfilament = list(actin_angle_data.values())

        filaments_min_d_lst, filaments_angle_lst, = self.actinangle2dist(filament_data_new, imgsize, voxel_size_xyz)
        
        actin2actin_distangle_savepd = pd.DataFrame(data= {'Distances':filaments_min_d_lst,
                                          'Angles': filament_angle_byfilament
                                          })

        distangle_savepd = pd.DataFrame(data= {'Distances':filaments_min_d_lst,
                                          'Angles': filaments_angle_lst ,
                                          'coord1': self.filaments_coord1_lst,
                                            'coord2': self.filaments_coord2_lst,
                                          })
        distangle_savepd.to_csv(self.curdatapath, index=False)
        actin2actin_distangle_savepd.to_csv(self.curdatapath2, index=False)


    def actinangle2dist(self, filament_coord_extend_dict, imgsize, voxel_size_xyz):
        'input: filaments distance dict'
        edge = 0.
        x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
        y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
        z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

        filament_coord_extend_lst = []
        filamentnames = list(filament_coord_extend_dict.keys())

        for idx in range(len(filamentnames)):
            # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
            curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
            filament_coord_extend_lst.append(curfilamentcoordslst)

        filament_coord_lst = filament_coord_extend_lst
        filaments_vect_lst = []  # vect for every point in every filament

        for single_filament_coord in filament_coord_lst:
            filament_vectlst = []
            for jj, point in enumerate(single_filament_coord):
                if jj < 1:
                    tempvect = single_filament_coord[jj+1] - point
                else:
                    tempvect = point- single_filament_coord[jj-1]
                # unit_vect = vect / np.linalg.norm(vect)
                filament_vectlst.append(tempvect)  
            filaments_vect_lst.append(filament_vectlst)  # not unit vect 


        filaments_min_d_lst = []
        filaments_angle_lst = []
        filaments_coord1_lst = []
        filaments_coord2_lst = []


        if len(filament_coord_extend_lst) < 2:
            return  filaments_min_d_lst, filaments_angle_lst
        else:
            for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
                curfilaments_vect = filaments_vect_lst[i]

                rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
                rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)

                _ = rest_filament_coord_extend_lst.pop(i)
                _ = rest_filament_vect_extend_lst.pop(i)

                # curfilaments_coords = [point1, point2, ...]
                filament_dist_lst = []
                filament_angle_lst = []
                filament_coord1_lst = []
                filament_coord2_lst = []

                for ll, point in enumerate(curfilaments_coords):
                    np.array(curfilaments_coords).shape 
                    len (curfilaments_coords)
                    vect1 = curfilaments_vect[ll]
                    min_dist = np.inf
                    cur_angle = np.inf
                    for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                        distss = cdist([point], filaments_coords, metric='euclidean')
                        min_d = np.min(distss)
                        min_d_idx = np.where(distss == min_d)[1][0]
                        coord1 = point

                        if min_d < min_dist:
                            '''
                            update min_dist, cur_angle
                            '''
                            temp_d = min_d
                            vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                            coord2 = rest_filament_coord_extend_lst[kk][min_d_idx]

                            if (vect1 == vect2).all():
                                continue
                            else:
                                vector_dot_product = np.dot(vect1, vect2)  
                                arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                                angle = np.degrees(arccos)
                                if angle >= 180:
                                    temp_angle = 0
                                elif angle > 90:
                                    temp_angle = 180 - angle
                                else:
                                    temp_angle = angle
                                    if angle >0 and angle <=90:
                                        pass

                                if temp_angle>=0 or temp_angle<=180:  
                                    min_dist = temp_d
                                    cur_angle = temp_angle
         

                    filament_dist_lst.append(min_dist)
                    filament_angle_lst.append(cur_angle)
                    filament_coord1_lst.append(coord1)
                    filament_coord2_lst.append(coord2)

                filaments_min_d_lst.extend(filament_dist_lst)
                filaments_angle_lst.extend(filament_angle_lst)
                filaments_coord1_lst.extend(filament_coord1_lst)
                filaments_coord2_lst.extend(filament_coord2_lst)

                self.filaments_coord1_lst = filaments_coord1_lst
                self.filaments_coord2_lst = filaments_coord2_lst

        return  filaments_min_d_lst, filaments_angle_lst

    def actinangle2dist2(self, filament_coord_extend_dict, imgsize, voxel_size_xyz):
        'input: filaments distance dict'
        edge = 0.
        x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
        y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
        z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

        filament_coord_extend_lst = []
        filamentnames = list(filament_coord_extend_dict.keys())

        for idx in range(len(filamentnames)):
            # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
            curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
            filament_coord_extend_lst.append(curfilamentcoordslst)

        filament_coord_lst = filament_coord_extend_lst
        filaments_vect_lst = []  # vect for every point in every filament

        for single_filament_coord in filament_coord_lst:
            filament_vectlst = []
            for jj, point in enumerate(single_filament_coord):
                if jj < 1:
                    tempvect = single_filament_coord[jj+1] - point
                else:
                    tempvect = point- single_filament_coord[jj-1]
                # unit_vect = vect / np.linalg.norm(vect)
                filament_vectlst.append(tempvect)  
            filaments_vect_lst.append(filament_vectlst)  # not unit vect 

        filament_coords_modified ={}    
        
        for idx in range(len(filamentnames)):
            curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
            filament_coords_modified[f'{filamentnames[idx]}'] =[]
            for i, coord_new in enumerate(curfilamentcoordslst):
                coord_new = [coord_new[0]/voxel_size_xyz[0], coord_new[1]/voxel_size_xyz[1], coord_new[2]/voxel_size_xyz[2]]
                filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)  ## coord in axis


        filament_coords_byfilament = list(filament_coords_modified.values())

        actin_actin_min_angle_lst =[]
        actin_actin_actinidx = []

        


        filaments_min_d_lst = []
        filaments_angle_lst = []
        filaments_coord1_lst = []
        filaments_coord2_lst = []



        if len(filament_coord_extend_lst) < 2:
            return  filaments_min_d_lst, filaments_angle_lst
        else:
            for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
                curfilaments_vect = filaments_vect_lst[i]

                rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
                rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)

                _ = rest_filament_coord_extend_lst.pop(i)
                _ = rest_filament_vect_extend_lst.pop(i)

                actin_actin_angle_lst =[]

                filament_dist_lst = []
                filament_angle_lst = []
                filament_coord1_lst = []
                filament_coord2_lst = []

                for ll, point in enumerate(curfilaments_coords):
                    np.array(curfilaments_coords).shape
                    len (curfilaments_coords)
                    vect1 = curfilaments_vect[ll]
                    min_dist = np.inf
                    cur_angle = np.inf
                    for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                        distss = cdist([point], filaments_coords, metric='euclidean')
                        min_d = np.min(distss)
                        min_d_idx = np.where(distss == min_d)[1][0]
                        coord1 = point
                        coord_00 = rest_filament_coord_extend_lst[i][0]
                        coord_01 = rest_filament_coord_extend_lst[i][-1]

                        if min_d < min_dist:
                            '''
                            update min_dist, cur_angle
                            '''
                            temp_d = min_d
                            vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                            coord2 = rest_filament_coord_extend_lst[kk][min_d_idx]
                            coord_10 = rest_filament_coord_extend_lst[kk][0]
                            coord_11 = rest_filament_coord_extend_lst[kk][-1]

                            if (vect1 == vect2).all():
                                continue
                            else:
                                vector_dot_product = np.dot(vect1, vect2)  # when nan, 1600
                                arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                                angle = np.degrees(arccos)
                                if angle >= 180:
                                    temp_angle = 0
                                elif angle > 90:
                                    temp_angle = 180 - angle
                                else:
                                    temp_angle = angle
                                    if angle >0 and angle <=90:
                                        pass

                                if temp_angle>=0 or temp_angle<=180:  # exclude nan data
                                    min_dist = temp_d
                                    cur_angle = temp_angle


                    vect1 = np.array(coord_00) - np.array(coord_01)
                    vect2 = np.array(coord_10) - np.array(coord_11)
                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)
                    if angle >= 180:  # include nan
                        cur_angle = 0
                    elif angle > 90:
                        cur_angle = 180 - angle
                    else:
                        cur_angle = angle
                    
                    actin_actin_angle_lst.append(cur_angle)

                    filament_dist_lst.append(min_dist)
                    filament_angle_lst.append(cur_angle)
                    filament_coord1_lst.append(coord1)
                    filament_coord2_lst.append(coord2)

                actin_actin_min_angle_lst.extend(actin_actin_angle_lst)
                filaments_min_d_lst.extend(filament_dist_lst)
                filaments_angle_lst.extend(filament_angle_lst)
                filaments_coord1_lst.extend(filament_coord1_lst)
                filaments_coord2_lst.extend(filament_coord2_lst)

                self.filaments_coord1_lst = filaments_coord1_lst
                self.filaments_coord2_lst = filaments_coord2_lst

        return  filaments_min_d_lst, filaments_angle_lst, actin_actin_min_angle_lst


class pm_aggrated_actin_angles_pairedangle(process_organelles):
    '''
    
    '''

    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'pm_aggrated_actin_angles_pairedangle'
        self.dist_threshold = 60 
        self.actin2actin_threshold = 120 
    def step(self):
        organelletypelst = [self.pm_down_type, self.pm_up_type, self.actin_type]

        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.pm_aggrated_actin_angles_pairedangle_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')



    def pm_aggrated_actin_angles_pairedangle_plot(self, overlap):

        self.focal_aggregated_actin_angles_pair_filepath = f'{self.datanewpath}/{self.dataid}_{self.pm_aggrated_actin_angles_pairedangle_filename}'
        self.curdatapath = self.focal_aggregated_actin_angles_pair_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_aggregated_actin_angles_pair_map()
     
        # self.angle_to_mem_plot()

     
    def focal_aggregated_actin_angles_pair_map(self):
        
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(self.actinanglefilepath, 'r') as f:
            actin_angle_data = json.load(f)        

        pm_down_data = mrcfile.open(self.pmdownfilepath, permissive=True).data
        pm_up_data = mrcfile.open(self.pmupfilepath, permissive=True).data


        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)
        
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)     
        
        
        anglelst1_up, anglelst2_up, anglelst1_down, anglelst2_down = distances_generator.pm_to_actin_angle_angle_pair_generator(filament_data_new, actin_angle_data, pm_up_data, pm_down_data, voxel_size_xyz, self.dist_threshold, self.actin2actin_threshold)
        num = distances_generator.pm_to_actin_number(filament_data_new, actin_angle_data, pm_up_data, pm_down_data, voxel_size_xyz, self.dist_threshold, self.actin2actin_threshold)
        
        angle_up_savepd = pd.DataFrame({'PM2ActinAngle':anglelst1_up, 'Actin2ActinAngle':anglelst2_up})
        angle_down_savepd = pd.DataFrame({'PM2ActinAngle':anglelst1_down, 'Actin2ActinAngle':anglelst2_down})
       
        np.savetxt(fname = f'{self.datanewpath}/{self.dataid}_pm_to_actin_number.csv', X=(num,), delimiter=',', encoding='utf-8')

        angle_up_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_up_{self.pm_aggrated_actin_angles_pairedangle_filename}', index=False)
        angle_down_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_down_{self.pm_aggrated_actin_angles_pairedangle_filename}', index=False)
        
    def actinangle2dist2(self, filament_coord_extend_dict, imgsize, voxel_size_xyz):
        'input: filaments distance dict'
        edge = 0.
        x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
        y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
        z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

        filament_coord_extend_lst = []
        filamentnames = list(filament_coord_extend_dict.keys())

        for idx in range(len(filamentnames)):
            # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
            curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
            filament_coord_extend_lst.append(curfilamentcoordslst)

        filament_coord_lst = filament_coord_extend_lst
        filaments_vect_lst = []  # vect for every point in every filament

        for single_filament_coord in filament_coord_lst:
            filament_vectlst = []
            for jj, point in enumerate(single_filament_coord):
                if jj < 1:
                    tempvect = single_filament_coord[jj+1] - point
                else:
                    tempvect = point- single_filament_coord[jj-1]
                # unit_vect = vect / np.linalg.norm(vect)
                filament_vectlst.append(tempvect)  
            filaments_vect_lst.append(filament_vectlst)  # not unit vect 

        filament_coords_modified ={}    
        
        for idx in range(len(filamentnames)):
            curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
            filament_coords_modified[f'{filamentnames[idx]}'] =[]
            for i, coord_new in enumerate(curfilamentcoordslst):
                coord_new = [coord_new[0]/voxel_size_xyz[0], coord_new[1]/voxel_size_xyz[1], coord_new[2]/voxel_size_xyz[2]]
                filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)  ## coord in axis


        filament_coords_byfilament = list(filament_coords_modified.values())

        actin_actin_min_angle_lst =[]
        filaments_min_d_lst = []
        filaments_angle_lst = []
        filaments_coord1_lst = []
        filaments_coord2_lst = []



        if len(filament_coord_extend_lst) < 2:
            return  filaments_min_d_lst, filaments_angle_lst
        else:
            for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
                curfilaments_vect = filaments_vect_lst[i]

                rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
                rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)
        
                _ = rest_filament_coord_extend_lst.pop(i)
                _ = rest_filament_vect_extend_lst.pop(i)

                # curfilaments_coords = [point1, point2, ...]
                actin_actin_angle_lst =[]

                filament_dist_lst = []
                filament_angle_lst = []
                filament_coord1_lst = []
                filament_coord2_lst = []

                for ll, point in enumerate(curfilaments_coords):
                    np.array(curfilaments_coords).shape 
                    len (curfilaments_coords)
                    vect1 = curfilaments_vect[ll]
                    min_dist = np.inf
                    cur_angle = np.inf
                    for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                        distss = cdist([point], filaments_coords, metric='euclidean')
                        min_d = np.min(distss)
                        min_d_idx = np.where(distss == min_d)[1][0]
                        coord1 = point
                        coord_00 = rest_filament_coord_extend_lst[i][0]
                        coord_01 = rest_filament_coord_extend_lst[i][-1]

                        if min_d < min_dist:
                            '''
                            update min_dist, cur_angle
                            '''
                            temp_d = min_d
                            vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                            coord2 = rest_filament_coord_extend_lst[kk][min_d_idx]
                            coord_10 = rest_filament_coord_extend_lst[kk][0]
                            coord_11 = rest_filament_coord_extend_lst[kk][-1]

                            if (vect1 == vect2).all():
                                continue
                            else:
                                vector_dot_product = np.dot(vect1, vect2)  # when nan, 1600
                                arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                                angle = np.degrees(arccos)
                                if angle >= 180:
                                    temp_angle = 0
                                elif angle > 90:
                                    temp_angle = 180 - angle
                                else:
                                    temp_angle = angle
                                    if angle >0 and angle <=90:
                                        pass
    

                                if temp_angle>=0 or temp_angle<=180:  
                                    min_dist = temp_d
                                    cur_angle = temp_angle


                    vect1 = np.array(coord_00) - np.array(coord_01)
                    vect2 = np.array(coord_10) - np.array(coord_11)
                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)
                    if angle >= 180:  # include nan
                        cur_angle = 0
                    elif angle > 90:
                        cur_angle = 180 - angle
                    else:
                        cur_angle = angle
                    
                    actin_actin_angle_lst.append(cur_angle)

                    filament_dist_lst.append(min_dist)
                    filament_angle_lst.append(cur_angle)
                    filament_coord1_lst.append(coord1)
                    filament_coord2_lst.append(coord2)

                actin_actin_min_angle_lst.extend(actin_actin_angle_lst)
                filaments_min_d_lst.extend(filament_dist_lst)
                filaments_angle_lst.extend(filament_angle_lst)
                filaments_coord1_lst.extend(filament_coord1_lst)
                filaments_coord2_lst.extend(filament_coord2_lst)

                self.filaments_coord1_lst = filaments_coord1_lst
                self.filaments_coord2_lst = filaments_coord2_lst

        return  filaments_min_d_lst, filaments_angle_lst, actin_actin_min_angle_lst


def main():
    '''
    get dataid
    get catagories for calculation
    get list if catagory available
    check if filtered data exist
    generate figure 
    '''
    dataidlst = preprocess.read_dataid()

    for dataid in dataidlst: # check all data 
        aaa = check_data(dataid)
        aaa.step()

    for dataid in dataidlst:
        a2a = actin_to_actin(dataid)
        a2a.step()

        a2v = actin_to_vesicle(dataid) 
        a2v.step()

        a2v_d2a = actin_to_vesicle_dist_angle_pair(dataid) 
        a2v_d2a.step()

        mt2a = microtube_to_actin(dataid)
        mt2a.step()

        mt2i = microtube_to_vesicle(dataid) 
        mt2i.step()

        acta2dp = actin_angle_distance_pair(dataid) 
        acta2dp.step() 
        acta2dp.step2() 

        pm2actin = pm_aggrated_actin_angles_pairedangle(dataid)
        pm2actin.plotdata_overlap = True
        pm2actin.step()

    print('Done.')

if __name__ == '__main__':
    main()
 

