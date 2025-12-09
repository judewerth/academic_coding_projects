# %% Import Modules
import os
if os.path.basename(os.getcwd()) == "notebooks":
    os.chdir("..")
import datajoint as dj
from datetime import datetime
from pathlib import Path
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import workflow
from matplotlib import pyplot as plt
import spikeinterface as si
from spikeinterface import widgets, exporters, postprocessing, qualitymetrics, sorters
from spikeinterface import preprocessing
import spikeinterface.extractors as se
import probeinterface as pi
from probeinterface import plotting
from workflow.pipeline import *
from workflow.utils.ingestion_utils import El2ROW
from workflow.utils.paths import (
    get_ephys_root_data_dir,
    get_raw_root_data_dir,
    get_processed_root_data_dir,
)
from element_interface.utils import dict_to_uuid, find_full_path, find_root_directory, _to_Path

#%% Functions
class spike:

    # Global Parameters
    chan = 32; # channels (per organoid)
    org = 4; # number of organoids
    # el_in = { # number electrodes inside of the organoid 
    #     "O09":32 , "O10":16 , "O11":20 , "O12":14 , 
    #     "O13":25 , "O14":13 , "O15":11 , "O16":11 ,
    #     "O17":22 , "O18":19 , "O19":20 , "O20":17
    # }
    el_out = { # number electrodes outside of the organoid 
        "O09":[] , "O10":16 , "O11":12 , "O12":18 , 
        "O13":7 , "O14":19 , "O15":21 , "O16":21 ,
        "O17":10 , "O18":13 , "O19":12 , "O20":15
    }


    def __init__(self , session_info , session_probe_info):
        # Update query // Datajoint database

        # Make sure organoid_id in both are equal
        assert session_info["organoid_id"] == session_probe_info["organoid_id"] , "organoid_id values need to agree"

        # Check to make sure duration is appropriate length
        SPIKE_SORTING_DURATION = 120 # minutes
        start_time = datetime.strptime(session_info["start_time"], '%Y-%m-%d %H:%M:%S')
        end_time = datetime.strptime(session_info["end_time"], '%Y-%m-%d %H:%M:%S')
        duration = (end_time - start_time).total_seconds() / 60
        assert session_info["session_type"] == "spike_sorting" and duration <= SPIKE_SORTING_DURATION, \
                f"Session type must be 'spike_sorting' and duration must be less than {SPIKE_SORTING_DURATION} minutes"

        # Set Ephys Session Info
        # ephys.EphysSession.insert1(session_info, ignore_extra_fields=True, skip_duplicates=True)
        query = culture.Experiment().proj("drug_name") * ephys.EphysSession & {"session_type": "spike_sorting"} 
        self.key = (query & session_info).fetch1()

        # Set Ephys Session Probe Info
        #   Insert Probe Info into database
        # ephys.EphysSessionProbe.insert1(session_probe_info, ignore_extra_fields=True, skip_duplicates=True)
        self.probe_info = (ephys.EphysSessionProbe & self.key).fetch1()
        self.probe_type = ((probe.Probe * ephys.EphysSessionProbe()) & self.key).fetch1("probe_type")

        # Assign used electrodes and organoids to probe_info
        self.probe_info["unused_electrodes"] = session_probe_info["unused_electrodes"]   
        self.unused_organoids = session_info["unused_organoids"]
        

    def get_files(self):
    
        # Get File Names
        files = (
            ephys.EphysRawFile
            & f"file_time BETWEEN '{self.key['start_time']}' AND '{self.key['end_time']}'"
        ).fetch("file_path", order_by="file_time")

        # Prepare file names
        # Remove Folder before in files
        f = 0
        for file in files:
            file_list = file.split("/")
            name = file_list[-1]

            name = name.removeprefix("processed_")
            
            files[f] = name
            f += 1
        return files

    def get_recording(self , files):

        # Get and Save Recording Data
        # Find data root folderss
        data_roots = []
        for folder in os.listdir(get_ephys_root_data_dir()):
            data_roots.append(get_raw_root_data_dir() / folder)

        # Prepare to Record Files
        stream_name = "RHD2000 amplifier channel"
        datapath = get_processed_root_data_dir() / "recording"
        datapath.mkdir(exist_ok=True , parents=True)
        recording_list = []

        for file in files:

            filename = file[0:-4] + ".pkl"
            filedata = None

            # Get Recording Data for Single File
            if (datapath / filename).exists():
                filedata = si.load_extractor(datapath / filename)
            else:
                filepath = find_full_path(data_roots , file)
                print(f"Processsing {file} Recording.")

                filedata = si.extractors.read_intan(filepath , stream_name=stream_name)
                filedata.dump_to_pickle(file_path = datapath / filename)

            recording_list.append(filedata)

        recording = si.concatenate_recordings(recording_list = recording_list)

        return recording

    def filter_recording(self , raw_recording):
        
        # Find organoid batch    
        if self.probe_info["organoid_id"] in ["O09" , "O10" , "O11" , "O12"]:
            batch_organoids = ["O09" , "O10" , "O11" , "O12"]
        elif self.probe_info["organoid_id"] in ["O13" , "O14" , "O15" , "O16"]:
            batch_organoids = ["O13" , "O14" , "O15" , "O16"]
        elif  self.probe_info["organoid_id"] in ["O17" , "O18" , "O19" , "O20"]:
            batch_organoids = ["O17" , "O18" , "O19" , "O20"]
        else:
            raise TypeError("Invalid organoid_id")
        
        # Assign groups to channels based on organoid
        group_assignments = [organoid for organoid in batch_organoids for _ in range(self.chan)] # 128 str list 
        raw_recording.set_channel_groups(groups=group_assignments)
        group_recordings = raw_recording.split_by('group') # list of raw_recording seperated by organoid
        
        # Remove unused organoids 
        if self.unused_organoids:
            for organoid in self.unused_organoids:
                batch_organoids.pop(organoid)

        
        # Get organoid recordings and remove unused electrodes
        # recording_list = []
        self.probegroup = pi.ProbeGroup()
        for organoid in batch_organoids:
            
            # get recording for a single organoid
            organoid_recording = group_recordings[organoid]

            # find number of unused electrodes
            if not self.probe_info["unused_electrodes"]:
                el = self.el_out[organoid] 
            else:
                el = self.probe_info["unused_electrodes"]

            # Make list of unused electrodes
            if el:

                if isinstance(el , int):
                    unused_electrodes = El2ROW[el:]
                    unused_electrodes = unused_electrodes.tolist()
                elif isinstance(el , list):
                    unused_electrodes = el
                else:
                    raise TypeError("used_electrodes must be type, int or list")
                
            else:
                unused_electrodes = []

            # Advance unused electrodes (after org1 the channel ids are 65-128 but the unused electrodes are from 0-31)
            channel_advance = organoid_recording.get_channel_ids()[0].astype(int) # first value of "channel ids"
            unused_electrodes += channel_advance

            # Remove unused electrodes 
            raw_recording = raw_recording.remove_channels(
            remove_channel_ids=np.array([str(elec) for elec in unused_electrodes])
            )

            # Add probe to organoid recording
            probe = pi.generate_linear_probe(num_elec=self.chan-len(unused_electrodes) , ypitch=100 , contact_shape_params={'radius': 15})

            # Append onto global variables
            # recording_list.append(organoid_recording)
            self.probegroup.add_probe(probe)
        

        # Set global probe index
        self.probegroup.set_global_device_channel_indices(range(self.probegroup.get_contact_count()))
            
        # Filter the Recordings
        recording_f = si.preprocessing.bandpass_filter(recording=raw_recording, freq_min=300, freq_max=6000)
        processed_recording = si.preprocessing.common_reference(recording=recording_f, operator="median")
        processed_recording = processed_recording.set_probegroup(self.probegroup)

        return processed_recording


    def get_sorting(self , recording , sorter_name , savefolder):

        # Get and Save waveform data
        # Check if recording and files cover the same number of segments (minutes)

        datapath = get_processed_root_data_dir() / sorter_name

        if savefolder:
            savepath = (datapath / savefolder)
            savepath.mkdir(exist_ok=True , parents=True)

            if (savepath / "sorting.pkl").exists():
                sorting = si.load_extractor(savepath / "sorting.pkl")
            else:
                sorting = si.sorters.run_sorter_by_property(sorter_name=sorter_name , recording=recording , grouping_property='group' , working_folder=savepath , remove_existing_folder=True,
            verbose=True)
                sorting.dump_to_pickle(file_path= savepath / "sorting.pkl")
        else:
            savepath = (datapath / "test")
            savepath.mkdir(exist_ok=True , parents=True)

            sorting = si.sorters.run_sorter_by_property(sorter_name=sorter_name , recording=recording , grouping_property='group' , working_folder=savepath , remove_existing_folder=True,
            verbose=True)
            
        return sorting

    def get_waveforms(self , sorting, recording , savefolder):
        
        # Get and save waveform data
        datapath = get_processed_root_data_dir() / "waveforms"

        if savefolder:
            savepath = datapath / savefolder

            if savepath.exists():
                we = si.load_waveforms(savepath , with_recording=True)

            else:
                savepath.mkdir(exist_ok=True , parents=True)
                we = si.extract_waveforms(
                    recording,
                    sorting,
                    folder=savepath,
                    ms_before=1.5,
                    ms_after=2.,
                    max_spikes_per_unit=500,
                    # overwrite=True, 
                    )
        else:
            savepath = datapath / "test"
            savepath.mkdir(exist_ok=True , parents=True)

            we = si.extract_waveforms(
            recording,
            sorting,
            folder=savepath,
            ms_before=1.5,
            ms_after=2.,
            max_spikes_per_unit=500,
            overwrite=True, 
            )           

        return we
