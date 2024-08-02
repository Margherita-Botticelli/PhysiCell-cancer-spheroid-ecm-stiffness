from pyMCDS_ECM import *
import os, re, sys
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np
from spheroid_area_function import *
from delaunay_function import *

# Reduce memory usage
def reduce_mem_usage(df, verbose=True):
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    start_mem = df.memory_usage(deep=True).sum() / 1024**2    
    for col in df.columns:
        col_type = df[col].dtypes
        # print(f'{col_type}',flush=True)
        if col_type in numerics:
            c_min = df[col].min()
            c_max = df[col].max()
            if (str(col_type)[:3] == 'int') or (col in ['t', 'ID', 'maximum_number_of_attachments', 'persistence_time', 'overcrowding_threshold', 'attached_cells']):
                # print(f'{col=}',flush=True)
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)  
            else:  # for floats.
                if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                    df[col] = df[col].astype(np.float16)
                elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)    
    end_mem = df.memory_usage(deep=True).sum() / 1024**2
    if verbose: print('Mem. usage decreased to {:5.2f} Mb ({:.1f}% reduction)'.format(end_mem, 100 * (start_mem - end_mem) / start_mem), flush=True)
    return df

def snapshot_data(ribose, simulation, seed,snapshot,data_folder):
    #### Load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    #### Get current time
    t = ( mcds.get_time() )

    #### Cell data
    cell_df = mcds.get_cell_df()
    cell_var = mcds.get_cell_variables()

    # print(f"Remove cell variables from cell dataframe", flush=True)

    cell_var_remove = ['total_volume', 'cell_type', 'cycle_model', 'current_phase', 'elapsed_time_in_phase', 'nuclear_volume', 
        'cytoplasmic_volume', 'fluid_fraction', 'calcified_fraction', 'orientation_x', 'orientation_y', 'orientation_z', 'polarity', 'velocity_x', 'velocity_y', 
        'velocity_z', 'number_of_nuclei', 'damage', 'total_attack_time', 'contact_with_basement_membrane', 'dead', 'current_death_model', 'death_rates_x', 
        'death_rates_y', 'cytoplasmic_biomass_change_rate', 'nuclear_biomass_change_rate', 'fluid_change_rate', 'calcification_rate', 'target_solid_cytoplasmic', 
        'target_solid_nuclear', 'target_fluid_fraction', 'nuclear_radius', 'surface_area', 'cell_BM_adhesion_strength', 'cell_BM_repulsion_strength', 'cell_adhesion_affinities', 'attachment_elastic_constant', 'attachment_rate', 'detachment_rate', 'is_motile', 
        'migration_bias', 'chemotaxis_index', 'chemotaxis_direction', 'chemotactic_sensitivities', 'secretion_rates', 'uptake_rates', 
        'saturation_densities', 'net_export_rates', 'internalized_total_substrates', 'fraction_released_at_death', 'fraction_transferred_when_ingested', 
        'dead_phagocytosis_rate', 'live_phagocytosis_rates', 'attack_rates', 'damage_rate', 'fusion_rates', 'transformation_rates', 'ecm_stiffness','ribose_concentration']

    cell_df = cell_df.drop(cell_var_remove, axis=1)

    # print(f"Move cell dataframe columns", flush=True)

    dictionary = ['simulation','ribose','seed','t' ] + list(cell_df.columns)

    cell_df['t'] = t
    cell_df['simulation'] = simulation
    cell_df['seed'] = seed
    cell_df['ribose'] = ribose

    cell_df = cell_df[dictionary]   

    # print(f"Add xml settings file parameters to dataframe", flush=True)
    
    #### Add simulation parameters to dataframe
    tree = ET.parse(data_folder + 'PhysiCell_settings.xml')   
    root = tree.getroot()

    #### Simulation parameters indipendent of time
    user_parameters = root.find('user_parameters')
    cell_definitions = root.find('cell_definitions')
    cell_definition = cell_definitions.find('cell_definition')

    phenotype = cell_definition.find('phenotype')
    cycle = phenotype.find('cycle')
    motility = phenotype.find('motility')
    mechanics = phenotype.find('mechanics')
    phase_transition_rates = cycle.find('phase_transition_rates')

    alpha = float(user_parameters.find('alpha').text) # type: ignore
    beta = float(user_parameters.find('beta').text) # type: ignore
    cell_df['alpha'] = alpha
    cell_df['beta'] = beta


    cell_cell_adhesion_strength = mechanics.find('cell_cell_adhesion_strength')
    cell_adh = cell_cell_adhesion_strength.text
    cell_adh = float(cell_adh) # type: ignore
    cell_df['cell_adh'] = cell_adh
    cell_cell_repulsion_strength = mechanics.find('cell_cell_repulsion_strength')
    cell_rep = cell_cell_repulsion_strength.text
    cell_rep = float(cell_rep) # type: ignore
    cell_df['cell_rep'] = cell_rep

    prolif = float(phase_transition_rates.find('rate').text) # type: ignore
    max_mot_speed = float(motility.find('speed').text) # type: ignore
    
    cell_df['prolif'] = prolif
    cell_df['max_mot_speed'] = max_mot_speed

    # print(f"Compute spheroid area", flush=True)

    #### Call spheroid area function 
    spheroid_area = spheroid_area_function(cell_df)
    cell_df['spheroid_area'] = spheroid_area

    #### Call delaunday distance function
    delaunay_distance = delaunay_distance_function(cell_df)
    cell_df['delaunay_distance'] = delaunay_distance

    index = int(snapshot.replace("output", "") )
    # print(f'{index=}',flush=True)

    cell_df = cell_df.set_index(pd.Series([index]*len(cell_df.index)))

    # print('data frame \n', flush=True)
    # print(cell_df, flush = True)

    return cell_df


def simulation_data(data_folder_dir,simulation,ribose,seed,replace=False):

    print(f"\n#### {ribose=}, {simulation=}, {seed=} ####\n", flush=True)

    if((os.path.exists(data_folder_dir + f'output_rib{ribose}_{simulation}_{seed}/dataframe_rib{ribose}_{simulation}_{seed}.pkl')) and (replace==False)):
        print(f"\nDataframe exists\n", flush=True)
        return
    
    else:
        if(replace):
            print(f"\n!!! Replace {replace} !!!\n", flush=True)

        # Data folder
        data_folder = data_folder_dir + f'output_rib{ribose}_{simulation}_{seed}/'
        # print(f"{data_folder=}\n", flush=True)
        
        files = os.listdir(data_folder)
        files.sort()
        # print(files)

        df = pd.DataFrame()
        df_list = []
        
        for i in range(len(files)):
            if not re.search('_ECM\.mat', files[i]):
                continue
            file = files[i].split('_')[0]
            # print(f'Processing snapshot {file}', flush=True )
            # Call function to get time and average total speeds values
            df_new = snapshot_data(ribose, simulation, seed, files[i].split('_')[0], data_folder)

            df_new = reduce_mem_usage(df_new,verbose=False)

            # df = pd.read_pickle(f'../dataframes_pickle/dataframe_rib{ribose}_{simulation}_{seed}.pkl')
            df_list.append(df_new)


        df = pd.concat(df_list, copy=False, axis=0)
        
        # pd.set_option('display.max_rows', None)
        # pd.set_option('display.max_columns', None)
        
        # print('data frame \n', flush=True)
        # print(df, flush = True)

        df.to_pickle(data_folder + f'dataframe_rib{ribose}_{simulation}_{seed}.pkl')

        return 



