from pyMCDS_ECM import *
import os, re, sys
import pandas as pd
import xml.etree.ElementTree as ET
from spheroid_area_function import *



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
        'dead_phagocytosis_rate', 'live_phagocytosis_rates', 'attack_rates', 'damage_rate', 'fusion_rates', 'transformation_rates',  'previous_position_x', 
        'previous_position_y', 'previous_position_z', 'point_on_membrane_x', 'point_on_membrane_y', 'point_on_membrane_z', 'ecm_stiffness','ribose_concentration']

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
    if cell_definitions:
        cell_definition = cell_definitions.find('cell_definition')
    if cell_definition:
        phenotype = cell_definition.find('phenotype')
    if phenotype:
        cycle = phenotype.find('cycle')
        motility = phenotype.find('motility')
        mechanics = phenotype.find('mechanics')
    if cycle:
        phase_transition_rates = cycle.find('phase_transition_rates')

    if mechanics:
        cell_cell_adhesion_strength = mechanics.find('cell_cell_adhesion_strength')
        if cell_cell_adhesion_strength:
            cell_adh = cell_cell_adhesion_strength.text
            cell_adh = float(cell_adh) # type: ignore
            cell_df['cell_adh'] = cell_adh
        cell_cell_repulsion_strength = mechanics.find('cell_cell_repulsion_strength')
        if cell_cell_repulsion_strength:
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


    index = int(snapshot.replace("output", "") )
    # print(f'{index=}',flush=True)

    cell_df = cell_df.set_index(pd.Series([index]*len(cell_df.index)))

    # print('data frame \n', flush=True)
    # print(cell_df, flush = True)


    return cell_df


def simulation_data(data_folder_dir,simulation,ribose,seed):
    replace = False
 
    print(f"\n### Simulation {simulation}, ribose {ribose}, seed {seed} ###\n", flush=True)

    if((os.path.exists(data_folder_dir + f'/output_rib{ribose}_{simulation}_{seed}/dataframe_rib{ribose}_{simulation}_{seed}.pkl')) and (replace==False)):
        print(f"\nDataframe exists\n", flush=True)
        return
    
    else:
        if(replace):
            print(f"\n!!! Replace {replace} !!!\n", flush=True)

        # Data folder
        data_folder = data_folder_dir + f'/output_rib{ribose}_{simulation}_{seed}/'
        
        files = os.listdir(data_folder)
        files.sort()
        # print(files)

        df = pd.DataFrame()
        df_list = []
        
        for i in range(len(files)):
            if not re.search('_ECM\.mat', files[i]):
                continue
            if files[i].split('_')[0]!='output00000500':
                continue
            file = files[i].split('_')[0]
            # print(f'Processing snapshot {file}', flush=True )
            # Call function to get time and average total speeds values
            df_new = snapshot_data(ribose, simulation, seed, files[i].split('_')[0], data_folder)

            # df = pd.read_pickle(f'../dataframes_pickle/dataframe_rib{ribose}_{simulation}_{seed}.pkl')
            df_list.append(df_new)

            # df = pd.concat([df, df_new], copy=False, axis=0)

            # size = df.memory_usage().sum()
            # print(f'data frame size {size/1024**3} \n', flush=True)
            
        
        df = pd.concat(df_list, copy=False, axis=0)
        
        # pd.set_option('display.max_rows', None)
        # pd.set_option('display.max_columns', None)
        
        print('data frame \n', flush=True)
        print(df, flush = True)

        df.to_pickle(data_folder + f'dataframe_rib{ribose}_{simulation}_{seed}.pkl')

        return


