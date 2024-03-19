from cmath import phase
import xml.etree.ElementTree as ET
import os, sys
from pathlib import Path
from networkx import is_empty
import numpy as np
import subprocess

repeat_simulations = False

if repeat_simulations:
    seeds = [0]#[6,7,8,9,10,11]
    ribose = [0]
    simulations = [758]
    return_code = os.system(f'make reset; make clean; make load PROJ=simulations')
    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()


else:
    #### Define parameter values
    seeds = [0,1,2,3,4,5,6,7,8,9]
    ribose = [0]#, 50, 200]

    proliferations_val = [0.0007,0.00105,0.0014] #np.repeat([0.00035, 0.0007, 0.00105, 0.0014],6) # np.arange(0.02, 1, ).round(3)
    mot_speeds_val = [0.5,0.75,1.0] #[0.3, 0.5, 0.7, 0.9, 1.1] #* 4

    repulsions_val = [2] * 3 + [4] * 4 + [8] * 5 + [16] * 6 + [32] * 7 + [64] * 8 + [128] * 8
                      
    adhesions_val = [0.25, 0.5, 1,
                     0.25, 0.5, 1, 2, 
                     0.25, 0.5, 1, 2, 4, 
                     0.25, 0.5, 1, 2, 4, 8, 
                     0.25, 0.5, 1, 2, 4, 8, 16,
                     0.25, 0.5, 1, 2, 4, 8, 16, 32,
                     0.25, 0.5, 1, 2, 4, 8, 16, 32] 
    
    if len(adhesions_val) != len(repulsions_val):
        print(f'Adhesion and repulsion lists have different lengths', flush=True)
        sys.exit() 

    r_density_val_temp = [0.01] 
    # r_orientation_val_temp = [0]
    # r_anisotropy_val_temp = [0] 

    r_density_val = r_density_val_temp #* len(r_orientation_val_temp) * len(r_anisotropy_val_temp)
    # r_orientation_val = np.repeat(r_orientation_val_temp * len(r_anisotropy_val_temp), len(r_density_val_temp))
    # r_anisotropy_val = np.repeat(r_anisotropy_val_temp, len(r_density_val_temp) * len(r_orientation_val_temp))

    mot_speeds = np.repeat(mot_speeds_val * len(proliferations_val), len(adhesions_val)) #* len(r_orientation_val) #np.repeat(proliferations_val,len(mot_speeds_val))
    proliferations = np.repeat(proliferations_val,  len(adhesions_val) * len(mot_speeds_val)) #* len(r_orientation_val))
    repulsions = repulsions_val * len(proliferations_val) * len(mot_speeds_val) #* len(r_orientation_val)
    adhesions = adhesions_val * len(proliferations_val) * len(mot_speeds_val) #* len(r_orientation_val)
    r_densities = r_density_val * len(proliferations_val) * len(mot_speeds_val) * len(adhesions_val) 
    # r_orientations = r_orientation_val * len(proliferations_val) * len(mot_speeds_val) * len(adhesions_val) 
    # r_anisotropies = r_anisotropy_val * len(proliferations_val) * len(mot_speeds_val) * len(adhesions_val) 

    length_list_parameters = [len(proliferations),len(repulsions),len(adhesions),len(mot_speeds),len(r_densities)] #,len(r_orientations),len(r_anisotropies)]
    print(f'Parameter arrays lengths {length_list_parameters}', flush=True)

    length_list_parameters = [i for i in length_list_parameters if i != 0]
    length_list_parameters = list(np.unique(length_list_parameters))
    print(f'Mot speeds: {mot_speeds}\n')
    print(f'Prolif: {proliferations}\n')
    print(f'Rep: {repulsions}\n')
    print(f'Adh: {adhesions}\n')

    if len(length_list_parameters)>1:
        print(f'Parameter arrays length no match', flush=True)
        sys.exit() 

    num_simulations = len(mot_speeds)
    # sys.exit() 

    ### Start the experiment with the given ribose and seed in the settings file ###
    return_code = os.system(f'make reset; make clean; make load PROJ=simulations')

    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()

    directories = []

    #### Check what is the last simulation number
    for p in Path('./user_projects/ribose_0/output_density/').iterdir():
        if p.is_dir():
            directories.append(p.name)
    if not directories:
        simulation_id = 0
    else:
        simulation_ids = [int(name.split('_')[1]) for name in directories]
        simulation_ids.sort()
        simulation_id = 1#max(simulation_ids) + 1

    print(f'simulations: {simulation_id} to {simulation_id + num_simulations -1} ')

    simulations = range(simulation_id, simulation_id + num_simulations)

results = []
i = 0

for simulation in simulations:
    
    if repeat_simulations == False:
        proliferation = proliferations[i]
        motility_speed = mot_speeds[i] 
        adhesion = adhesions[i]
        repulsion = repulsions[i]
        # r_orientation = 0
        r_density = r_densities[i]
        # r_anisotropy = 0

    for rib in ribose: 
        for seed in seeds:
            print(f'\n####Simulation {simulation}, ribose {rib}, random seed {seed}', flush=True)
            
            return_code = os.system(f'mkdir -p ./user_projects/ribose_{rib}/output_density/output_{simulation}_{seed}')
            if return_code != 0:
                print(f'Failed with exit code: {return_code}')
                sys.exit()  

            if repeat_simulations:
                ### Start the experiment with the given ribose and seed in the settings file ###
                
                #### Copy PhysiCell_settings file from existing simulation with ribose 0 and seed 0 into config folder
                return_code = os.system(f'cp ./user_projects/ribose_0/output_density/output_{simulation}_0/PhysiCell_settings.xml ./config/')
                if return_code != 0:
                    print(f'Failed with exit code: {return_code}')
                    sys.exit()
                
                #### Copy custom.cpp file from existing simulation with ribose 0 and seed 0 into custom_modules folder
                return_code = os.system(f'cp ./user_projects/simulations/custom_modules/* ./custom_modules/ ')
                if return_code != 0:
                    print(f'Failed with exit code: {return_code}')
                    sys.exit()
            else:
                print(f'#### {motility_speed=}, {proliferation=}, {adhesion=}, {repulsion=}, {r_density=} ####\n',flush=True)

            #### Select PhysiCell_settings file from config folder 
            tree = ET.parse('./config/PhysiCell_settings.xml')  
            root = tree.getroot()

            #### BUILD THE SETTINGS FILE ####
            options = root.find('options')
            user_parameters = root.find('user_parameters')
            save = root.find('save')
            cell_definitions = root.find('cell_definitions')
            virtual_wall_at_domain_edge = options.find('virtual_wall_at_domain_edge')
            random_seed = user_parameters.find('random_seed')
            ribose_concentration = user_parameters.find('ribose_concentration')
            ecm_orientation_setup = user_parameters.find('ecm_orientation_setup')
            initial_anisotropy = user_parameters.find('initial_anisotropy')
            folder = save.find('folder')   
            cell_definition = cell_definitions.find('cell_definition')
            phenotype = cell_definition.find('phenotype')
            custom_data = cell_definition.find('custom_data')
            mechanics = phenotype.find('mechanics')
            cycle = phenotype.find('cycle')
            motility = phenotype.find('motility')
            cell_cell_adhesion_strength = mechanics.find('cell_cell_adhesion_strength')
            cell_cell_repulsion_strength = mechanics.find('cell_cell_repulsion_strength')
            phase_transition_rates = cycle.find('phase_transition_rates')
            prolif_rate = phase_transition_rates.find('rate')
            mot_speed = motility.find('speed')
            fiber_realignment_rate = custom_data.find('fiber_realignment_rate')
            ecm_density_rate = custom_data.find('ecm_density_rate')
            anisotropy_increase_rate = custom_data.find('anisotropy_increase_rate')
            ecm_sensitivity = custom_data.find('ecm_sensitivity')
        
            if(repeat_simulations):
                virtual_wall_at_domain_edge.text = 'true'
                random_seed.text = str(seed) 
                ribose_concentration.text = str(rib)
                folder.text = f'./user_projects/ribose_{rib}/output_density/output_{simulation}_{seed}'

            else:
                random_seed.text = str(seed) 
                folder.text = f'./user_projects/ribose_{rib}/output_density/output_{simulation}_{seed}'
                ribose_concentration.text = str(rib)
                cell_cell_adhesion_strength.text = str(adhesion)
                cell_cell_repulsion_strength.text = str(repulsion)
                prolif_rate.text = str(proliferation)
                mot_speed.text = str(motility_speed)
                # fiber_realignment_rate.text = str(r_orientation)
                ecm_density_rate.text = str(r_density)
                # anisotropy_increase_rate.text = str(r_anisotropy)
                # ecm_sensitivity.text = '0'
                # initial_anisotropy.text = '0'
                # rate.text = str(r_f0)
                # rate.text = str(r_density)

                # ecm_orientation_setup.text =   # 'horizontal' # 'starburst' # 'random' # 'circular' # 
            
            #### WRITE THE XML FILE FOR THE CURRENT SIMULATIONS ####
            tree.write('./config/PhysiCell_settings.xml')    # this should be the path too where pyhsicell wants it's settings
            
            print(f'#### Xml settings file has been written\n',flush=True)

            #### Copy custom.cpp file from custom_modules folder into new simulation folder
            os.system(f'cp ./custom_modules/custom.cpp ./user_projects/ribose_{rib}/output_density/output_{simulation}_{seed}')

            #### Copy PhysiCell_settings file from config folder into new simulation folder
            os.system(f'cp ./config/PhysiCell_settings.xml ./user_projects/ribose_{rib}/output_density/output_{simulation}_{seed}')

            #### START SIMULATION ####
            os.system('make') 
            return_code = subprocess.Popen(f'./project')
            results.append(return_code)

            os.system('sleep 5') 

    i += 1

while True:
    temp_results = []
    for x in results:
        res = x.poll()
        if not(res == None or res == 0):
            sys.exit(1)
        temp_results.append(res == 0)
    if all(temp_results):
        break

