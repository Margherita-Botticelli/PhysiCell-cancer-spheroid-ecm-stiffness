import xml.etree.ElementTree as ET
import os, sys
from pathlib import Path
import numpy as np
import subprocess

#### Script to run multiple simulations automatically by giving the parameter values

#### Project name corresponding to folder in user_projects
proj = 'ecm_density' # 'test' # 'ecm_fibres'

#### True if we want to repeat an existing simulation with different random seeds or ribose concentrations 
repeat_simulations = False

if repeat_simulations:
    seeds = [0,1,2,3,4,5,6,7,8,9]
    riboses = [50,200]
    simulations = [0]

    return_code = os.system(f'make reset; make clean; make load PROJ={proj}')
    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()
    # sys.exit()

else:
    #### Define parameter values for the simulations
    seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] # [0] # 
    riboses = [0] # [0, 50, 200] # 
    proliferation_vals = [0.0002, 0.0003, 0.0004] # [0.00032] # [0.00037] # 
    mot_speed_vals = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6] # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] # [0.1] # [0.6] # [0.6] # 
    repulsion_vals = [10]
    adhesion_vals = [0.4]
    r_density_vals = [0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032, 0.0064, 0.0128] # [0.001] # [0.01] # 
    sigma_vals = [0] # [0, 0.001]
    delta_vals = [0] # [0, 0.001]

    #### Combine parameters
    sigmas = np.repeat(sigma_vals, len(delta_vals))
    deltas = delta_vals * len(sigma_vals)
    proliferations = np.repeat(proliferation_vals * len(sigmas), len(mot_speed_vals) * len(adhesion_vals) * len(r_density_vals))
    mot_speeds = np.repeat(mot_speed_vals * len(proliferation_vals)  * len(sigmas) ,  len(adhesion_vals)  * len(r_density_vals)) 
    repulsions = repulsion_vals * len(proliferation_vals) * len(mot_speed_vals) * len(r_density_vals) * len(sigmas)
    adhesions = adhesion_vals * len(proliferation_vals) * len(mot_speed_vals) * len(r_density_vals)  * len(sigmas) 
    r_densities = np.repeat(r_density_vals * len(proliferation_vals) * len(mot_speed_vals) * len(sigmas), len(adhesion_vals) ) 
    sigmas = np.repeat(sigmas,len(proliferations)* len(delta_vals))
    deltas = np.repeat(deltas * len(sigma_vals),len(proliferations))
    
    #### Check if the parameters all have same length
    length_list_parameters = [len(proliferations),len(repulsions),len(adhesions),len(mot_speeds),len(r_densities), len(sigmas), len(deltas)] #,len(r_orientations),len(r_anisotropies)]
    print(f'Parameter arrays lengths {length_list_parameters}', flush=True)

    length_list_parameters = [i for i in length_list_parameters if i != 0]
    length_list_parameters = list(np.unique(length_list_parameters))
    
    print(f'Mot speeds: {mot_speeds}\n')
    print(f'Prolif: {proliferations}\n')
    print(f'Rep: {repulsions}\n')
    print(f'Adh: {adhesions}\n')
    print(f'Degradation: {r_densities}\n')
    print(f'sigma  : {sigmas}\n')
    print(f'delta: {deltas}\n')

    #### If not exit
    if len(length_list_parameters)>1:
        print(f'Parameter arrays length no match', flush=True)
        sys.exit() 

    #### Get total number of simulations
    num_simulations = len(mot_speeds)

    #### Reset, clean and initiate the project (make sure the project folder exists in data folder too)
    return_code = os.system(f'make reset; make clean; make load PROJ={proj}')
    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()

    #### Simulations ID number: write number manually or write negative number to find last simulation number in the data folder 
    simulation_id = 0

    #### If the simulation_id number is negative find the last simulation number in the data folder
    if simulation_id < 0:
        directories = []

        #### Check what is the last simulation number in the data folder
        for p in Path(f'./data/{proj}').iterdir():
            if p.is_dir():
                directories.append(p.name)
        if not directories:
            simulation_id = 0
        else:
            simulation_ids = [int(name.split('_')[2]) for name in directories]
            simulation_ids.sort()
            simulation_id = max(simulation_ids) + 1

    #### List of all simulations
    simulations = range(simulation_id, simulation_id + num_simulations)
    print(f'simulations: {simulation_id} to {simulation_id + num_simulations -1} ')

#### List to check when all of the simulations are finished
results = []

#### Initiate parameter to go through parameter lists every new simulation loop
i = 0

for sim in simulations:
    if repeat_simulations == False:
        #### Select parameter values
        proliferation = proliferations[i]
        motility_speed = mot_speeds[i] 
        adhesion = adhesions[i]
        repulsion = repulsions[i]
        r_density = r_densities[i]
        sigma_text = sigmas[i]
        delta_text = deltas[i]

    for rib in riboses: 
        for seed in seeds:
            print(f'\n####Simulation {sim}, ribose {rib}, random seed {seed}', flush=True)
            
            #### Create the output folder for the current simulation
            return_code = os.system(f'mkdir -p ./data/{proj}/output_rib{rib}_{sim}_{seed}')
            if return_code != 0:
                print(f'Failed with exit code: {return_code}')
                sys.exit()  

            if repeat_simulations:
                ### Start the experiment with the given ribose and seed in the settings file ###
                
                #### Copy PhysiCell_settings file from existing simulation with ribose 0 into config folder
                return_code = os.system(f'cp ./data/{proj}/output_rib0_{sim}_{seed}/PhysiCell_settings.xml ./config/')
                if return_code != 0:
                    print(f'Failed with exit code: {return_code}')
                    sys.exit()
                
                #### Copy custom.cpp file from existing simulation with ribose 0 into custom_modules folder
                return_code = os.system(f'cp ./user_projects/{proj}/custom_modules/* ./custom_modules/ ')
                if return_code != 0:
                    print(f'Failed with exit code: {return_code}')
                    sys.exit()
            else:
                print(f'#### {motility_speed=}, {proliferation=}, {adhesion=}, {repulsion=}, {r_density=} ####\n',flush=True)

            #### Select PhysiCell_settings file from config folder 
            tree = ET.parse('./config/PhysiCell_settings.xml')  
            root = tree.getroot()

            #### Build the PhysiCell_settings file with the given parameter values
            options = root.find('options')
            user_parameters = root.find('user_parameters')
            save = root.find('save')
            cell_definitions = root.find('cell_definitions')
            virtual_wall_at_domain_edge = options.find('virtual_wall_at_domain_edge') # type: ignore
            random_seed = user_parameters.find('random_seed') # type: ignore
            ribose_concentration = user_parameters.find('ribose_concentration') # type: ignore
            ecm_orientation_setup = user_parameters.find('ecm_orientation_setup') # type: ignore
            initial_anisotropy = user_parameters.find('initial_anisotropy') # type: ignore
            tumor_radius = user_parameters.find('tumor_radius') # type: ignore
            folder = save.find('folder')    # type: ignore
            cell_definition = cell_definitions.find('cell_definition') # type: ignore
            phenotype = cell_definition.find('phenotype') # type: ignore
            custom_data = cell_definition.find('custom_data') # type: ignore
            mechanics = phenotype.find('mechanics') # type: ignore
            cycle = phenotype.find('cycle') # type: ignore
            volume = phenotype.find('volume') # type: ignore
            motility = phenotype.find('motility') # type: ignore
            total = volume.find('total') # type: ignore
            cell_cell_adhesion_strength = mechanics.find('cell_cell_adhesion_strength') # type: ignore
            cell_cell_repulsion_strength = mechanics.find('cell_cell_repulsion_strength') # type: ignore
            phase_transition_rates = cycle.find('phase_transition_rates') # type: ignore
            prolif_rate = phase_transition_rates.find('rate') # type: ignore
            mot_speed = motility.find('speed') # type: ignore
            fiber_realignment_rate = custom_data.find('fiber_realignment_rate') # type: ignore
            ecm_density_rate = custom_data.find('ecm_density_rate') # type: ignore
            anisotropy_increase_rate = custom_data.find('anisotropy_increase_rate') # type: ignore
            ecm_sensitivity = custom_data.find('ecm_sensitivity') # type: ignore
            sigma = user_parameters.find('sigma') # type: ignore
            delta = user_parameters.find('delta') # type: ignore
        
            if(repeat_simulations):
                virtual_wall_at_domain_edge.text = 'true' # type: ignore
                random_seed.text = str(seed)  # type: ignore
                ribose_concentration.text = str(rib) # type: ignore
                folder.text = f'data/{proj}/output_rib{rib}_{sim}_{seed}/' # type: ignore

            else:
                random_seed.text = str(seed)  # type: ignore
                folder.text = f'data/{proj}/output_rib{rib}_{sim}_{seed}/' # type: ignore
                ribose_concentration.text = str(rib) # type: ignore
                cell_cell_adhesion_strength.text = str(adhesion) # type: ignore
                cell_cell_repulsion_strength.text = str(repulsion) # type: ignore
                prolif_rate.text = str(proliferation) # type: ignore
                mot_speed.text = str(motility_speed) # type: ignore
                fiber_realignment_rate.text = '0' # type: ignore
                ecm_density_rate.text = str(r_density) # type: ignore
                anisotropy_increase_rate.text = '0' # type: ignore
                # ecm_sensitivity.text = '0'
                sigma.text = str(sigma_text) # type: ignore
                delta.text = str(delta_text) # type: ignore
                initial_anisotropy.text = '0' # type: ignore
                tumor_radius.text = '100' # type: ignore
                ecm_orientation_setup.text =  'random' # 'horizontal' # 'starburst' # 'random' # 'circular' #  # type: ignore
            
            #### Write the xml file for the current simulation
            tree.write('./config/PhysiCell_settings.xml')    
            print(f'#### Xml settings file has been written\n',flush=True)

            #### Copy custom.cpp file from custom_modules folder into new simulation folder
            os.system(f'cp ./custom_modules/custom.cpp ./data/{proj}/output_rib{rib}_{sim}_{seed}/')

            #### Copy PhysiCell_settings file from config folder into new simulation folder
            os.system(f'cp ./config/PhysiCell_settings.xml ./data/{proj}/output_rib{rib}_{sim}_{seed}/')

            #### Start simulation (simulations run in parallel)
            os.system('make') 
            if len(results)==0:
                return_code = subprocess.Popen(f'./project')
                results.append(return_code)
            elif len(results)%19==0:
                return_code = subprocess.run(f'./project')
                results.append(return_code)
            else:
                return_code = subprocess.Popen(f'./project')
                results.append(return_code)

            os.system('sleep 5') 

    i += 1

# Iterate through all return codes, only exit when all are done (subprocess.run)
# returns True because we know it waits until finishing
while True:
    temp_results = []
    for x in results:
        if isinstance(x,subprocess.CompletedProcess):
            res = x.returncode
        else:
            res = x.poll()

        # If encounter a return code of not 0, exit program immedietly
        if not(res == None or res == 0):
            sys.exit(1)
            print('Non-zero exit code found, exiting...', flush=True)

        temp_results.append(res == 0)
    if all(temp_results):
        break

