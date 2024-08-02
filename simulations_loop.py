import xml.etree.ElementTree as ET
import os, sys
from pathlib import Path
import numpy as np
import subprocess

repeat_simulations = False
proj = 'ecm_density' # 'test' # 

if repeat_simulations:
    seeds = [0,1,2,3,4,5,6,7,8,9]
    riboses = [50,200]
    simulations = [18]

    return_code = os.system(f'make reset; make clean; make load PROJ={proj}')
    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()
    # sys.exit()


else:
    simulation_id = 292
    #### Define parameter values
    seeds = [0,1,2,3,4,5,6,7,8,9]
    riboses = [0,50,200] # [0] # 

    proliferations_val = [0.00037] # [0.00032] # [0.0002,0.0003,0.0004] # 
    mot_speeds_val = [0.6] # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] # [0.6] # [0.1] # [0.5]# 
    
    repulsions_val = [10]#,20,40,80] #[10] # 

    adhesions_val = [0.4]#,0.8,1.6,3.2] # [0.4] # 
    
    if len(adhesions_val) != len(repulsions_val):
        print(f'Adhesion and repulsion lists have different lengths', flush=True)
        sys.exit() 

    r_density_val_temp = [0.0002] # [0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064,0.0128] # [0.005] # [0.0004] # [0.0004,0.01] # 
    # r_orientation_val_temp = [0]
    # r_anisotropy_val_temp = [0] 

    r_density_val = r_density_val_temp #* len(r_orientation_val_temp) * len(r_anisotropy_val_temp)
    # r_orientation_val = np.repeat(r_orientation_val_temp * len(r_anisotropy_val_temp), len(r_density_val_temp))
    # r_anisotropy_val = np.repeat(r_anisotropy_val_temp, len(r_density_val_temp) * len(r_orientation_val_temp))

    # alpha_val = [1,2,4,8,16,32,64]
    # beta_val = [1,2,4,8,16,32,64]

    # alphas = np.repeat(alpha_val, len(beta_val))
    # betas = beta_val * len(alpha_val)

    alphas = [32]
    betas = [8]
    
    proliferations = np.repeat(proliferations_val * len(alphas), len(mot_speeds_val) * len(adhesions_val) * len(r_density_val))
    mot_speeds = np.repeat(mot_speeds_val * len(proliferations_val)  * len(alphas) ,  len(adhesions_val)  * len(r_density_val)) 
    repulsions = repulsions_val * len(proliferations_val) * len(mot_speeds_val) * len(r_density_val) * len(alphas)
    adhesions = adhesions_val * len(proliferations_val) * len(mot_speeds_val) * len(r_density_val)  * len(alphas) 
    r_densities = np.repeat(r_density_val * len(proliferations_val) * len(mot_speeds_val) * len(alphas), len(adhesions_val) ) 
    # r_orientations = r_orientation_val * len(proliferations_val) * len(mot_speeds_val) * len(adhesions_val) 
    # r_anisotropies = r_anisotropy_val * len(proliferations_val) * len(mot_speeds_val) * len(adhesions_val) 

    alphas = np.repeat(alphas,len(proliferations))
    betas = np.repeat(betas,len(proliferations))
    

    length_list_parameters = [len(proliferations),len(repulsions),len(adhesions),len(mot_speeds),len(r_densities), len(alphas), len(betas)] #,len(r_orientations),len(r_anisotropies)]
    print(f'Parameter arrays lengths {length_list_parameters}', flush=True)

    length_list_parameters = [i for i in length_list_parameters if i != 0]
    length_list_parameters = list(np.unique(length_list_parameters))
    print(f'Mot speeds: {mot_speeds}\n')
    print(f'Prolif: {proliferations}\n')
    print(f'Rep: {repulsions}\n')
    print(f'Adh: {adhesions}\n')
    print(f'Degradation: {r_densities}\n')
    print(f'Alpha  : {alphas}\n')
    print(f'Beta: {betas}\n')


    if len(length_list_parameters)>1:
        print(f'Parameter arrays length no match', flush=True)
        sys.exit() 

    num_simulations = len(mot_speeds)
    # sys.exit() 

    ### Start the experiment with the given ribose and seed in the settings file ###
    return_code = os.system(f'make reset; make clean; make load PROJ={proj}')

    if return_code != 0:
        print(f'Failed with exit code: {return_code}')
        sys.exit()

    if simulation_id < 0:
        directories = []

        #### Check what is the last simulation number
        for p in Path(f'./data/{proj}').iterdir():
            if p.is_dir():
                directories.append(p.name)
        if not directories:
            simulation_id = 0
        else:
            simulation_ids = [int(name.split('_')[2]) for name in directories]
            simulation_ids.sort()
            simulation_id = max(simulation_ids) + 1

    print(f'simulations: {simulation_id} to {simulation_id + num_simulations -1} ')

    simulations = range(simulation_id, simulation_id + num_simulations)

results = []
i = 0

for sim in simulations:
    
    if repeat_simulations == False:
        proliferation = proliferations[i]
        motility_speed = mot_speeds[i] 
        adhesion = adhesions[i]
        repulsion = repulsions[i]
        # r_orientation = 0
        r_density = r_densities[i]
        # r_anisotropy = 0
        alpha_text = alphas[i]
        beta_text = betas[i]

    for rib in riboses: 
        for seed in seeds:
            print(f'\n####Simulation {sim}, ribose {rib}, random seed {seed}', flush=True)
            
            return_code = os.system(f'mkdir -p ./data/{proj}/output_rib{rib}_{sim}_{seed}')
            if return_code != 0:
                print(f'Failed with exit code: {return_code}')
                sys.exit()  

            if repeat_simulations:
                ### Start the experiment with the given ribose and seed in the settings file ###
                
                #### Copy PhysiCell_settings file from existing simulation with ribose 0 and seed 0 into config folder
                return_code = os.system(f'cp ./data/{proj}/output_rib0_{sim}_{seed}/PhysiCell_settings.xml ./config/')
                if return_code != 0:
                    print(f'Failed with exit code: {return_code}')
                    sys.exit()
                
                #### Copy custom.cpp file from existing simulation with ribose 0 and seed 0 into custom_modules folder
                return_code = os.system(f'cp ./user_projects/{proj}/custom_modules/* ./custom_modules/ ')
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
            tumor_radius = user_parameters.find('tumor_radius')
            folder = save.find('folder')   
            cell_definition = cell_definitions.find('cell_definition')
            phenotype = cell_definition.find('phenotype')
            custom_data = cell_definition.find('custom_data')
            mechanics = phenotype.find('mechanics')
            cycle = phenotype.find('cycle')
            volume = phenotype.find('volume')
            motility = phenotype.find('motility')
            total = volume.find('total')
            cell_cell_adhesion_strength = mechanics.find('cell_cell_adhesion_strength')
            cell_cell_repulsion_strength = mechanics.find('cell_cell_repulsion_strength')
            phase_transition_rates = cycle.find('phase_transition_rates')
            prolif_rate = phase_transition_rates.find('rate')
            mot_speed = motility.find('speed')
            fiber_realignment_rate = custom_data.find('fiber_realignment_rate')
            ecm_density_rate = custom_data.find('ecm_density_rate')
            anisotropy_increase_rate = custom_data.find('anisotropy_increase_rate')
            ecm_sensitivity = custom_data.find('ecm_sensitivity')
            alpha = user_parameters.find('alpha')
            beta = user_parameters.find('beta')
        
            if(repeat_simulations):
                virtual_wall_at_domain_edge.text = 'true'
                random_seed.text = str(seed) 
                ribose_concentration.text = str(rib)
                folder.text = f'data/{proj}/output_rib{rib}_{sim}_{seed}/'

            else:
                random_seed.text = str(seed) 
                folder.text = f'data/{proj}/output_rib{rib}_{sim}_{seed}/'
                ribose_concentration.text = str(rib)
                cell_cell_adhesion_strength.text = str(adhesion)
                cell_cell_repulsion_strength.text = str(repulsion)
                prolif_rate.text = str(proliferation)
                mot_speed.text = str(motility_speed)
                fiber_realignment_rate.text = '0'#str(r_orientation)
                ecm_density_rate.text = str(r_density)
                anisotropy_increase_rate.text = '0'#str(r_anisotropy)
                # ecm_sensitivity.text = '0'
                alpha.text = str(alpha_text)
                beta.text = str(beta_text)
                initial_anisotropy.text = '0'
                # rate.text = str(r_f0)
                # rate.text = str(r_density)
                tumor_radius.text = '100'
                # total.text = '1900' # '4027' # 

                ecm_orientation_setup.text =  'random' # # 'horizontal' # 'starburst' # 'random' # 'circular' # 
            
            #### WRITE THE XML FILE FOR THE CURRENT SIMULATIONS ####
            tree.write('./config/PhysiCell_settings.xml')    # this should be the path too where pyhsicell wants it's settings
            
            print(f'#### Xml settings file has been written\n',flush=True)

            #### Copy custom.cpp file from custom_modules folder into new simulation folder
            os.system(f'cp ./custom_modules/custom.cpp ./data/{proj}/output_rib{rib}_{sim}_{seed}/')

            #### Copy PhysiCell_settings file from config folder into new simulation folder
            os.system(f'cp ./config/PhysiCell_settings.xml ./data/{proj}/output_rib{rib}_{sim}_{seed}/')

            #### START SIMULATION ####
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

# Iterate through all return codes, only exit when all are done. sp.run always
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

