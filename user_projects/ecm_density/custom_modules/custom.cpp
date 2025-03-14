/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include "./extracellular_matrix.h"
#include <cmath>

using namespace BioFVM; 
using namespace PhysiCell;

// extern ECM ecm;

ECM ecm;

void create_cell_types( void )
{
	// set the random seed 
	//SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.update_migration_bias = cell_ecm_adhesion_direction; 

	cell_defaults.functions.update_phenotype = NULL;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	// cell_defaults.functions.update_migration_bias = NULL; 
	// cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	// cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	// cell_defaults.functions.custom_cell_rule = custom_function; 
	// cell_defaults.functions.contact_function = contact_function; 
	// cell_defaults.functions.update_migration_bias = cell_ecm_adhesion_direction; 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	// cell_defaults.functions.update_velocity = custom_update_cell_velocity;
	
		
	// display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 

	initialize_microenvironment(); 

	return; 
}
void setup_extracellular_matrix(void)
{
	// DEPENDS ON MICROENVIRONMENT - CALL SETUP MICROENVIRONEMNT FIRST!!!!!

    // Resize the ECM mesh according to the microenvironment options
    ecm.ecm_mesh.resize(default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1], 
                        default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1], 
                        default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1], 
                        parameters.doubles("ecm_dx"), parameters.doubles("ecm_dy"), parameters.doubles("ecm_dz"));
    // Adjust the ECM units based on the resized ECM mesh
	ecm.resize_ecm_units_from_ecm_mesh();

	// ecm.ecm_mesh.display_information(std::cout );

    // Set up ECM alignment

    // <ecm_orientation_setup description="Specifies the initial ECM orientation: random, circular, starburt, oriented to the right, or oriented to the top" type="string" units="NA">circular</ecm_orientation_setup> parameters.string("ecm_orientation_setup")
	// std::vector<double> fiber_direction = { 1.0 , 0.0, 0.0 }; 
    // ecm_fiber_alignment.resize(microenvironment.mesh.voxels.size(), fiber_direction);  

    // Retrieve parameters for spacing and tumor radius
    double spacing = parameters.doubles("ecm_dx");
    double tumor_radius = parameters.doubles("tumor_radius");

    for(int n = 0; n < ecm.ecm_mesh.voxels.size(); n++)
	{
        // For random 2-D initialization
        if(parameters.strings("ecm_orientation_setup") == "random")
		{
            double theta = 6.2831853071795864769252867665590 * uniform_random(); // Random angle in radians
            ecm.ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0}; // Set ECM fiber alignment

            double x = spacing / 2; // Initialize x-coordinate
            double y = spacing / 2; // Initialize y-coordinate
            double x_outer; // Outer boundary of x-coordinate within the tumor

            double ecm_density = 0; // Initial ECM density

            int inner_n = 0; // Counter for loop iterations
            int voxel_index; // Index of the voxel in the microenvironment
            std::vector<double> pos; // Position vector for voxel

            // Iterate over y-coordinates within the tumor radius
            while(y < tumor_radius)
			{
                x = spacing / 2; // Reset x-coordinate for each y-level
            
                // Adjust x-coordinate step for odd/even iterations
                if(inner_n % 2 == 1)
                { 
                    x = 0.5 * spacing; 
                }
                x_outer = sqrt(tumor_radius * tumor_radius - y * y); // Calculate outer boundary for x-coordinate
                
                // Iterate over x-coordinates within the current y-level
                while(x < x_outer)
				{
                    // Set ECM density for the current position
                    pos = {x, y, 0.0};
                    voxel_index = microenvironment.nearest_voxel_index(pos);
                    // Check if the voxel index is valid before accessing
                    if (voxel_index >= 0 && voxel_index < ecm.ecm_voxels.size())
                    {
					ecm.ecm_voxels[voxel_index].density = ecm_density;
                    }
                    
                    // Set ECM density for the mirrored y-position if needed
                    if(fabs(y) > 0.01)
					{
                        pos = {x, -y, 0.0};
                        voxel_index = microenvironment.nearest_voxel_index(pos);
                        if (voxel_index >= 0 && voxel_index < ecm.ecm_voxels.size())
                        {
						ecm.ecm_voxels[voxel_index].density = ecm_density;
					}
                    }
                    
                    // Set ECM density for the mirrored x-position if needed
                    if(fabs(x) > 0.01)
					{ 
                        pos = {-x, y, 0.0};
                        voxel_index = microenvironment.nearest_voxel_index(pos);
                        if (voxel_index >= 0 && voxel_index < ecm.ecm_voxels.size())
                        {
						ecm.ecm_voxels[voxel_index].density = ecm_density;	
                        }
                        
                        // Set ECM density for the mirrored x and y positions if needed
                        if(fabs(y) > 0.01)
						{
                            pos = {-x, -y, 0.0};
                            voxel_index = microenvironment.nearest_voxel_index(pos);
                            if (voxel_index >= 0 && voxel_index < ecm.ecm_voxels.size())
                            {
							ecm.ecm_voxels[voxel_index].density = ecm_density;
						}
					}
				}
                    x += spacing; // Move to the next x-coordinate
			}
                y += spacing; // Move to the next y-coordinate
                inner_n++; // Increment the counter
		}
        }
		else
		{
            // If the orientation setup is not specified, halt the program
            std::cout << "WARNING: NO ECM ORIENTATION SPECIFIED. FIX THIS!!!" << std::endl;
            std::cout << "Halting program!!!" << std::endl;
			abort();
			return;
		}

        // Ensure index is within bounds before accessing
        if(n >= 0 && n < ecm.ecm_voxels.size())
        {
		ecm.ecm_voxels[n].density = parameters.doubles("initial_ecm_density");
		ecm.ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");
        }

        // Check if the initial density and anisotropy values are within valid bounds
        if(parameters.doubles("initial_ecm_density") > 1 || parameters.doubles("initial_anisotropy") > 1 ||
           parameters.doubles("initial_ecm_density") < 0 || parameters.doubles("initial_anisotropy") < 0)
		{
            std::cout << "WARNING: INITIAL DENSITY OR ANISOTROPY OUT OF BOUNDS! FIX THIS!" << std::endl;
            std::cout << "Halting program!!!" << std::endl;
			abort();
			return;
		}		
	}
}

void setup_tissue(void)
{
	// double Xmin = microenvironment.mesh.bounding_box[0]; 
	// double Ymin = microenvironment.mesh.bounding_box[1]; 
	// double Zmin = microenvironment.mesh.bounding_box[2]; 

	// double Xmax = microenvironment.mesh.bounding_box[3]; 
	// double Ymax = microenvironment.mesh.bounding_box[4]; 
	// double Zmax = microenvironment.mesh.bounding_box[5]; 
	
    // If simulating in 2D, adjust the Z bounds to be zero
    // if(default_microenvironment_options.simulate_2D == true)
	// {
    //     Zmin = 0.0; 
    //     Zmax = 0.0; 
	// }
	
	// double Xrange = Xmax - Xmin; 
	// double Yrange = Ymax - Ymin; 
	// double Zrange = Zmax - Zmin; 
	
    // Create cells of each type defined
	// Cell* pC;
	
    // for(int k = 0; k < cell_definitions_by_index.size(); k++)
	// {
    //     Cell_Definition* pCD = cell_definitions_by_index[k]; 
    //     std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
    //     for(int n = 0; n < parameters.ints("number_of_cells"); n++)
    //     {
    //         std::vector<double> position = {0, 0, 0}; 
    //         position[0] = Xmin + UniformRandom() * Xrange; 
    //         position[1] = Ymin + UniformRandom() * Yrange; 
    //         position[2] = Zmin + UniformRandom() * Zrange; 
    //         
    //         pC = create_cell(*pCD); 
    //         pC->assign_position(position);
    //     }
	// }
	// std::cout << std::endl; 
	
    // Load cells from your CSV file (if enabled)
	// load_cells_from_pugixml(); 	
	
    // Setting seed so cells always start with the same initial configuration
    SeedRandom(parameters.ints("random_seed"));  

    /******************************Single cell initialization*****************************/
	if(parameters.strings("cell_setup") == "single")
	{
		Cell* pCell;
		pCell = create_cell(cell_defaults);
		pCell->assign_position(0.0, 0.0, 0.0);
        // std::cin.get(); 
	}

	/******************************************2D Spheroid initialization***************************************/
	else if(parameters.strings("cell_setup") == "lesion")
	{
        // Place a cluster of tumor cells at the center
	
        // Get tumor radius from XML parameters
		double tumor_radius = parameters.doubles("tumor_radius"); 

        // Calculate cell spacing based on adhesion and repulsion parameters
		double cell_radius = cell_defaults.phenotype.geometry.radius;
		double relative_maximum_adhesion_distance = cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance;
		double cell_adhesion = cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength;
		double cell_repulsion = cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength;
        double sqrt_adhesion_to_repulsion_ratio = sqrt(cell_adhesion / cell_repulsion);
		double cell_spacing = (1 - sqrt_adhesion_to_repulsion_ratio);
        cell_spacing /= (0.5 * 1 / cell_radius - 0.5 * sqrt_adhesion_to_repulsion_ratio / (relative_maximum_adhesion_distance * cell_radius));

		Cell* pCell = NULL; 
		
        double x = 0.0; // Initialize x-coordinate
        double x_outer = tumor_radius; // Outer boundary for x-coordinate within the tumor
        double y = 0.0; // Initialize y-coordinate

        int n = 0; // Counter for loop iterations
        int voxel_index; // Index of the voxel in the microenvironment
        std::vector<double> pos; // Position vector for voxel

        // Iterate over y-coordinates within the tumor radius
        while(y < tumor_radius)
		{
			x = 0.0; 
            if(n % 2 == 1)
            { 
                x = 0.5 * cell_spacing; 
            }
            x_outer = sqrt(tumor_radius * tumor_radius - y * y); // Calculate outer boundary for x-coordinate
            
            // Iterate over x-coordinates within the current y-level
            while(x < x_outer)
			{
				pCell = create_cell(cell_defaults);
                pCell->assign_position(x, y, 0.0);
                set_single_behavior(pCell, "cycle entry", 0); // Set initial cell behavior
                
                // Set cells at mirrored y-position if needed
                if(fabs(y) > 0.01)
				{
					pCell = create_cell(cell_defaults);
                    pCell->assign_position(x, -y, 0.0);
                    set_single_behavior(pCell, "cycle entry", 0);
				}
				
                // Set cells at mirrored x-position if needed
                if(fabs(x) > 0.01)
				{ 
					pCell = create_cell(cell_defaults); 
                    pCell->assign_position(-x, y, 0.0);
                    set_single_behavior(pCell, "cycle entry", 0);
					
                    // Set cells at mirrored x and y positions if needed
                    if(fabs(y) > 0.01)
					{
						pCell = create_cell(cell_defaults);
                        pCell->assign_position(-x, -y, 0.0);
                        set_single_behavior(pCell, "cycle entry", 0);
					}
				}

                x += cell_spacing; // Move to the next x-coordinate
			}
			
            y += cell_spacing * sqrt(3.0) / 2.0; // Adjust y-coordinate for hexagonal packing
            n++; // Increment the counter
		}

        std::cout << "Cells placed in 2D lesion at the center of domain" << std::endl;
	}

	else
	{
        // If no cell setup is specified, halt the program
        std::cout << "WARNING!!! NO CELL SETUP SPECIFIED. SEE DOCUMENTATION and FIX" << std::endl;
        std::cout << "Halting program!!!" << std::endl;
		abort();
		return;
	}

	return;
}


std::vector<std::string> AMIGOS_invasion_coloring_function( Cell* pCell )
{
	// leaders are blue 
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		output[0] = "blue"; 
		output[2] = "blue"; 
		return output; 
	} 

	// followers are yellow except for the marker cells when running testing mode (then 20 % of cells are red)
    
	if( pCell->type == 2 )
    {
        output[2] = "yellow";	
		output[0] = "yellow";
		
		// Return yellow for followers and exit statement

		if(parameters.ints("unit_test_setup")==0)
		{return output;}

		// Return red for 20% of followers if unit test is called for	

		if( pCell->ID % 5 == 0)
		{
			output[0] = "red";
        	output[2] = "red";
        	return output;
		}

		// If doing unit testing AND cell not selected as marker cell, return blue. 
		output[2] = "blue";	
		output[0] = "blue";


        return output;
    }
	
	// dead are red
	if( pCell->phenotype.death.dead == false )
	{
		output[0] = "red";
		output[2] = "red";
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

double dot_product_ext(const std::vector<double>& v, const std::vector<double>& w)
{
    // Initialize the result of the dot product to zero
	double out = 0.0; 
    
    // Ensure both vectors are of the same size
    if (v.size() != w.size()) {
        std::cerr << "Error: Vectors must be of the same size." << std::endl;
        return 0.0; // Return zero or handle the error as needed
    }

    // Compute the dot product of vectors v and w
    for (unsigned int i = 0; i < v.size(); i++)
    {
        out += (v[i] * w[i]);
    }

    // If the absolute value of the result is very small, set it to zero
    if (fabs(out) < 1e-10)
    {
        out = 0.0;
    }

    return out; // Return the computed dot product
}

double sign_function (double number)
{
	// double sign = 0.0
	if (number<0)
	{ return -1.0;}

	else
	{ return 1.0;}
}

void cell_ecm_adhesion_speed(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Function that computes the cell-ECM adhesion speed and stores it as motility speed (migration speed)
	
    // Check if the cell is dead
    if (phenotype.death.dead == true)
	{
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
        std::cout << "Cell is dead" << std::endl;
        return;
	}

    // Determine which voxel to use based on the specified method
    int voxel_index = 0;

	if (parameters.strings("nearest_voxel_speed") == "front")
	{
        // Pick voxel closest to the cell's membrane in the direction of movement
        std::vector<double> direction = {pCell->custom_data["total_velocity_x"], 
                                         pCell->custom_data["total_velocity_y"], 
                                         pCell->custom_data["total_velocity_z"]};
        normalize(&direction); // Normalize the direction vector
		std::vector<double> scaled_direction = pCell->phenotype.geometry.radius * direction;

        // Calculate the position of the membrane in the direction of the cell
        std::vector<double> position_membrane = pCell->position;
        for (size_t i = 0; i < position_membrane.size(); ++i) {
            position_membrane[i] += scaled_direction[i];
        }

        // Compute the index of the voxel at the position of the membrane
        voxel_index = microenvironment.nearest_voxel_index(position_membrane);
	}
	else if (parameters.strings("nearest_voxel_speed") == "position")
	{
        // Compute the index of the voxel at the cell's current position
        voxel_index = microenvironment.nearest_voxel_index(pCell->position);
	}
	else
	{
        // Handle case where the voxel selection method is not defined
        std::cout << "Nearest voxel for computing migration speed not defined!" << std::endl;
        abort(); // Terminate the program if voxel method is not defined
		return;
	}

    // Retrieve ECM density at the computed voxel index
	double ecm_density = ecm.ecm_voxels[voxel_index].density;

    // If the motility vector is to be normalized, set migration speed to 1.0
    if (parameters.bools("normalize_ecm_influenced_motility_vector") == true)
	{
		pCell->phenotype.motility.migration_speed = 1.0;
	}
	else
	{
        // Otherwise, set migration speed based on the magnitude of the motility vector
        pCell->phenotype.motility.migration_speed = norm(phenotype.motility.migration_bias_direction);
        // std::cout << "Magnitude of motility vector is " << pCell->phenotype.motility.migration_speed << std::endl;
	}
	
    // Set migration bias to 1.0 to ensure no additional random motion is added
	phenotype.motility.migration_bias = 1.0; 

    // Retrieve ribose concentration from parameters
	double ribose_concentration = parameters.doubles("ribose_concentration");
	
    // Get the maximum migration speed at ribose concentration of 0 (S_0)
    double max_migration_speed = get_single_base_behavior(pCell, "migration speed");

    // Calculate maximum migration speed (S_rib) as a function of ribose concentration
    double sigma = parameters.doubles("sigma");
    max_migration_speed *= exp(-sigma * ribose_concentration);
		
    // Compute cell-ECM adhesion speed, which is stored in migration speed
	pCell->phenotype.motility.migration_speed *= 4 * max_migration_speed * ecm_density;

	return; 
}


void cell_ecm_adhesion_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Function to compute cell-ECM adhesion direction, stored in motility direction

	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		std::cout<<"Cell is dead"<<std::endl;
				
		return; 
	}

	// Find nearest voxel to cell's front or position
    int voxel_index = 0;
    if (parameters.strings("nearest_voxel_direction") == "front")
    {
        //Pick voxel closest to cell's membrane in direction of movement
        std::vector<double> direction = {pCell->custom_data["total_velocity_x"], pCell->custom_data["total_velocity_y"], pCell->custom_data["total_velocity_z"]};
		normalize(&direction); // Normalize the direction vector
        std::vector<double> scaled_direction = pCell->phenotype.geometry.radius * direction;

        // Calculate the position of the membrane in the direction of the cell
        std::vector<double> position_membrane = pCell->position;
        for (size_t i = 0; i < position_membrane.size(); ++i) {
            position_membrane[i] += scaled_direction[i];
        }

        // Computing the index of the voxel at that position
        voxel_index = microenvironment.nearest_voxel_index( position_membrane );
    }
    
    else if (parameters.strings("nearest_voxel_direction") == "position")
    {
        // Computing the index of the voxel at cell position
        voxel_index = microenvironment.nearest_voxel_index( pCell->position );
    }

    else
    {
        std::cout<<"Nearest voxel for computing motility direction not defined!"<<std::endl;
        abort();
        return;
    }

	// Get random vector 
    std::vector<double> d_random(3,0.0);
    if( phenotype.motility.restrict_to_2D == true )
    { 
        d_random = UniformOnUnitCircle(); 
    }
    else
    { 
        d_random = UniformOnUnitSphere(); 
    }

	// Initiate preferred direction 
	std::vector<double> d_pref;

	// Preferred direction given by random vector with no chemotaxis
	double chemotaxis_bias = pCell->custom_data["chemotaxis_bias"];

	if(chemotaxis_bias == 0)
	{
		d_pref = d_random;
	}

	// Preferred direction with chemotaxis
	else
	{
		// Get vector for chemotaxis 
		static int substrate_index = microenvironment.find_density_index( "substrate" ); 
		std::vector<double> chemotaxis_grad = pCell->nearest_gradient(substrate_index);
		normalize( &chemotaxis_grad ); 

		// combine cell chosen random direction and chemotaxis direction (like standard update_motility function)
		d_pref = (1 - chemotaxis_bias) * d_random + chemotaxis_bias * chemotaxis_grad;
	}

	normalize( &d_pref ); 

	// Compute fiber orientation and anisotropy contribution to motility preferred direction
	std::vector<double> motility_direction;

    /************************* cell-ECM adhesion direction with density ************************/
	if (parameters.strings("ecm_definition") == "ecm_density")
	{
		motility_direction = d_pref;
	}

	/************************* cell-ECM adhesion direction with fibres ************************/
	else if(parameters.strings("ecm_definition") == "ecm_fibers")
	{
		// Get cell sensitivity for ECM 
		double ecm_sensitivity = pCell->custom_data["ecm_sensitivity"];

		double anisotropy = ecm.ecm_voxels[voxel_index].anisotropy; 
		std::vector<double> fiber_orientation = ecm.ecm_voxels[voxel_index].ecm_fiber_alignment;

		// to determine direction along fiber_orientation, find part of d_pref that is perpendicular to fiber_orientation; 
		std::vector<double> d_perp = d_pref - dot_product_ext(d_pref, fiber_orientation) * fiber_orientation; 
		normalize( &d_perp ); 
		
		// find constants to span d_pref with d_perp and fiber_orientation
		double c_1 = dot_product_ext( d_pref , d_perp ); 
		double c_2 = dot_product_ext( d_pref, fiber_orientation ); 

		// calculate bias away from directed motility - combination of sensitity to ECM and anisotropy
		double gamma = ecm_sensitivity * anisotropy; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.
		
		// Compute motility direction
		motility_direction = (1.0 - gamma) * c_1 * d_perp + c_2 * fiber_orientation;
	}
	else
	{
		std::cout<<"ECM definition no declared!!!"<<std::endl;
		abort();
		return;
	}

	pCell->phenotype.motility.migration_bias_direction = motility_direction;

	return;
}


void cell_ecm_repulsion(Cell* pCell , Phenotype& phenotype , double dt) 
{ 
    // Function that computes the cell-ECM repulsion velocity

    // Get cell's direction of movement
    std::vector<double> direction = pCell->velocity;
    normalize(&direction); // Normalize the direction vector

    double cell_radius = pCell->phenotype.geometry.radius;

    // Pick voxel closest to cell's front
    std::vector<double> scaled_direction = cell_radius * direction;
    
    // Calculate the position of the membrane in the direction of the cell
    std::vector<double> position_membrane = pCell->position;
    for (size_t i = 0; i < position_membrane.size(); ++i) {
        position_membrane[i] += scaled_direction[i];
    }

    // Compute the index of the voxel at that position
    int voxel_index_memb = microenvironment.nearest_voxel_index(position_membrane);

    // Ensure voxel index is valid
    if (voxel_index_memb < 0 || voxel_index_memb >= ecm.ecm_voxels.size()) {
        std::cerr << "Invalid voxel index: " << voxel_index_memb << std::endl;
        return;
    }

    double ecm_density = ecm.ecm_voxels[voxel_index_memb].density;
    ecm_density = std::min(ecm_density, 1.0); 
    // std::cout<<"ecm_density: "<<ecm_density<<std::endl;

    // If ECM density is positive, compute repulsion
    if (ecm_density > 0.000001)
    {
        // Get total cell's velocity
        std::vector<double> velocity = pCell->velocity;

        // Assign resulting velocity to cell's velocity
        for (size_t i = 0; i < pCell->velocity.size(); ++i) {
            pCell->velocity[i] -= ecm_density * velocity[i];    
        }
    }

    return; 
}


void ecm_remodelling(Cell* pCell , Phenotype& phenotype , double dt) 
{
	/**************** Change in ECM density***************/

	// Get rate of ECM degradation 
	double r_density = pCell->custom_data["ecm_density_rate"]; 
	
	// Ribose concentration 
	double ribose_concentration = parameters.doubles("ribose_concentration");

	// Degradation rate (r_deg,rib) function wrt ribose concentration 
	double delta = parameters.doubles("delta");
	r_density *= exp(-delta * ribose_concentration);

	// Set density target 
	double density_target = 0;

	// Get nearest voxel to cell's position
	int voxel_index = microenvironment.nearest_voxel_index( pCell->position );

	//Pick voxel closest to cell's membrane in direction of movement
	std::vector<double> direction = {pCell->custom_data["total_velocity_x"], pCell->custom_data["total_velocity_y"], pCell->custom_data["total_velocity_z"]};

	std::vector<double> scaled_direction;
    normalize(&direction);
	scaled_direction = pCell->phenotype.geometry.radius * direction;

	// Calculating the position of the membrane in the direction of the cell
	std::vector<double> position_membrane = pCell->position + scaled_direction;

	if(parameters.strings("nearest_voxel_remodeling") == "front")
	{
		// Computing the index of the voxel at cell front
		voxel_index = microenvironment.nearest_voxel_index( position_membrane );
	}

    // Cell-ECM density interaction
    double ecm_density = ecm.ecm_voxels[voxel_index].density;
	
	ecm.ecm_voxels[voxel_index].density = ecm_density + r_density * dt  * (density_target - ecm_density);
	
	/************************* ECM remodeling with fibers ************************/
	if(parameters.strings("ecm_definition") == "ecm_fibers")
	{
		/**************** Cell-ECM Fiber realingment***************/

		// Computing the index of the voxel at cell position
		voxel_index = microenvironment.nearest_voxel_index( pCell->position );

		// Get index for accessing the ecm_fiber_alignment data structure and then copy the correct value
		std::vector<double> fiber_orientation = ecm.ecm_voxels[voxel_index].ecm_fiber_alignment; 
		double anisotropy = ecm.ecm_voxels[voxel_index].anisotropy;

		// Get cell migration speed
		double migration_speed = pCell->phenotype.motility.migration_speed;

		// Compute fibre realignment rate
		double r_f0 = pCell->custom_data["fiber_realignment_rate"];
		double r_fiber = r_f0 * migration_speed * (1 - anisotropy);

		double ddotf;
		std::vector<double> norm_cell_motility; // = pCell->velocity;
		norm_cell_motility.resize(3,0.0);
		norm_cell_motility = phenotype.motility.motility_vector;
		// std::cout<<"Motility vector: "<< norm_cell_motility<<std::endl;
		normalize(&norm_cell_motility);
		
		ddotf = dot_product_ext(fiber_orientation, norm_cell_motility);

		// flips the orientation vector so that it is aligned correctly with the moving cell for proper reoirentation later.
		fiber_orientation = sign_function(ddotf) * fiber_orientation; 

		std::vector<double> f_minus_d;
		f_minus_d.resize(3,0.0);

		for(int i = 0; i < 3; i++)
		{
			// if (ddotf<0.0)
			// {
			// 	fiber_orientation = -1.0 * fiber_orientation;
			// }
			f_minus_d[i] = fiber_orientation[i] - norm_cell_motility[i]; 
			ecm.ecm_voxels[voxel_index].ecm_fiber_alignment[i] -= dt * r_fiber * f_minus_d[i]; 
		}
		
		normalize(&(ecm.ecm_voxels[voxel_index].ecm_fiber_alignment)); 
		// if(PhysiCell_globals.current_time >600)
		// {
		// 	std::cout<<"f_minus_d: "<< f_minus_d<<std::endl;
		// 	std::cout<<"fiber orientation after: "<< ecm.ecm_voxels[voxel_index].ecm_fiber_alignment<<std::endl<<std::endl;
		// }
	
		/**************************** Cell-ECM Anisotropy Modification ******************/

		// Compute anisotropy update rate 
		double r_a0 = pCell->custom_data["anisotropy_increase_rate"]; // min-1
		
		double r_anisotropy = r_a0 * migration_speed;
		
		ecm.ecm_voxels[voxel_index].anisotropy = anisotropy + r_anisotropy * dt * (1- anisotropy);
	}
	
    return;
}

void proliferation_inhibition( Cell* pCell, Phenotype& phenotype, double dt )
{
	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		//std::cout<<2<<std::endl;
		std::cout<<"Cell is dead"<<std::endl;
	}

	// Proliferation depends on number of cells in contact
	double proliferation_rate = get_single_base_behavior(pCell,"cycle entry");

	// Get threshold neighbours parameter
	double overcrowding_threshold = pCell->custom_data["overcrowding_threshold"];

	// Get cell's neighbours 
	double n_attached =  pCell->state.neighbors.size();

    if(n_attached >= overcrowding_threshold)
	{
		proliferation_rate = 0;
	}

	// Computing the index of the voxel at cell position
	int voxel_index = microenvironment.nearest_voxel_index( pCell->position );

    // Dependence on ECM density
	proliferation_rate *= (1 - ecm.ecm_voxels[voxel_index].density);
	
	// Set new proliferation rate
	set_single_behavior(pCell,"cycle entry", proliferation_rate);
}


void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	// Update cell-ECM adhesion speed
	cell_ecm_adhesion_speed(pCell, phenotype, dt);

	// Assign speed to motility vector with same motility direction
	pCell->phenotype.motility.motility_vector = pCell->phenotype.motility.migration_bias_direction; 
	normalize( &(pCell->phenotype.motility.motility_vector) ); 

	pCell->phenotype.motility.motility_vector *= pCell->phenotype.motility.migration_speed;

	// Calling the standard update velocity of PhysiCell (updates cell-ECM adhesion direction and cell-cell adhesion and repulsion)
	standard_update_cell_velocity(pCell, phenotype, dt);

	// Update cell-ECM repulsion
	cell_ecm_repulsion(pCell, phenotype, dt);

	/********************************************SAVE TOTAL SPEED AND VELOCITY DATA*******************************************/
	// Compute total cell's speed and velocity
	double total_speed = sqrt(pow(pCell->velocity[0],2) + pow(pCell->velocity[1],2) + pow(pCell->velocity[2],2));
	// std::cout<<"total_speed: "<<total_speed<<std::endl;
	pCell->custom_data["total_speed"] = total_speed;
	pCell->custom_data["total_velocity_x"] = pCell->velocity[0];
	pCell->custom_data["total_velocity_y"] = pCell->velocity[1];
	pCell->custom_data["total_velocity_z"] = pCell->velocity[2];

    pCell->custom_data["attached_cells"] = pCell->state.neighbors.size();

	return; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 
	ecm_remodelling(pCell,phenotype,dt);

	proliferation_inhibition(pCell,phenotype,dt);
	
	return;
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void write_ecm_data_matlab( std::string filename )
{
    int number_of_data_entries = ecm.ecm_mesh.voxels.size();
    int size_of_each_datum = 8;

	// static int ecm_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	// static int ecm_density_index = microenvironment.find_density_index( "ECM" ); 

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "ECM_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )
    {
	    fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); // 1
        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp ); // 2
        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); //3
		fwrite( (char*) &( ecm.ecm_voxels[i].anisotropy), sizeof(double) , 1 , fp ); // 4
        fwrite( (char*) &( ecm.ecm_voxels[i].density), sizeof(double) , 1 , fp ); // 5
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[0]), sizeof(double) , 1 , fp ); // 6
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[1]), sizeof(double) , 1 , fp ); // 7
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[2]), sizeof(double) , 1 , fp ); // 8
		// This will only work if the diffusion and ECM meshes are the same size. Commenting out for actualrunning. To do a direct comparison, leave them in and change length. Will have to change the vizualization to get this from the regular BioFVM outputs.
		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][0]), sizeof(double) , 1 , fp ); // 9
		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][1]), sizeof(double) , 1 , fp ); // 10
		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][2]), sizeof(double) , 1 , fp ); // 11
    }

    fclose( fp );

    return;
}