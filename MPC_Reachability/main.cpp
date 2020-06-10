/*


LICENSE : Please see the license file, in the main directory

*/

#include "./src/propagate_intervals.h"
// #include "./include/sherlock.h"


using namespace std;
using namespace std::chrono;

int test_encoding(network_handler system_network,
	GRBModel * milp_model,
	GRBEnv * milp_env,
	vector< GRBVar > input_variables,
	GRBVar output_variable
);

datatype pick_random_data_point(
	datatype current_data_point,
	datatype physiological_constant,
	unsigned int seed
);
void pick_linearly_interpolated_data_points(
	datatype current_data_point,
	datatype taret_data_point,
	int no_of_data_points,
	vector< datatype >& data_points
);
int main(int argc, char ** argv)
{

	bool do_with_bounds = true;
	char * file_name_1 = (char *)"./network_files/AP_predict.nt"; // Prediction
	char * file_name_2 = (char *)"./network_files/AP_high.nt"; // Upper
	char * file_name_3 = (char *)"./network_files/AP_low.nt"; // Lower
	char * file_name_4 = (char *)"./network_files/AP_plant.nt"; // Simulation

	char * init_cond_file_name = (char *)"./network_files/input_1.csv";
    char * output_file = (char *)"testOut"; //base name for output

	datatype min_insulin_level = 0.008;
	datatype insulin_granularity = 0.001;
	datatype total_insulin_constraint = 5;

	bool print_detail = true;
	datatype upper_bound_network_error = 10.5;
	datatype lower_bound_network_error = 12.5;


	vector< datatype > target_val(3);
	target_val[0] = 140;
	target_val[1] = 220;
	target_val[2] = 100;

	unsigned int total_time_steps = 120;
	datatype look_ahead_period_for_MPC_MILP = 60; // 120

	assert((argc-1)%2 == 0);

	int arg_index = 1;
	int no_of_violations = 0;
	int ret_val;
	while( arg_index < argc )
	{
		if(strcmp(argv[arg_index], "-predict_network") == 0)
		{
			file_name_1 = argv[arg_index+1];
		}
        else if(strcmp(argv[arg_index], "-O")==0)
        {
            output_file=argv[arg_index+1];
        }
		else if(strcmp(argv[arg_index], "-upper_bound_network") == 0)
		{
			file_name_2 = argv[arg_index + 1];
		}
		else if(strcmp(argv[arg_index], "-lower_bound_network") == 0)
		{
			file_name_3 = argv[arg_index + 1];
		}
		else if(strcmp(argv[arg_index], "-sim_network") == 0)
		{
			file_name_4 = argv[arg_index + 1];
		}
		else if(strcmp(argv[arg_index], "-init_trace_file") == 0)
		{
			init_cond_file_name = argv[arg_index + 1];
		}
		else if(strcmp(argv[arg_index], "-mean_gluc_val") == 0)
		{
			sscanf(argv[arg_index + 1], "%lf", & target_val[0]);
		}
		else if(strcmp(argv[arg_index], "-upper_gluc_val") == 0)
		{
			sscanf(argv[arg_index + 1], "%lf", & target_val[1]);
		}
		else if(strcmp(argv[arg_index], "-lower_gluc_val") == 0)
		{
			sscanf(argv[arg_index + 1], "%lf", & target_val[2]);
		}
		else if(strcmp(argv[arg_index], "-total_time_steps") == 0)
		{
			sscanf(argv[arg_index + 1], "%d", & total_time_steps);
		}
		else if(strcmp(argv[arg_index], "-mpc_time") == 0)
		{
			sscanf(argv[arg_index + 1], "%lf", &look_ahead_period_for_MPC_MILP);
		}
		else if (strcmp(argv[arg_index], "-apply_bounds") == 0)
		{
			int some;
			sscanf(argv[arg_index + 1], "%d", & some);
			do_with_bounds = some?true:false;
		}
		else if (strcmp(argv[arg_index], "-min_insulin_level") == 0)
		{
			sscanf(argv[arg_index+1], "%lf", & min_insulin_level);
		}
		else if (strcmp(argv[arg_index], "-insulin_granularity") == 0)
		{
			sscanf(argv[arg_index+1], "%lf", & insulin_granularity);
		}
		else if (strcmp(argv[arg_index], "-total_insulin") == 0)
		{
			sscanf(argv[arg_index+1], "%lf", & total_insulin_constraint);
		}
		else if (strcmp(argv[arg_index], "-verbosity") == 0)
		{
			int some;
			sscanf(argv[arg_index + 1], "%d", & some);
			print_detail = some?true:false;
		}
		else
		{
			cout << "Argument received doesn't make sense !! " << endl;
			cout << "For argument value = " << argv[arg_index] << endl;
			cout << "Exiting... " << endl;
			exit(0);
		}
		arg_index += 2;
	}

	// Checking if the command line arguments make sense :

	assert(ifstream(file_name_1));
	assert(ifstream(file_name_2));
	assert(ifstream(file_name_3));
	assert(ifstream(file_name_4));
	assert(ifstream(init_cond_file_name));
	assert((target_val[0] > 0) && (target_val[0] < 500)) ;
	assert((target_val[1] > 0) && (target_val[1] < 500)) ;
	assert((target_val[2] > 0) && (target_val[2] < 500)) ;
	assert((total_time_steps > 0) && (total_time_steps < 300)) ;
	assert((look_ahead_period_for_MPC_MILP > 0) && (look_ahead_period_for_MPC_MILP < 200));
	assert((min_insulin_level >= 0.0) && (min_insulin_level < 5.0)) ;
	assert((insulin_granularity > 0.0) && (insulin_granularity <= 1.0)) ;
	assert((total_insulin_constraint > 0.0) && (total_insulin_constraint <= 30.0)) ;

	if(print_detail)
	{
		cout <<"\n \n \nData  details ============>  \n " << endl;
		cout << "Files being used ---------- " << endl << endl;
		cout << "|--------- Simulation file = " << file_name_4 << endl;
		cout << "|--------- Prediction file  = " << file_name_1 << endl;
		cout << "|--------- Upper bound file = " << file_name_2 << endl;
		cout << "|--------- Lower bound file = " << file_name_3 << endl;
		cout << "|--------- Init trace file = " << init_cond_file_name << endl;
		cout << "|------------------------------------------------------------------------------------------- " << endl << endl;

		cout << "\nValues being used ---------- \n"  << endl;
		cout << "|--------- Target glucose value = " << target_val[0] << endl;
		cout << "|--------- Upper bound of glucose value = " << target_val[1] << endl;
		cout << "|--------- Lower bound of glucose value = " << target_val[2] << endl;
		cout << "|--------- Total time steps to run for = " << total_time_steps << endl;
		cout << "|--------- MPC Time = " << look_ahead_period_for_MPC_MILP << " mins " << endl;
		cout << "|--------- Minimum Insulin value set = " << min_insulin_level << endl;
		cout << "|--------- Granularity of the insulin input = " << insulin_granularity << endl;
		cout << "|--------- Total insulin constraint = " << total_insulin_constraint << endl;
		cout << "|----------------------------------------------------------------------------------------- " << endl << endl;
		// cout << file_name_3 << endl;
	}

    char glucose_file_name[50] = "gluc_trace_";
    char insulin_file_name[50] = "ins_trace_";
    char plot_file_name[50] = "simData_";
    char statistics_file_name[50] = "stats_";

    strcat(glucose_file_name, output_file);
    strcat(insulin_file_name, output_file);
    strcat(plot_file_name, output_file);
    strcat(statistics_file_name, output_file);



	int dimension = 2;
	unsigned int i, no_of_infeasibilities = 0;
	network_handler prediction_network(file_name_1);
	network_handler upper_bound_network(file_name_2);
	network_handler lower_bound_network(file_name_3);
	network_handler simulation_network(file_name_4);

	plotting_data system_plots(dimension);

	// GRBEnv * milp_env = new GRBEnv();
	// milp_env->set(GRB_IntParam_OutputFlag, 0);
	// GRBModel * milp_model = new GRBModel(* milp_env);
	//
	// int no_of_network_inputs = prediction_network.no_of_inputs;
	// vector< GRBVar > input_variables(no_of_network_inputs);
	// GRBVar output_variable;
	// prediction_network.return_GUROBI_handle_of_network(milp_model, milp_env, input_variables, output_variable);
	// cout << "Reaches here... " << endl;
	// test_encoding(prediction_network, milp_model, milp_env, input_variables, output_variable);

	// Code for simulating the system

	// The sizes need to be adjusted for the correct values
	unsigned int model_step_index, control_step_index, time_index;



	// All times in the following setting is in minutes.
	datatype time_period_of_model_predictions = 30;
	datatype time_period_of_control_actions = 5;
	datatype physiological_constant = 10;
	clock_t begin, end;
	double current_time, average_time, max_time = 0, total_time = 0;


	unsigned int MPC_MILP_time = (unsigned int) (look_ahead_period_for_MPC_MILP / time_period_of_model_predictions);
	// 2 hours, of total look up and each prediction step is 30 mins out, thus 4.
	unsigned int look_ahead_steps = (unsigned int )(time_period_of_model_predictions / time_period_of_control_actions);
	// 30 mins of prediction step, with 5 mins for each of the control steps

	assert(prediction_network.no_of_inputs == upper_bound_network.no_of_inputs);
	assert(upper_bound_network.no_of_inputs == lower_bound_network.no_of_inputs);
	// assert(lower_bound_network.no_of_inputs == no_of_system_inputs);

	unsigned int no_of_glucose_inputs = (unsigned int) prediction_network.no_of_inputs/2;
	unsigned int no_of_insulin_inputs = (unsigned int) prediction_network.no_of_inputs/2;
	unsigned int no_of_system_inputs = no_of_glucose_inputs + no_of_insulin_inputs;


	vector< datatype > G_vector;
	vector< datatype > I_vector;
	vector< datatype > combined_vector;
	datatype new_G_value;
	datatype new_I_value;

	vector< vector< datatype > > trace_data;
	vector< datatype > system_point(dimension);

	vector< datatype > system_trace_G;
	vector< datatype > system_trace_I;

	vector< vector< datatype > > plot_data;
	vector< datatype > plot_point;
	datatype upper_G_value, lower_G_value;

	// Initialize the system trace data by the trace you want to
	// start it with :

	read_glucose_and_insulin_from_file(init_cond_file_name, system_trace_G, system_trace_I);

	vector< datatype > insulin_inputs;



	assert(system_trace_G.size() == system_trace_I.size());
	unsigned int no_of_bootstrap_data_points = system_trace_G.size();

	combined_vector.clear();
	combined_vector = system_trace_G;
	combined_vector.insert(combined_vector.end(), system_trace_I.begin(), system_trace_I.end());




	new_G_value = simulation_network.get_network_output(combined_vector);
	if(print_detail)
	{
		cout << "Test values for initial trace = " << endl;
		cout << "New G value from simulation network = " << new_G_value << endl;
		cout << "New G value from prediction network = " <<
		prediction_network.get_network_output(combined_vector) << endl;
		cout << "Upper bound = " << upper_bound_network.get_network_output(combined_vector) << endl;
		cout << "Lower bound = " << lower_bound_network.get_network_output(combined_vector) << endl;
		cout << endl << endl;
	}
	// Populate the initial numbers randomly to fill in the space,
	// between the curent values and the prediction time step.

	bool pick_linear_interpolation = true;
	// cout << "New set of G values added = " << endl;

	if(!pick_linear_interpolation)
	{
		i = 0;
		while(i < look_ahead_steps-1)
		{
			system_trace_G.push_back(
				pick_random_data_point(system_trace_G.back(), physiological_constant, i)
			);
			if(print_detail)
			{
				cout << system_trace_G.back() << " ";
			}
			i++;
		}
	}
	else
	{
		vector< datatype > new_data_points;
		pick_linearly_interpolated_data_points(system_trace_G.back(), new_G_value, look_ahead_steps-1, new_data_points);
		system_trace_G.insert(system_trace_G.end(), new_data_points.begin(), new_data_points.end());
		if(print_detail)
		{
			cout << "Values picked by the simulator for the initial setting : " << endl;
			print_vector(new_data_points);
		}
	}
	if(print_detail)
	{
		cout << endl;
	}

	system_trace_G.push_back(new_G_value);


	if(print_detail)
	{
		cout << " Initial settings : " << endl;
		cout << " MPC MILP look ahead period in terms of no of prediction steps = " << MPC_MILP_time << endl;
		cout << " No of time steps in model prediction = " << look_ahead_steps << endl;
		cout << " No of elements in Glucose Trace = " << system_trace_G.size() << endl;
		cout << " No of elements in Insulin Trace = " << system_trace_I.size() << endl;
		// cout << " G value predicted = " << system_trace_G.back() << endl;
		// cout << " Elements in system trace G are " << endl;
		// print_vector(system_trace_G);
		//
		// cout << " Elements in combined_vector are " << endl;
		// print_vector(combined_vector);
		// exit(0);
		// cout << "Doing sanity checks for the functions : " << endl;
		// GRBEnv * milp_env = new GRBEnv();
		// milp_env->set(GRB_IntParam_OutputFlag, 0);
		// GRBModel * milp_model = new GRBModel(* milp_env);
		//
		// int no_of_network_inputs = prediction_network.no_of_inputs;
		// vector< GRBVar > input_variables;
		// GRBVar var;
		// i = 0;
	  // while(i < no_of_network_inputs)
	  // {
	  //   var = milp_model->addVar(-GRB_INFINITY,
	  //                           GRB_INFINITY,
	  //                           0.0,
	  //                           GRB_CONTINUOUS,
	  //                           "input_val");
		//
	  //   input_variables.push_back(var);
	  //   i++;
	  // }
		// GRBVar output_var;
		// output_var = milp_model->addVar(-GRB_INFINITY,
		// 												GRB_INFINITY,
		// 												0.0,
		// 												GRB_CONTINUOUS,
		// 												"output_val");
		//
		// vector< vector< vector< datatype > > > weights;
		// vector< vector< datatype > > biases;
		// prediction_network.return_network_information(weights, biases);
		//
		// cout << "Reaches here... 1 " << endl;
		// cout << "combined_vector size = " << combined_vector.size() << endl;
		//
		// GRBLinExpr expr;
		// datatype data;
		// i = 0;
		// while(i < combined_vector.size())
		// {
		// 	expr = 0;
		// 	data = 1.0;
		// 	expr.addTerms(& data, & input_variables[i], 1);
		// 	milp_model->addConstr(expr, GRB_EQUAL, combined_vector[i],  "input_equality_constraint");
		// 	i++;
		// }
		//
		// relate_input_output_variables_through_network(weights, biases, milp_model, milp_env, input_variables, output_var);
		// GRBLinExpr objective_expr;
	  // objective_expr = 0;
		//
		// milp_model->write("debug_network.lp");
		//
	  // milp_model->setObjective(objective_expr, GRB_MINIMIZE);
	  // milp_model->optimize();
		//
		// cout << "Output variable value = " << output_var.get(GRB_DoubleAttr_X) << endl;
		// cout << "Input variables are =  ";
		// i = 0;
		// while(i < no_of_network_inputs)
		// {
		// 	cout << input_variables[i].get(GRB_DoubleAttr_X) << "  " ;
		// 	i++;
		// }
		// cout << endl;
		// exit(0);

	}

	if(print_detail)
	{
		cout <<" ---------- Simulation starts ----------------- " << endl;
		cout << endl;
	}

	vector< network_handler > all_networks;
	all_networks.push_back(prediction_network);
	// all_networks.push_back(lower_bound_network);
	// all_networks.push_back(upper_bound_network);
	all_networks.push_back(upper_bound_network);
	all_networks.push_back(lower_bound_network);

	time_index = no_of_bootstrap_data_points - 1;
	while(time_index < total_time_steps)
	{
		if(print_detail)
		{
			cout << "Time index = " << time_index << endl;
		}
		// Grab the last #no_of_glucose_inputs from the system_trace_G starting from the current location
		// and moving back

		// cout << "G values grabbed start from " << time_index - no_of_glucose_inputs + 1 ;
		G_vector.clear();
		i = time_index - no_of_glucose_inputs + 1;
		while(i < (time_index + 1) )
		{
			G_vector.push_back(system_trace_G[i]);
			i++;
		}
		// cout << " to " << i - 1 << endl;

		// Similarly grab the last few inputs from the insulin trace data
		// cout << "I values grabbed start from " << time_index - no_of_insulin_inputs + 1 ;

		I_vector.clear();
		i = time_index - no_of_insulin_inputs + 1;
		while(i < (time_index + 1))
		{
			I_vector.push_back(system_trace_I[i]);
			i++;
		}
		// cout << " to " << i - 1 << endl;

		// Get the control inputs for all the steps in the model prediction horizon , BUT just grab the first
		// one

		begin = clock();

		if(!do_with_bounds)
		{
			ret_val = compute_control_action(prediction_network, G_vector, I_vector,  MPC_MILP_time , look_ahead_steps ,
				physiological_constant, insulin_inputs, target_val, total_insulin_constraint);
			if(!ret_val)
				{
					cout << "Breaks from the loop ... " << endl;
					break;
				}
		}
		else
		{
			ret_val = compute_control_action(all_networks, G_vector, I_vector,  MPC_MILP_time , look_ahead_steps ,
				physiological_constant, insulin_inputs, target_val, total_insulin_constraint);
			if(!ret_val)
				{
					cout << "Breaks from the loop ... " << endl;
					break;
				}
			if(ret_val < 0)
			{
				no_of_infeasibilities ++ ;
			}
		}

		end = clock();
		current_time = (double)((end - begin) / (double) CLOCKS_PER_SEC);
		total_time += current_time;
		if(current_time > max_time )
		{
			max_time = current_time;
		}

		new_I_value = adjust_insulin_input(insulin_inputs[0], min_insulin_level, insulin_granularity);

		if(print_detail)
		{
			cout << "Control action taken = " << new_I_value << endl;
		}

		// Add this insulin input into the system_trace_I
		system_trace_I.push_back(new_I_value);

		// Adjust the G_Vector by throwing away the first value in the time series data
		// and appending something to it

		G_vector.erase(G_vector.begin());
		G_vector.push_back(system_trace_G[time_index + 1]);
		I_vector.erase(I_vector.begin());
		I_vector.push_back(new_I_value);


		combined_vector.clear();
		combined_vector = G_vector;
		combined_vector.insert(combined_vector.end(), I_vector.begin(), I_vector.end());

		// Push the new control input (INSULIN) and the Glucose value into the same trace data
		new_G_value = simulation_network.get_network_output(combined_vector);
		upper_G_value = upper_bound_network.get_network_output(combined_vector) + upper_bound_network_error ;
		lower_G_value = lower_bound_network.get_network_output(combined_vector) - lower_bound_network_error ;

		system_trace_G.push_back(new_G_value);

		plot_point.clear();
		plot_point.push_back(upper_G_value);
		plot_point.push_back(new_G_value);
		plot_point.push_back(lower_G_value);
		plot_point.push_back(new_I_value);

		plot_data.push_back(plot_point);

		if(
					(target_val[1] < upper_G_value ) ||
		  		(70 > lower_G_value)
		  )
		{
			if(print_detail)
			{
				cout << "Predicted value = " << new_G_value << endl;
				cout << "Upper bound value = " << upper_G_value << endl;
				cout << "Lower bound value = " << lower_G_value << endl;
				cout << "Bounds violated........................ " << endl;
			}
			no_of_violations++;
		}

		if((upper_G_value < new_G_value) || (lower_G_value > new_G_value))
		{
			if(print_detail)
			{
				cout << "Networks trained for upper and lower limits are misbehaving " << endl;
				cout << "Predicted value = " << new_G_value << endl;
				cout << "Upper bound value = " << upper_G_value << endl;
				cout << "Lower bound value = " << lower_G_value << endl;
			}
		}


		if(print_detail)
		{
			cout << "New G value added = " << new_G_value << endl;
		}

		time_index++;

	}

	if(print_detail)
	{
		cout << " ---------- Simulation ends -----------------  " << endl;
	}

	while(system_trace_I.size() < system_trace_G.size())
	{
		system_trace_I.push_back(0.0);
	}
	assert(system_trace_G.size() == system_trace_I.size());


	// For plotting purpose
	trace_data.clear();
	i = 0;
	while(i < system_trace_G.size())
	{
		system_point.clear();
		system_point.push_back(system_trace_G[i]);
		system_point.push_back(system_trace_I[i]);
		trace_data.push_back(system_point);
		i++;
	}


	vector< datatype > count_inf;
	count_inf.push_back(no_of_infeasibilities);

	system_plots.add_system_trace(trace_data);

	system_plots.plot(1);

	save_1D_vector_to_file(system_trace_G, glucose_file_name);
	save_1D_vector_to_file(system_trace_I, insulin_file_name);
	// save_1D_vector_to_file(count_inf, infeasibility_count);
	print_pancreas_data(plot_data, plot_file_name);

	average_time = total_time / (total_time_steps - no_of_bootstrap_data_points) ;
	ofstream file;
	file.open(statistics_file_name);
	file << "Average_time = " << average_time << "\n";
	file << "Max_time = " << max_time << "\n";
	file << "No_of_infeasibilities = " << no_of_infeasibilities << "\n";
	file.close();

	if(print_detail)
	{
		cout << "No of violations = " << no_of_violations << endl;
		cout << "\n \n" ;
		cout << "Glucose value written to --------------- " << glucose_file_name << endl;
		cout << "Insulin value written to --------------- " << insulin_file_name << endl;
		cout << "All plotting information written to ---- " << plot_file_name << endl;
		cout << "Running statistics written to ---------- " << statistics_file_name << endl;
		cout << endl << endl;

	}
	return 0;
}

int test_encoding(network_handler system_network,
	GRBModel * milp_model,
	GRBEnv * milp_env,
	vector< GRBVar > input_variables,
	GRBVar output_variable
)
{
	vector< vector< double > > input_region_bounds(system_network.no_of_inputs);
	vector< vector< double > > input_region_constraints;

	unsigned int i, j, no_of_points;
	i = 0;
	while(i < system_network.no_of_inputs)
	{
		input_region_bounds[i].push_back(-1.0);
		input_region_bounds[i].push_back(1.0);
		i++;
	}
	create_constraint_from_interval(input_region_constraints, input_region_bounds);
	vector< vector< vector < double > > > weights;
	vector< vector< double > > biases;
	system_network.return_network_information(weights, biases);

	double data;
	GRBLinExpr expr_buffer_0(0.0);
	GRBLinExpr expr_buffer_1(0.0);
	GRBLinExpr expr_buffer_2(0.0);
	GRBLinExpr expr_buffer_3(0.0);
	string const_name = "constant";
	GRBVar const_var = milp_model->addVar(1.0, 1.0, 0.0, GRB_CONTINUOUS, const_name);

	// Putting the constraints imposed by the input region constraints
	i = 0;
	while(i < input_region_constraints.size())
	{
		expr_buffer_0 = 0.0;
		j = 0;
		while(j < system_network.no_of_inputs)
		{
			data = input_region_constraints[i][j];
			expr_buffer_0.addTerms(& data, & input_variables[j], 1);
			j++;
		}
		data = input_region_constraints[i][j];
		expr_buffer_0.addTerms(& data, & const_var, 1);

		milp_model->addConstr(expr_buffer_0, GRB_GREATER_EQUAL, 0.0, "some_constraint_name");
		i++;
	}

	GRBLinExpr objective_expr;

	objective_expr = 0;
	data = 1;
	objective_expr.addTerms( & data, & output_variable, 1);
	milp_model->setObjective(objective_expr, GRB_MAXIMIZE);
	milp_model->optimize();

	cout << "Milp maxima " << milp_model->get(GRB_DoubleAttr_ObjVal) << endl;



	no_of_points = 1000;
	double max_val = -1, current_val;
	vector< vector< unsigned int > > active_weights ;
	vector< double > input_point;
	i = 0;
	while(i < no_of_points)
	{
		find_random_sample_with_seed(input_region_constraints, input_point, i);
		current_val = compute_network_output(input_point, weights, biases, active_weights);
		if(max_val < current_val)
		{
			max_val = current_val;
		}
		i++;
	}
	cout << "Max val computed experimentally = " << max_val << endl;

	return 0;
}

datatype pick_random_data_point(
	datatype current_data_point,
	datatype physiological_constant,
	unsigned int seed
)
{
	datatype return_data_point;
	srand(seed);
	datatype change = ((datatype)(rand() % 10)) / ((datatype)10.0) ;
	change *= 2 * (physiological_constant/2.0) ;
	return_data_point = current_data_point - (physiological_constant/2.0) + change ;

	// return_data_point = current_data_point;
	return return_data_point;
}
void pick_linearly_interpolated_data_points(
	datatype current_data_point,
	datatype target_data_point,
	int no_of_data_points,
	vector< datatype >& data_points
)
{
	assert(no_of_data_points > 0);
	data_points.clear();
	datatype increments = (target_data_point - current_data_point)/(no_of_data_points + 1);
	int i = 1;
	while(i <= no_of_data_points)
	{
		data_points.push_back(current_data_point + increments * i);
		i++;
	}

}
