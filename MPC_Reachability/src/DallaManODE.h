//
// Created by Sriram Sankaranarayanan on 4/18/18.
//

#ifndef DALLA_MAN_MODEL_CPP_DALLAMANODE_H
#define DALLA_MAN_MODEL_CPP_DALLAMANODE_H

#include <vector>
#include <map>
#include <string>

using std::string;
using std::map;
using std::vector;

namespace odesimulation {

    /* Function: read_parameters_from_file
     *
     * Read the patient parameters needed for simulation from the provided txt file.
     * Note: Sriram will generate the files for you corresponding to the patients we simulated for.
     *
     * Parameters:
     *    filename: name of the file with full path
     *    params: the result structure where we store the parameters
     */
    void read_parameters_from_file(const char * filename, map<string, double> & params);

    /*
     * Function: simulate_meals_to_obtain_rate_of_appearance
     *
     * Given a set of meal times and carbohydrate amounts in grams, this function
     * runs the meal gut absorption model to obtain a map of how much carbohydrates are released
     * into the blood, and when.
     *
     * Parameters:
     *    meal_time_and_carbs: a map from meal times -> carb amounts. The meal times must be in
     *    ascending order.
     *    patient_params: obtained by reading them from file. See read_parameters_file function.
     *    meal_ra: Result structure that can be passed directly into the Dalla man ODE model.
     */
    void simulate_meals_to_obtain_rate_of_appearance(const map<double, double> & meal_time_and_carbs,
                                                     const map<string, double>  & patient_params,
                                                     map< double, double> & meal_ra);

    class DallaManODE {
    protected:
        map<string, double> params; // This should be loaded in from matlab
        std::vector<double> _cur_state;
        map<double, double> meal_rate_of_appearance;
        map<double, double> insulin_input;
        map<double, double> gluc_val;
        double cur_time;
    public:
        /*
         * Function: DallaManODE constructor
         *
         * Constructs the object given patient parameters and the meal gut absorption results.
         *
         * Parameters:
         *   patient_params: the patient parameters we read in from the file.
         *   meal_ra: the meal rate of appearance.
         *
         * See Also:
         *   read_parameters_from_file
         *   simulate_meals_to_obtain_rate_of_appearance.
         */

        DallaManODE(map<string, double> const &patient_params, map<double, double> const &meal_ra);

        /*
         * Function: time_step
         *
         * Advance the time step of the simulation by specified number of minutes while
         * providing insulin at the specified rate in that time period.
         *
         * Parameters:
         *   num_minutes: how much to advance time by.
         *   insulin_U_per_hour: the insulin infused measured in U/hr. Caution: in the model the insulin
         *   may be specified as U delivered over 5 minutes. This must be multiplied by 12 to obtain U/hr.
         */

        double time_step(int num_minutes, double insulin_U_per_hour);

        /*
         * Function: get_current_glucose
         *
         * Gets the current value of BG levels.
         */
        double get_current_glucose() const;

        /*  Function: set_initial_state
         *
         *  Initialize the simulation by setting the initial state.
         */
        void set_initial_state(double G0);

        /*
         * Function: get_glucose_values
         *
         * Get a map of glucose values over time.
         * 
         */
        map<double, double> & get_glucose_values(){
            return gluc_val;
        };

    };

}


#endif //DALLA_MAN_MODEL_CPP_DALLAMANODE_H
