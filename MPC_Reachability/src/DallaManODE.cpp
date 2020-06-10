//
// Created by Sriram Sankaranarayanan on 4/18/18.
//

#include "DallaManODE.h"
#include <boost/numeric/odeint.hpp>
#include <boost/assert.hpp>
#include <fstream>
#include <sstream>
#include <cmath>


using boost::numeric::odeint::integrate;

namespace odesimulation {
    bool debug = true;

    double get_parameter_from_map(const std::map<string, double> & patient_params, const char *what) {
        string s(what);
        auto vi = patient_params.find(s);
        BOOST_ASSERT(vi != patient_params.end());
        return vi->second;
    };

    struct dm_ode_int {
        double t0;
        std::map<double, double> &meal_ra;
        double current_insulin_U_per_hr;
        std::map<string, double> &patient_params;

        dm_ode_int(double start_time,
                   map<double, double> &meal_rate_appear,
                   map<string, double> &params, double ins) : t0(start_time),
                                                              meal_ra(meal_rate_appear),
                                                              current_insulin_U_per_hr(ins),
                                                              patient_params(params) {

        };

        double interpolate_meal_at_time(double t){
            double prev_time = 0;
            double prev_value = 0.0;
            for (auto it: meal_ra){
                if (prev_time <= t && t <= it.first ){
                    double cur_time = (double) it.first;
                    double cur_value = it.second;
                    if (cur_time > prev_time + 0.1)
                        return prev_value + (cur_value - prev_value) * (t - prev_time) / (cur_time - prev_time);
                } else {
                    prev_time = it.first;
                    prev_value = it.second;
                }
            }
            return 0.0;
        };

        double get_parameter(const char *what) {
            return get_parameter_from_map(patient_params, what);
        };

        void operator()(const std::vector<double> &x, std::vector<double> & dxdt, const double t) {

            // first get the meal at current time
            double meal_input = interpolate_meal_at_time(t0 + t);
            double IIR = 100.0 * current_insulin_U_per_hr / get_parameter("BW");
            double X = x[0];
            double Isc1 = x[1];
            double Isc2 = x[2];
            double Gt = x[3];
            double Gp = x[4];
            double Il = x[5];
            double Ip = x[6];
            double I1 = x[7];
            double Id = x[8];
            double Gs = x[9];

            double G = Gp / get_parameter("Vg");
            double I = Ip / get_parameter("Vi");

            double EGP = get_parameter("kp1") - get_parameter("kp2") * Gp - get_parameter("kp3") * Id;
            double S = 0;
            double insulin_ra = get_parameter("ka1") * Isc1 + get_parameter("ka2") * Isc2;
            double HE = get_parameter("HEb");
            double m3 = HE * get_parameter("m1") / (1.0 - HE);
            double Uii = get_parameter("Fsnc");
            double Vm = get_parameter("Vm0") + get_parameter("Vmx") * X;
            double Km = get_parameter("Km0");

            double Uid = Vm * Gt / (Km + Gt);

            double d_X = -get_parameter("p2u") * X + get_parameter("p2u") * (I - get_parameter("Ib"));
            double d_Isc1 = -(get_parameter("kd") + get_parameter("ka1")) * Isc1 + IIR;
            double d_Isc2 = get_parameter("kd") * Isc1 - get_parameter("ka2") * Isc2;
            double E = 0;
            if (Gp > get_parameter("ke2")) {
                E = get_parameter("ke1") * (Gp - get_parameter("ke2"));
            }

            double d_Gp = EGP + meal_input - Uii - E - get_parameter("k1") * Gp + get_parameter("k2") * Gt;
            double d_Gt = -Uid + get_parameter("k1") * Gp - get_parameter("k2") * Gt;
            double d_Gs = get_parameter("ksc") * (G - Gs);
            double d_Il = -(get_parameter("m1") + m3) * Il + get_parameter("m2") * Ip + S;
            double d_Ip = -(get_parameter("m2") + get_parameter("m4")) * Ip + get_parameter("m1") * Il + insulin_ra;

            double d_I1 = -get_parameter("ki") * (I1 - I);
            double d_Id = -get_parameter("ki") * (Id - I1);

            dxdt[0] = d_X;
            dxdt[1] = d_Isc1;
            dxdt[2] = d_Isc2;
            dxdt[3] = d_Gt;
            dxdt[4] = d_Gp;
            dxdt[5] = d_Il;
            dxdt[6] = d_Ip;
            dxdt[7] = d_I1;
            dxdt[8] = d_Id;
            dxdt[9] = d_Gs;



        }

    };


    struct dm_observer {
        std::map<double, double> & gp_values;
        double t0;
        double last_time;
        std::vector<double> & cur_state;
        double Vg;


        dm_observer(std::map<double, double> & gluc_val,
                double start_time,
                double param_vg,
                std::vector<double> & cur_st):
                gp_values(gluc_val),
                t0(start_time),
                last_time(0.0),
                cur_state(cur_st),
                Vg(param_vg){};

        void operator() (const std::vector<double> & x, double t){
            gp_values.insert(std::make_pair(t + t0, x[4]/Vg));

            for (int i = 0; i < 10; ++i) {
                cur_state[i] = x[i];
                if (debug){
                    std::cout << i << "," << x[i] << std::endl;
                }
            }
            last_time = t;

        };

//        const std::vector<double> & get_last_state() const { return cur_state; };

    };


    DallaManODE::DallaManODE(map<string, double> const &patient_params, map<double, double> const &meal_ra) :
            params(patient_params),
            _cur_state(10, 0.0),
            meal_rate_of_appearance(meal_ra),
            cur_time(0.0)
            {};

    double DallaManODE::time_step(int num_minutes, double insulin_U_per_hour) {
        // 1. Set up the integrator struct to integrate the ODE starting from current state.
        // 2. Set up an observer struct
        // 3. Integrate
        if (debug){
            std::cout << "STATE before the time step at : " << cur_time << std::endl;
            for (int i = 0; i < 10; ++i){
                std::cout << i << " --> " << _cur_state[i] << std::endl;

            }
        }
        dm_ode_int ode_instance(cur_time, meal_rate_of_appearance, params, insulin_U_per_hour);
        dm_observer ode_observe(gluc_val, cur_time, ode_instance.get_parameter("Vg"), _cur_state);

        integrate(ode_instance, _cur_state, 0.0, (double) num_minutes, 0.1, ode_observe);
        cur_time += num_minutes; // Advance Time


        if (debug){
            std::cout << "STATE after the time step at : " << cur_time << std::endl;
            for (int i = 0; i < 10; ++i){
                std::cout << i << " --> " << _cur_state[i] << std::endl;

            }
        }
        return _cur_state[4]/get_parameter_from_map(params, "Vg");
    }

    double DallaManODE::get_current_glucose() const {
        return _cur_state[4]/get_parameter_from_map(params,"Vg");
    }

    void DallaManODE::set_initial_state(double G0) {
        //patientParams.Gb = GStart;
        double Gb = G0;
        double Gpb = get_parameter_from_map(params, "Vg") * G0;
        double Gtb = (get_parameter_from_map(params, "k1") * Gpb
                      - get_parameter_from_map(params, "EGPb")
                      + get_parameter_from_map(params, "Fsnc"))/get_parameter_from_map(params, "k2");

        _cur_state[0] = 0.0;
        _cur_state[1] = get_parameter_from_map(params, "isc1ss");
        _cur_state[2] = get_parameter_from_map(params, "isc2ss");
        _cur_state[3] = Gtb;
        _cur_state[4] = Gpb;
        _cur_state[5] = get_parameter_from_map(params, "Ilb");
        _cur_state[6] = get_parameter_from_map(params, "Ipb");
        _cur_state[7] = get_parameter_from_map(params, "Ib");
        _cur_state[8] = get_parameter_from_map(params, "Ib");
        _cur_state[9] = Gb;


        return;
    }

    /* -- Read Patient Parameters -- */

    void read_parameters_from_file(const char *filename, map<string, double> &params) {
        std::ifstream fhandle(filename, std::ios_base::in);
        string param_name;
        double param_value;
        string line;
        while (std::getline(fhandle, line)){
            std::istringstream iss(line);
            iss >> param_name >> param_value;
            std::cout << ">"<< param_name << "=" << param_value << std::endl;
            params.insert(std::make_pair(param_name, param_value));
        }
        fhandle.close();
    }

    /*-- Simulate meals --*/
    struct meal_ode_model {
        const map<double, double> & meal_data;
        const map<string, double> & params;
        const double meal_duration = 15.0;

        meal_ode_model(const map<double, double> & meal_times_and_carbs,
                       const map<string, double> & patient_params): meal_data(meal_times_and_carbs),
                                                        params(patient_params){};


        double get_parameter(const char *what) {
            return get_parameter_from_map(params, what);
        };

        std::pair<double, double> get_relevant_meal(const double t){
            for(auto const p: meal_data){
                if (p.first <= t && t <= p.first + meal_duration){
                    return p;
                }
            }
            return std::make_pair(-1.0, -1.0);
        };


        void operator() (const std::vector<double> & x, std::vector<double> & dxdt, const double t){
            double relevant_meal_time, relevant_meal_carb;
            std::tie(relevant_meal_time, relevant_meal_carb) = get_relevant_meal(t);
            double D = 0.0, dImp = 0.0;
            if (relevant_meal_time >= 0.0){
                dImp = relevant_meal_carb/meal_duration;
                D = relevant_meal_carb * ( t - relevant_meal_time)/meal_duration;
            }
            double bb = get_parameter("b");
            double d = get_parameter("d");
            double kmin = get_parameter("kmin");
            double kmax = get_parameter("kmax");
            double kabs = get_parameter("kabs");

            double Qsto1 = x[0];
            double Qsto2 = x[1];
            double Qgut = x[2];

            double Qsto = Qsto1 + Qsto2;
            double alpha = 5.0/ (2.0 * (0.01+ relevant_meal_carb)* (1- bb));
            double beta =  5.0/ (2.0 * (0.01+ relevant_meal_carb)* d);

            double kempt =  kmin + (kmax - kmin)/2 * ( 2+  tanh( alpha * (Qsto - bb * D) ) - tanh( beta * (Qsto - d * D)));
            double d_Qsto1 = - kmax * Qsto1 + dImp;
            double d_Qsto2 = - kempt * Qsto2 + kmax * Qsto1;
            double d_Qgut = - kabs * Qgut + kempt * Qsto2;
            dxdt[0] = d_Qsto1;
            dxdt[1] = d_Qsto2;
            dxdt[2] = d_Qgut;


        };
    };

    struct meal_ra_observer {
        double f;
        double kabs;
        map<double, double> & result;
        meal_ra_observer(double f_param, double kabs_param, map<double, double> & meal_ra_map):f(f_param),
                                                                                               kabs(kabs_param),
                                                                                               result(meal_ra_map){};

        void operator() (const std::vector<double> & x, const double t){
            double ra = f * kabs * x[2];
            result.insert(std::make_pair(t, ra));
        }
    };

    void simulate_meals_to_obtain_rate_of_appearance(const map<double, double> & meal_time_and_carbs,
                                                     const map<string, double>  & patient_params,
                                                     map< double, double> & meal_ra){
        // 1. find how much time to simulate

        double max_time = (meal_time_and_carbs.rbegin()) -> first + 400.0;
        // 2. Next set up the ODE simulator
        meal_ode_model ode_model(meal_time_and_carbs, patient_params);
        // 3. Initial state
        std::vector<double> init_state ({0.0, 0.0, 0.0});
        // 4. Observer
        meal_ra_observer ode_observer( ode_model.get_parameter("f"), ode_model.get_parameter("kabs"), meal_ra);
        // 5. Simulate and get the result
        integrate(ode_model, init_state, 0.0, max_time, 0.1, ode_observer);
    }

};