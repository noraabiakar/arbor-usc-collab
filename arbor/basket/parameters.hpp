#include <iostream>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

#include <arbor/cable_cell.hpp>
#include <common/json_params.hpp>

std::vector<double> read_spike_times();

struct basket_params {
    std::vector<float> spikes;
    double temp, v_init;
    double tau1_syn, tau2_syn, e_syn;

    double gnatbar_ichan2, gkfbar_ichan2, gkabar_borgka, gncabar_nca, glcabar_lca;
    double gskbar_gskch, gkbar_cagk, gl_ichan2, catau_ccanl, caiinf_ccanl;

    double enat, ekf, ek, elca, esk, el_ichan2, cao;

    double ra, cm;

    double seg_res, weight;
    double dt, run_time;

    std::string morph_file;
};

basket_params read_params(int argc, char** argv) {
    basket_params p;

    using sup::param_from_json;

    if (argc<2) {
        throw std::runtime_error("No input parameter file provided.");
    }
    if (argc>2) {
        throw std::runtime_error("More than command line one option not permitted.");
    }

    std::string fname = argv[1];
    std::cout << "Loading parameters from file: " << fname << "\n";
    std::ifstream f(fname);

    if (!f.good()) {
        throw std::runtime_error("Unable to open input parameter file: "+fname);
    }

    nlohmann::json json;
    json << f;

    param_from_json(p.temp, "temp", json);
    param_from_json(p.v_init, "vinit", json);
    param_from_json(p.dt, "dt_arbor", json);
    param_from_json(p.run_time, "run_time", json);
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);

    param_from_json(p.spikes, "spikes", json);

    param_from_json(p.gnatbar_ichan2, "gnatbar_ichan2", json);
    param_from_json(p.gkfbar_ichan2, "gkfbar_ichan2", json);
    param_from_json(p.gkabar_borgka, "gkabar_borgka", json);
    param_from_json(p.gncabar_nca, "gncabar_nca", json);
    param_from_json(p.glcabar_lca, "glcabar_lca", json);
    param_from_json(p.gskbar_gskch, "gskbar_gskch", json);
    param_from_json(p.gkbar_cagk, "gkbar_cagk", json);
    param_from_json(p.gl_ichan2, "gl_ichan2", json);
    param_from_json(p.catau_ccanl, "catau_ccanl", json);
    param_from_json(p.caiinf_ccanl, "caiinf_ccanl", json);

    param_from_json(p.ra, "ra", json);
    param_from_json(p.cm, "cm", json);

    param_from_json(p.enat, "enat", json);
    param_from_json(p.ekf, "ekf", json);
    param_from_json(p.ek, "ek", json);
    param_from_json(p.elca, "elca", json);
    param_from_json(p.esk, "esk", json);
    param_from_json(p.el_ichan2, "el_ichan2", json);
    param_from_json(p.cao, "cao", json);

    param_from_json(p.seg_res, "seg_res", json);
    param_from_json(p.weight, "weight", json);
    param_from_json(p.morph_file, "morph_file", json);


    for (auto it=json.begin(); it!=json.end(); ++it) {
        std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
    }
    std::cout << "\n";

    return p;

}
