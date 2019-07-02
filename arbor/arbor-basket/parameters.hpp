#include <iostream>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

#include <arbor/cable_cell.hpp>
#include <common/json_params.hpp>

std::vector<double> read_spike_times();

struct basket_params {
    std::vector<double> spikes = {
            25.269724183039855, 29.37076391451496, 58.472477010286546, 93.80268485203328,
            112.71090127018375, 142.6472406502223, 293.3318516075217, 456.96763867081177,
            559.0808556112847, 1091.4908947982385, 1265.7627055896965, 1286.3817213526308,
            1726.173576102434, 1744.7268859340948, 2072.358562894649, 2429.33700123744,
            2438.882187257708, 2444.85687651729, 2500.3783411639783, 2523.5435646207175,
            2633.0843793058734, 2663.3213690478333, 3081.2091382891876, 3104.234316785872,
            3209.158889778191, 3311.7494479555003, 3628.0607334944084, 3892.4078916268645,
            3905.382779351989, 3972.478937325283, 4039.5190445966164, 4275.579471624872,
            4761.533462201488, 4875.327265653268, 4946.068519184462, 5186.38947000671,
            5250.193973512949, 5405.921064743746, 6075.548637890089, 6106.605889233615,
            6392.503123563068, 6484.87209757147, 6622.667183819736, 7132.979244248274,
            7214.140067854784, 7632.383314077037, 7662.664989661292, 7663.029726732657,
            8205.521919274834, 8514.66346930178, 8998.325190823954, 9387.223218469393,
            9453.933657798472, 9544.328220467469, 9858.711284584584, 9955.045230553718,
            9956.054906300105
    };
    double temp, v_init;
    double tau1_syn, tau2_syn, e_syn;

    double gnatbar_ichan2, gkfbar_ichan2, gkabar_borgka, gncabar_nca, glcabar_lca;
    double gskbar_gskch, gkbar_cagk, gl_ichan2, catau_ccanl, caiinf_ccanl;

    double enat, ekf, ek, elca, esk, el_ichan2, cao;

    double ra, cm;

    double seg_res, dt, weight;

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
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);

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
