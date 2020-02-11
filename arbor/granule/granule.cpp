#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morph/locset.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/sample_tree.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/swcio.hpp>
#include <arbor/simulation.hpp>
#include <arbor/simple_sampler.hpp>

#include <arborenv/concurrency.hpp>
#include <arborenv/gpu_env.hpp>

#include "parameters.hpp"

#include <common/json_params.hpp>

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Writes voltage trace as a json file.
arb::sample_tree read_swc(const std::string& path);
void write_trace_json(const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell granule_cell(arb::sample_tree tree, const granule_params& params, cell_gid_type gid);

class granule_recipe: public arb::recipe {
public:
    granule_recipe(const granule_params& gparams):
        num_cells_(1),
        params_(gparams),
        event_weight_(params_.weight),
        catalogue_(arb::global_default_catalogue()) {

        cell_gprop_.catalogue = &catalogue_;

        cell_gprop_.default_parameters = arb::neuron_parameter_defaults;

        cell_gprop_.default_parameters.temperature_K = params_.temp + 273.15;
        cell_gprop_.default_parameters.init_membrane_potential = params_.v_init;

        cell_gprop_.default_parameters.ion_data["k"].init_reversal_potential = params_.ek;
        cell_gprop_.default_parameters.ion_data["ca"].init_ext_concentration = params_.cao;
    }

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        tree_ = read_swc(params_.morph_file);
        return granule_cell(tree_, params_, gid);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 1;
    }

    // The cell has one target synapse, which will be connected to cell gid-1.
    cell_size_type num_targets(cell_gid_type gid) const override {
        return 1;
    }

    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        arb::pse_vector svec;
        for (auto s: params_.spikes) {
            if (s > params_.run_time) {
                break;
            }
            // Add spikes to the first (and only) synapse on the cell
            svec.push_back({{0, 0}, s, event_weight_});
        }
        return {arb::explicit_generator(svec)};
    }


    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid) const override {
        return 1;
    }

    arb::probe_info get_probe(cell_member_type id) const override {
        arb::mlocation mid_soma = {0, 0.5};
        arb::cell_probe_address probe = {mid_soma, arb::cell_probe_address::membrane_voltage};
        return {id, 0, probe};
    }

    arb::util::any get_global_properties(cell_kind k) const override {
        return cell_gprop_;
    }

    void add_ion(const std::string& ion_name, int charge, double init_iconc, double init_econc, double init_revpot) {
        cell_gprop_.add_ion(ion_name, charge, init_iconc, init_econc, init_revpot);
    }

private:
    cell_size_type num_cells_;
    mutable arb::sample_tree tree_;
    granule_params params_;
    float event_weight_;

    arb::cable_cell_global_properties cell_gprop_;
    arb::mechanism_catalogue catalogue_;
};


int main(int argc, char** argv) {
    try {
        bool root = true;

        arb::proc_allocation resources;
        if (auto nt = arbenv::get_env_num_threads()) {
            resources.num_threads = nt;
        }
        else {
            resources.num_threads = arbenv::thread_concurrency();
        }

#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
        resources.gpu_id = arbenv::find_private_gpu(MPI_COMM_WORLD);
        auto context = arb::make_context(resources, MPI_COMM_WORLD);
        root = arb::rank(context) == 0;
#else
        resources.gpu_id = arbenv::default_gpu();
        auto context = arb::make_context(resources);
#endif

#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(context);
#endif

        // Print a banner with information about hardware configuration
        std::cout << "gpu:      " << (has_gpu(context)? "yes": "no") << "\n";
        std::cout << "threads:  " << num_threads(context) << "\n";
        std::cout << "mpi:      " << (has_mpi(context)? "yes": "no") << "\n";
        std::cout << "ranks:    " << num_ranks(context) << "\n" << std::endl;

        arb::profile::meter_manager meters;
        meters.start(context);

        auto params = read_params(argc, argv);

        // Create an instance of our recipe.
        granule_recipe recipe(params);
        recipe.add_ion("nca", 2, 1.0, 1.0, 0);
        recipe.add_ion("lca", 2, 1.0, 1.0, params.elca);
        recipe.add_ion("tca", 2, 1.0, 1.0, params.etca);
        recipe.add_ion("nat", 1, 1.0, 1.0, params.enat);
        recipe.add_ion("kf",  1, 1.0, 1.0, params.ekf);
        recipe.add_ion("ks",  1, 1.0, 1.0, params.eks);
        recipe.add_ion("sk",  1, 1.0, 1.0, params.esk);

        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};

        auto sched = arb::regular_schedule(params.dt);

        // This is where the voltage samples will be stored as (time, value) pairs
        arb::trace_data<double> voltage;

        // Now attach the sampler at probe_id, with sampling schedule sched, writing to voltage
        sim.add_sampler(arb::one_probe(probe_id), sched, arb::make_simple_sampler(voltage));

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                    [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                        recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                    });
        }

        meters.checkpoint("model-init", context);

        std::cout << "running simulation" << std::endl;

        sim.run(params.run_time, params.dt);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated at rate of "
                      << params.run_time << " ms between spikes\n";
            std::ofstream fid("spikes.gdf");
            if (!fid.good()) {
                std::cerr << "Warning: unable to open file spikes.gdf for spike output\n";
            }
            else {
                char linebuf[45];
                for (auto spike: recorded_spikes) {
                    auto n = std::snprintf(
                            linebuf, sizeof(linebuf), "%u %.4f\n",
                            unsigned{spike.source.gid}, float(spike.time));
                    fid.write(linebuf, n);
                }
            }
        }

        auto profile = arb::profile::profiler_summary();
        std::cout << profile << "\n";

        // Write the samples to a json file.
        if (root) write_trace_json(voltage);

        auto report = arb::profile::make_meter_report(meters, context);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in granule miniapp: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

void write_trace_json(const arb::trace_data<double>& trace) {
    std::string path = "./voltages.json";

    nlohmann::json json;
    json["name"] = "ring demo";
    json["units"] = "mV";
    json["cell"] = "0.0";
    json["probe"] = "0";

    auto& jt = json["data"]["time"];
    auto& jy = json["data"]["voltage"];

    for (const auto& sample: trace) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
    }

    std::ofstream file(path);
    file << std::setw(1) << json << "\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arb::cable_cell granule_cell(arb::sample_tree tree, const granule_params& params, cell_gid_type gid) {
    using namespace arb::reg;
    using namespace arb::ls;
    using mech = arb::mechanism_desc;
    using reg = arb::region;

    arb::label_dict dict;

    dict.set("soma", tagged(1));
    dict.set("dend", tagged(3));
    dict.set("granule_layer", intersect("dend", z_dist_from_root_gt(0)    , z_dist_from_root_le(6.68)));
    dict.set("inner_layer",   intersect("dend", z_dist_from_root_gt(6.68) , z_dist_from_root_le(20.04)));
    dict.set("middle_layer",  intersect("dend", z_dist_from_root_gt(20.04), z_dist_from_root_le(40.08)));
    dict.set("outer_layer",   intersect("dend", z_dist_from_root_gt(40.08), z_dist_from_root_le(66.9)));


    // Create the cell
    arb::morphology morpho(tree);
    arb::cable_cell cell(morpho, dict);


    // Draw 100 random points on the selected layer
    // Attach a synapse to one of the random locations
    mech exp2syn("exp2syn");
    exp2syn["tau1"] = params.tau1_syn;
    exp2syn["tau2"] = params.tau2_syn;
    exp2syn["e"] = params.e_syn;

    arb::mlocation selected;
    if (params.syn_layer == "soma") {
        auto soma_syn = thingify(uniform("soma", 0, 0, gid), cell.provider());
        selected = soma_syn.front();
    } else if (params.syn_layer == "granuleCellLayer") {
        auto granule_syns = thingify(uniform("granule_layer", 1,  100, gid), cell.provider());
        selected = granule_syns[params.syn_id%granule_syns.size()];
    } else if (params.syn_layer == "innerThird") {
        auto inner_syns  = thingify(uniform("inner_layer", 101, 200, gid), cell.provider());
        selected = inner_syns[params.syn_id%inner_syns.size()];
    } else if (params.syn_layer == "middleThird") {
        auto middle_syns  = thingify(uniform("middle_layer",201, 300, gid), cell.provider());
        selected = middle_syns[params.syn_id%middle_syns.size()];
    } else if (params.syn_layer == "outerThird") {
        auto outer_syns   = thingify(uniform("outer_layer", 301, 400, gid), cell.provider());
        selected = outer_syns[params.syn_id%outer_syns.size()];
    }
    std::cout << "selected synapse = {" << selected.branch << ", " << selected.pos << "}" << std::endl;
    cell.place(selected, exp2syn);


    // Add density mechanisms
    cell.default_parameters.discretization = arb::cv_policy_max_extent(params.seg_res, arb::cv_policy_flag::single_root_cv);
    cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";

    cell.paint("soma", mech("ichan2").set
                           ("gnatbar", 0.120*params.gnatbar_ichan2).set
                           ("gkfbar", 0.016*params.gkfbar_ichan2).set
                           ("gksbar", 0.006*params.gksbar_ichan2).set
                           ("gl", 0.00004*params.gl_ichan2).set
                           ("el",params.el_ichan2));
    cell.paint("soma", mech("borgka").set("gkabar", 0.001   *params.gkabar_borgka));
    cell.paint("soma", mech("nca").set("gncabar",   0.001   *params.gncabar_nca));
    cell.paint("soma", mech("lca").set("glcabar",   0.005   *params.glcabar_lca));
    cell.paint("soma", mech("cat").set("gcatbar",   0.000037*params.gcatbar_cat));
    cell.paint("soma", mech("gskch").set("gskbar",  0.001   *params.gskbar_gskch));
    cell.paint("soma", mech("cagk").set("gkbar",    0.0006  *params.gkbar_cagk));
    cell.paint("soma", mech("ccanl").set("catau",   10      *params.catau_ccanl).set("caiinf",5.0e-6*params.caiinf_ccanl));
    cell.paint("soma", arb::membrane_capacitance{1.0 * params.cm_mult/100});
    cell.paint("soma", arb::axial_resistivity{410 * params.ra_mult});

    cell.paint("granule_layer", mech("ichan2").set
                ("gnatbar", 0.018*params.gnatbar_ichan2).set
                ("gkfbar", 0.004).set
                ("gksbar", 0.006).set
                ("gl", 0.00004*params.gl_ichan2).set
                ("el",params.el_ichan2));
    cell.paint("granule_layer", mech("nca").set("gncabar",   0.003   *params.gncabar_nca));
    cell.paint("granule_layer", mech("lca").set("glcabar",   0.0075));
    cell.paint("granule_layer", mech("cat").set("gcatbar",   0.000075));
    cell.paint("granule_layer", mech("gskch").set("gskbar",  0.0004));
    cell.paint("granule_layer", mech("cagk").set("gkbar",    0.0006  *params.gkbar_cagk));
    cell.paint("granule_layer", mech("ccanl").set("catau",   10      *params.catau_ccanl).set("caiinf",5.0e-6*params.caiinf_ccanl));
    cell.paint("granule_layer", arb::membrane_capacitance{1.0 * params.cm_mult/100});
    cell.paint("granule_layer", arb::axial_resistivity{410 * params.ra_mult});

    cell.paint("inner_layer", mech("ichan2").set
            ("gnatbar", 0.013*params.gnatbar_ichan2).set
            ("gkfbar", 0.004).set
            ("gksbar", 0.006).set
            ("gl", 0.000063*params.gl_ichan2).set
            ("el",params.el_ichan2));
    cell.paint("inner_layer", mech("nca").set("gncabar",   0.001   *params.gncabar_nca));
    cell.paint("inner_layer", mech("lca").set("glcabar",   0.0075));
    cell.paint("inner_layer", mech("cat").set("gcatbar",   0.00025));
    cell.paint("inner_layer", mech("gskch").set("gskbar",  0.0002));
    cell.paint("inner_layer", mech("cagk").set("gkbar",    0.001   *params.gkbar_cagk));
    cell.paint("inner_layer", mech("ccanl").set("catau",   10      *params.catau_ccanl).set("caiinf",5.0e-6*params.caiinf_ccanl));
    cell.paint("inner_layer", arb::membrane_capacitance{1.6 * params.cm_mult/100});
    cell.paint("inner_layer", arb::axial_resistivity{410 * params.ra_mult});

    cell.paint("middle_layer", mech("ichan2").set
            ("gnatbar", 0.008*params.gnatbar_ichan2).set
            ("gkfbar", 0.001).set
            ("gksbar", 0.006).set
            ("gl", 0.000063*params.gl_ichan2).set
            ("el",params.el_ichan2));
    cell.paint("middle_layer", mech("nca").set("gncabar",   0.001   *params.gncabar_nca));
    cell.paint("middle_layer", mech("lca").set("glcabar",   0.0005));
    cell.paint("middle_layer", mech("cat").set("gcatbar",   0.0005));
    cell.paint("middle_layer", mech("gskch").set("gskbar",  0.0));
    cell.paint("middle_layer", mech("cagk").set("gkbar",    0.0024  *params.gkbar_cagk));
    cell.paint("middle_layer", mech("ccanl").set("catau",   10      *params.catau_ccanl).set("caiinf",5.0e-6*params.caiinf_ccanl));
    cell.paint("middle_layer", arb::membrane_capacitance{1.6 * params.cm_mult/100});
    cell.paint("middle_layer", arb::axial_resistivity{410 * params.ra_mult});

    cell.paint("outer_layer", mech("ichan2").set
            ("gnatbar", 0.0*params.gnatbar_ichan2).set
            ("gkfbar", 0.001).set
            ("gksbar", 0.008).set
            ("gl", 0.000063*params.gl_ichan2).set
            ("el",params.el_ichan2));
    cell.paint("outer_layer", mech("nca").set("gncabar",   0.001   *params.gncabar_nca));
    cell.paint("outer_layer", mech("lca").set("glcabar",   0.0));
    cell.paint("outer_layer", mech("cat").set("gcatbar",   0.001));
    cell.paint("outer_layer", mech("gskch").set("gskbar",  0.0));
    cell.paint("outer_layer", mech("cagk").set("gkbar",    0.0024  *params.gkbar_cagk));
    cell.paint("outer_layer", mech("ccanl").set("catau",   10      *params.catau_ccanl).set("caiinf",5.0e-6*params.caiinf_ccanl));
    cell.paint("outer_layer", arb::membrane_capacitance{1.6 * params.cm_mult/100});
    cell.paint("outer_layer", arb::axial_resistivity{410 * params.ra_mult});

    return cell;

}

arb::sample_tree read_swc(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("unable to open SWC file: "+path);

    return arb::swc_as_sample_tree(arb::parse_swc_file(f));
}

