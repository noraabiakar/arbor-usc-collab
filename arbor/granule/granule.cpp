/*
 * A miniapp that demonstrates how to make a ring model
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <arbor/assert_macro.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morphology.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/spike_source_cell.hpp>
#include <arbor/swcio.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

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
void write_trace_json(const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell granule_cell(std::string filename, cell_layers layer_info, const granule_params& params);

class granule_recipe: public arb::recipe {
public:
    granule_recipe(const granule_params& gparams, const cell_layers& layer_info):
        num_cells_(1),
        params_(gparams),
        layer_info_(layer_info),
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
        auto cell = granule_cell(params_.morph_file, layer_info_, params_);
        return std::move(cell);
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
        unsigned glob_syn_id = 0;
        for (auto layer: layer_info_.synapses.map) {
            if (layer.first == params_.syn_layer) {
                break;
            }
            glob_syn_id += layer.second.size();
        }
        glob_syn_id += params_.syn_id;

        std::vector<arb::event_generator> gens;
        arb::pse_vector svec;
        for (auto s: params_.spikes) {
            if (s > params_.run_time) {
                break;
            }
            svec.push_back({{0, glob_syn_id}, s, event_weight_});
        }
        gens.push_back(arb::explicit_generator(svec));
        return gens;
    }


    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid) const override {
        return 1;
    }

    arb::probe_info get_probe(cell_member_type id) const override {
        // Get the appropriate kind for measuring voltage.
        cell_probe_address::probe_kind kind = cell_probe_address::membrane_voltage;
        // Measure at the soma.
        arb::segment_location loc(0, 0.5);

        return arb::probe_info{id, kind, cell_probe_address{loc, kind}};
    }

    arb::util::any get_global_properties(cell_kind k) const override {
        return cell_gprop_;
    }

    void add_ion(const std::string& ion_name, int charge, double init_iconc, double init_econc, double init_revpot) {
        cell_gprop_.add_ion(ion_name, charge, init_iconc, init_econc, init_revpot);
    }

private:
    cell_size_type num_cells_;
    granule_params params_;
    cell_layers layer_info_;
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

        auto layer_info = get_layer_info(params.morph_file, params.seg_res);

        // Create an instance of our recipe.
        granule_recipe recipe(params, layer_info);
        recipe.add_ion("nca", 2, 1.0, 1.0, 0);
        recipe.add_ion("nat", 1, 1.0, 1.0, params.enat);
        recipe.add_ion("kf",  1, 1.0, 1.0, params.ekf);
        recipe.add_ion("ks",  1, 1.0, 1.0, params.eks);
        recipe.add_ion("lca", 2, 1.0, 1.0, params.elca);
        recipe.add_ion("tca", 2, 1.0, 1.0, params.etca);
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

arb::cable_cell granule_cell(
        std::string filename,
        cell_layers layer_info,
        const granule_params& params) {

    layer_to_vec_unsigned synapse_ids;

    std::ifstream f(filename);
    if (!f) throw std::runtime_error("unable to open file");

    auto samples = arb::parse_swc_file(f);

    for (auto& s: samples) {
        if (s.type == arb::swc_record::kind::soma) {
            s.r *= std::sqrt(2);
            break;
        }
    }

    auto morph = arb::swc_as_morphology(samples);
    arb::cable_cell cell = arb::make_cable_cell(morph);

    int tot_comp = 1;
    for (auto& segment: cell.segments()) {
        if (!segment->as_soma()) {
            auto length = segment->as_cable()->length();
            auto n = (unsigned) std::ceil(length / params.seg_res);
            segment->set_compartments(n);
            tot_comp += n;
        }
    }

    std::cout << "total segments: " << tot_comp << std::endl;

    int tot_synapses = 0;
    for (auto layer: layer_info.synapses.map) {
        int c = 0;
        for (auto v: layer.second) {
            if(params.syn_layer == layer.first && params.syn_id == c) {
                arb::mechanism_desc exp2syn("exp2syn");
                exp2syn["tau1"] = params.tau1_syn;
                exp2syn["tau2"] = params.tau2_syn;
                exp2syn["e"] = params.e_syn;
                cell.add_synapse({v.segment, v.pos}, exp2syn);
                std::cout << "selected synapse: " << v.segment << " " << v.pos << std::endl;
            }
            else {
                arb::mechanism_desc exp2syn("exp2syn");
                exp2syn["tau1"] = params.tau1_reg;
                exp2syn["tau2"] = params.tau2_reg;
                exp2syn["e"] = params.e_reg;
                cell.add_synapse({v.segment, v.pos}, exp2syn);
            }
            c++;
            tot_synapses++;
        }
    }
    std::cout << "total synapses: " << tot_synapses << std::endl;

    cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";

    unsigned seg_id = 0;
    for (auto& segment: cell.segments()) {
        segment->parameters.axial_resistivity = 410 * params.ra_mult;

        if(segment->as_soma()) {
            arb::mechanism_desc ichan2("ichan2");
            ichan2["gnatbar"] = 0.120    * params.gnatbar_ichan2;
            ichan2["gkfbar"]  = 0.016    * params.gkfbar_ichan2;
            ichan2["gksbar"]  = 0.006    * params.gksbar_ichan2;
            ichan2["gl"]      = 0.00004  * params.gl_ichan2;
            ichan2["el"]      =            params.el_ichan2;

            arb::mechanism_desc borgka("borgka");
            borgka["gkabar"]  = 0.001    * params.gkabar_borgka;

            arb::mechanism_desc nca("nca");
            nca["gncabar"]    = 0.001    * params.gncabar_nca;

            arb::mechanism_desc lca("lca");
            lca["glcabar"]    = 0.005    * params.glcabar_lca;

            arb::mechanism_desc cat("cat");
            cat["gcatbar"]    = 0.000037 * params.gcatbar_cat;

            arb::mechanism_desc gskch("gskch");
            gskch["gskbar"]   = 0.001    * params.gskbar_gskch;

            arb::mechanism_desc cagk("cagk");
            cagk["gkbar"]     = 0.0006   * params.gkbar_cagk;

            arb::mechanism_desc ccanl("ccanl");
            ccanl["catau"]    = 10       * params.catau_ccanl;
            ccanl["caiinf"]   = 5.0e-6   * params.caiinf_ccanl;

            segment->parameters.membrane_capacitance = 1.0 * params.cm_mult/100;

            segment->add_mechanism(ichan2);
            segment->add_mechanism(borgka);
            segment->add_mechanism(nca);
            segment->add_mechanism(lca);
            segment->add_mechanism(cat);
            segment->add_mechanism(gskch);
            segment->add_mechanism(cagk);
            segment->add_mechanism(ccanl);

        } else {
            arb::mechanism_desc ichan2("ichan2");
            arb::mechanism_desc nca("nca");
            arb::mechanism_desc lca("lca");
            arb::mechanism_desc cat("cat");
            arb::mechanism_desc gskch("gskch");
            arb::mechanism_desc cagk("cagk");
            arb::mechanism_desc ccanl("ccanl");

            if(std::binary_search(layer_info.segments.map["granuleCellLayer"].begin(),
                                  layer_info.segments.map["granuleCellLayer"].end(), seg_id)) {
                ichan2["gnatbar"] = 0.018   * params.gnatbar_ichan2;
                ichan2["gkfbar"]  = 0.004;
                ichan2["gksbar"]  = 0.006;
                ichan2["gl"]      = 0.00004 * params.gl_ichan2;
                ichan2["el"]      =           params.el_ichan2;
                nca["gncabar"]    = 0.003   * params.gncabar_nca;
                lca["glcabar"]    = 0.0075;
                cat["gcatbar"]    = 0.000075;
                gskch["gskbar"]   = 0.0004;
                cagk["gkbar"]     = 0.0006  * params.gkbar_cagk;
                ccanl["catau"]    = 10       * params.catau_ccanl;
                ccanl["caiinf"]   = 5.0e-6   * params.caiinf_ccanl;
                segment->parameters.membrane_capacitance = 1.0 * params.cm_mult/100;
            }

            if(std::binary_search(layer_info.segments.map["innerThird"].begin(),
                                  layer_info.segments.map["innerThird"].end(), seg_id)) {
                ichan2["gnatbar"] = 0.013    * params.gnatbar_ichan2;
                ichan2["gkfbar"]  = 0.004;
                ichan2["gksbar"]  = 0.006;
                ichan2["gl"]      = 0.000063 * params.gl_ichan2;
                ichan2["el"]      =            params.el_ichan2;
                nca["gncabar"]    = 0.001    * params.gncabar_nca;
                lca["glcabar"]    = 0.0075;
                cat["gcatbar"]    = 0.00025;
                gskch["gskbar"]   = 0.0002;
                cagk["gkbar"]     = 0.001    * params.gkbar_cagk;
                ccanl["catau"]    = 10       * params.catau_ccanl;
                ccanl["caiinf"]   = 5.0e-6   * params.caiinf_ccanl;
                segment->parameters.membrane_capacitance = 1.6 * params.cm_mult/100;
            }

            if(std::binary_search(layer_info.segments.map["middleThird"].begin(),
                                  layer_info.segments.map["middleThird"].end(), seg_id)) {
                ichan2["gnatbar"] = 0.008    * params.gnatbar_ichan2;
                ichan2["gkfbar"]  = 0.001;
                ichan2["gksbar"]  = 0.006;
                ichan2["gl"]      = 0.000063 * params.gl_ichan2;
                ichan2["el"]      =            params.el_ichan2;
                nca["gncabar"]    = 0.001    * params.gncabar_nca;
                lca["glcabar"]    = 0.0005;
                cat["gcatbar"]    = 0.0005;
                gskch["gskbar"]   = 0.0;
                cagk["gkbar"]     = 0.0024   * params.gkbar_cagk;
                ccanl["catau"]    = 10       * params.catau_ccanl;
                ccanl["caiinf"]   = 5.0e-6   * params.caiinf_ccanl;
                segment->parameters.membrane_capacitance = 1.6 * params.cm_mult/100;
            }

            if(std::binary_search(layer_info.segments.map["outerThird"].begin(),
                                  layer_info.segments.map["outerThird"].end(), seg_id)) {
                ichan2["gnatbar"] = 0.0      * params.gnatbar_ichan2;
                ichan2["gkfbar"]  = 0.001;
                ichan2["gksbar"]  = 0.008;
                ichan2["gl"]      = 0.000063 * params.gl_ichan2;
                ichan2["el"]      =            params.el_ichan2;
                nca["gncabar"]    = 0.001    * params.gncabar_nca;
                lca["glcabar"]    = 0.0;
                cat["gcatbar"]    = 0.001;
                gskch["gskbar"]   = 0.0;
                cagk["gkbar"]     = 0.0024   * params.gkbar_cagk;
                ccanl["catau"]    = 10       * params.catau_ccanl;
                ccanl["caiinf"]   = 5.0e-6   * params.caiinf_ccanl;
                segment->parameters.membrane_capacitance = 1.6 * params.cm_mult/100;
            }

            segment->add_mechanism(ichan2);
            segment->add_mechanism(nca);
            segment->add_mechanism(lca);
            segment->add_mechanism(cat);
            segment->add_mechanism(gskch);
            segment->add_mechanism(cagk);
            segment->add_mechanism(ccanl);
        }
        seg_id++;
    }

    return cell;
}

