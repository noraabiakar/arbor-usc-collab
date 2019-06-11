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
arb::cable_cell granule_cell(std::string filename, synapse_layers syn_layers, const granule_params& params);

class granule_recipe: public arb::recipe {
public:
    granule_recipe(const granule_params& gparams,
                   const synapse_layers& syn_layers):
        num_cells_(1), params_(gparams), syn_layers_(syn_layers), event_weight_(params_.weight) {}

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        auto cell = granule_cell(params_.morph_file, syn_layers_, params_);
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
        for (auto layer: syn_layers_.map) {
            if (layer.first == params_.syn_layer) {
                break;
            }
            glob_syn_id += layer.second.size();
        }
        glob_syn_id += params_.syn_id;

        std::vector<arb::event_generator> gens;
        arb::pse_vector svec;
        for (auto s: params_.spikes) {
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
        arb::cable_cell_global_properties a;
        a.temperature_K = params_.temp + 273.15;
        a.init_membrane_potential_mV = params_.v_init;
        return a;
    }

private:
    cell_size_type num_cells_;
    granule_params params_;
    synapse_layers syn_layers_;
    float event_weight_;
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

        auto syn_pos = get_synapse_positions(params.morph_file, params.seg_res);
        //auto syn_ids = get_synapse_ids(syn_pos);

        // Create an instance of our recipe.
        granule_recipe recipe(params, syn_pos);

        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};
        // The schedule for sampling is 10 samples every 1 ms.
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
        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(200, params.dt);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated at rate of "
                      << 200 << " ms between spikes\n";
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
        synapse_layers syn_layers,
        const granule_params& params) {

    unsigned_layers synapse_ids;

    std::ifstream f(filename);
    if (!f) throw std::runtime_error("unable to open file");

    auto morph = arb::swc_as_morphology(arb::parse_swc_file(f));
    arb::cable_cell cell = arb::make_cable_cell(morph);

    int tot_seg = 1;
    for (auto& segment: cell.segments()) {
        if (!segment->as_soma()) {
            auto length = segment->as_cable()->length();
            auto n = (unsigned) std::ceil(length / params.seg_res);
            segment->set_compartments(n);
            tot_seg += n;
        }
    }

    std::cout << "total segments: " << tot_seg << std::endl;

    int tot = 0;
    for (auto layer: syn_layers.map) {
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
            tot++;
        }
    }
    std::cout << "total synapses: " << tot << std::endl;

    for (auto& segment: cell.segments()) {
        arb::mechanism_desc hh("hh");
        hh["gnabar"] = params.hh_gnabar;
        hh["gkbar"] = params.hh_gkbar;
        hh["gl"] = params.hh_gl;
        hh["ena"] = params.hh_ena;
        hh["ek"] = params.hh_ek;

        segment->add_mechanism(hh);
        segment->rL = params.ra;
        segment->cm = params.cm/100; //convert to right unit
    }

    return cell;
}

