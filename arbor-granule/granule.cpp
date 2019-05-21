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
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Writes voltage trace as a json file.
 void write_trace_json(const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell granule_cell(arb::cell_gid_type gid, std::string filename);

class granule_recipe: public arb::recipe {
public:
    granule_recipe(const std::vector<double>& spikes, granule_params gparams):
        num_cells_(1), spikes_(spikes), params_(gparams) {}

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return granule_cell(gid, params_.morph_file);
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
        std::vector<arb::event_generator> gens;
        arb::pse_vector svec;
        for (auto s: spikes_) {
            svec.push_back({{0, 0}, s, event_weight_});
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
        a.temperature_K = 308.15;
        a.init_membrane_potential_mV = -70;
        return a;
    }

private:
    cell_size_type num_cells_;
    std::vector<double> spikes_;
    granule_params params_;
    float event_weight_ = 1.438995e-04;
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

        // Create an instance of our recipe.
        granule_recipe recipe(read_spike_times(), read_params());

        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};
        // The schedule for sampling is 10 samples every 1 ms.
        auto sched = arb::regular_schedule(0.1);
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
        sim.run(10000, 0.025);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated at rate of "
                      << 10000 << " ms between spikes\n";
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

std::vector<double> linspace(double a, double b, unsigned num_pts) {

    std::vector<double> ret;
    if (num_pts == 0) {
        return ret;
    }

    ret.push_back(a);
    if (num_pts == 1) {
        return ret;
    }

    double inc = (b - a)/ (num_pts-1);
    for (unsigned i = 1; i < num_pts; i++) {
        ret.push_back(ret[i-1] + inc);
    }
    return ret;
}

std::vector<double> normalize(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& z,
                              const double length) {

    std::vector<double> norm = {0};
    double dist = 0;
    for (unsigned i = 1; i < x.size(); i++) {
        auto x_sq = std::pow((x[i]-x[i-1]),2);
        auto y_sq = std::pow((y[i]-y[i-1]),2);
        auto z_sq = std::pow((z[i]-z[i-1]),2);
        dist += std::sqrt(x_sq + y_sq + z_sq);
        norm.push_back(dist/length);
    }
    return norm;
}

std::vector<double> extent(const std::vector<double>& x,
                           const std::vector<double>& y,
                           const std::vector<double>& z,
                           arb::point<double> center) {

    std::vector<double> ext;
    for (unsigned i = 0; i < x.size(); i++) {
        arb::point<double> loc = {x[i], y[i], z[i]};
        auto diff = std::abs((loc-center).z);
        ext.push_back(diff);
    }
    return ext;
}

struct layers {
    std::unordered_map<std::string, std::vector<double>> map =
            { {"soma_layer",    {}},
              {"granule_layer", {}},
              {"inner_layer",   {}},
              {"middle_layer",  {}},
              {"outer_layer",  {}}
            };
};

arb::cable_cell granule_cell(arb::cell_gid_type gid, std::string filename) {
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("unable to open file");

    auto morph = arb::swc_as_morphology(arb::parse_swc_file(f));
    arb::cable_cell cell = arb::make_cable_cell(morph);

    auto soma_loc = cell.soma()->center();
    double max_extent = 0.0;
    for (auto& seg: cell.segments()) {
        if (!seg->is_soma()) {
            auto locs = seg->as_cable()->locations();
            for (auto l :locs) {
                auto diff = l - soma_loc;
                auto dist = diff.z;
                if (dist > max_extent) {
                    max_extent = dist;
                }
            }
        }
    }

    std::vector<layers> segment_layer_dist;

    std::unordered_map<std::string, std::pair<double, double>> layer_extents =
            { {"soma_layer",    {0,0}},
              {"granule_layer", {0,0.1*max_extent}},
              {"inner_layer",   {0.1*max_extent,0.3*max_extent}},
              {"middle_layer",  {0.3*max_extent,0.6*max_extent}},
              {"outer_layer",  {0.6*max_extent,max_extent}}
            };

    for (auto& seg: cell.segments()) {
        unsigned nseg = 1; //depends on segment length
        if (!seg->is_soma()) {
            std::vector<double> x_interp, y_interp, z_interp;
            auto locs = seg->as_cable()->locations();
            for (unsigned i = 0; i < locs.size()- 1; i++) {
                auto x0 = locs[i].x;
                auto y0 = locs[i].y;
                auto z0 = locs[i].z;
                auto x1 = locs[i + 1].x;
                auto y1 = locs[i + 1].y;
                auto z1 = locs[i + 1].z;
                auto x_t = linspace(x0, x1, nseg+2);
                auto y_t = linspace(y0, y1, nseg+2);
                auto z_t = linspace(z0, z1, nseg+2);
                x_interp.insert(x_interp.end(), x_t.begin(), x_t.end()-1);
                y_interp.insert(y_interp.end(), y_t.begin(), y_t.end()-1);
                z_interp.insert(z_interp.end(), z_t.begin(), z_t.end()-1);
            }
            x_interp.push_back(locs.back().x);
            y_interp.push_back(locs.back().y);
            z_interp.push_back(locs.back().z);

            auto norm = normalize(x_interp, y_interp, z_interp, seg->as_cable()->length());
            auto ext = extent(x_interp, y_interp, z_interp, cell.soma()->center());

            layers l;
            for (unsigned i = 0; i< ext.size(); i++) {
                for (auto e: layer_extents) {
                    if (ext[i] >= e.second.first && ext[i] < e.second.second) {
                        l.map[e.first].push_back(norm[i]);
                    }
                }
            }
            segment_layer_dist.push_back(std::move(l));
        }
    }

    std::vector<layers> segment_layer_pos;
    for (auto s: segment_layer_dist) {
        unsigned nseg = 1; // property of branch
        auto branch_pos = linspace(0.5/nseg, 1-0.5/nseg, nseg);
        std::vector<double> abs_diff;

        layers l;
        for (auto layer: s.map) {
            for (auto v : layer.second) {
                abs_diff = branch_pos;
                for (auto& a: abs_diff) {
                    a = std::abs(a - v);
                }
                unsigned idx = std::min_element(abs_diff.begin(), abs_diff.end() ) - abs_diff.begin();
                l.map[layer.first].push_back(branch_pos[idx]);
            }

            auto last = std::unique(l.map[layer.first].begin(), l.map[layer.first].end());
            l.map[layer.first].erase(last, l.map[layer.first].end());
        }
        segment_layer_pos.push_back(std::move(l));
    }

    layers soma;
    soma.map["soma_layer"].push_back(0.5);
    segment_layer_pos.insert(segment_layer_pos.begin(), soma);

    for (auto i: segment_layer_pos) {
        std::cout << "[";
        for (auto l : i.map["soma_layer"]) {
            std::cout << l << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "[";
        for (auto l : i.map["granule_layer"]) {
            std::cout << l << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "[";
        for (auto l : i.map["inner_layer"]) {
            std::cout << l << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "[";
        for (auto l : i.map["middle_layer"]) {
            std::cout << l << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "[";
        for (auto l : i.map["outer_layer"]) {
            std::cout << l << " ";
        }
        std::cout << "]" << std::endl << std::endl;
    }

    cell.add_detector({0, 0}, 10);

    arb::mechanism_desc exp2syn("exp2syn");
    exp2syn["tau1"] = 0.5;
    exp2syn["tau2"] = 0.6;
    exp2syn["e"] = 0;

    for (unsigned i = 0; i < segment_layer_pos.size(); i++) {
        for (auto layer: segment_layer_pos[i].map) {
            for (auto j: layer.second) {
                std::cout << i << " " << j <<std::endl;
                cell.add_synapse({i, j}, exp2syn);
            }
        }
        std::cout << std::endl;
    }

    return cell;
}

