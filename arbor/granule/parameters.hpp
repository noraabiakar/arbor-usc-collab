#include <iostream>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

#include <arbor/cable_cell.hpp>
#include <common/json_params.hpp>

std::vector<double> read_spike_times();

struct granule_params {
    std::string morph_file;
    std::vector<float> spikes;

    double temp, v_init;

    double tau1_reg, tau2_reg, e_reg;
    double tau1_syn, tau2_syn, e_syn;

    double gnatbar_ichan2, gkfbar_ichan2, gksbar_ichan2, gkabar_borgka, gncabar_nca, glcabar_lca;
    double gcatbar_cat, gskbar_gskch, gkbar_cagk, gl_ichan2, catau_ccanl, caiinf_ccanl;

    double enat, ekf, eks, ek, elca, etca, esk, el_ichan2, cao;

    double ra_mult, cm_mult;
    double ra, cm;

    std::string syn_layer;
    unsigned syn_id;
    double seg_res, weight;
    double run_time, dt;


};

granule_params read_params(int argc, char** argv) {
    granule_params p;

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

    param_from_json(p.spikes, "spikes", json);

    param_from_json(p.temp, "temp", json);
    param_from_json(p.v_init, "vinit", json);
    param_from_json(p.dt, "dt_arbor", json);
    param_from_json(p.run_time, "run_time", json);
    param_from_json(p.tau1_reg, "tau1_reg", json);
    param_from_json(p.tau2_reg, "tau2_reg", json);
    param_from_json(p.e_reg, "e_reg", json);
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);

    param_from_json(p.gnatbar_ichan2, "gnatbar_ichan2", json);
    param_from_json(p.gkfbar_ichan2, "gkfbar_ichan2", json);
    param_from_json(p.gksbar_ichan2, "gksbar_ichan2", json);
    param_from_json(p.gkabar_borgka, "gkabar_borgka", json);
    param_from_json(p.gncabar_nca, "gncabar_nca", json);
    param_from_json(p.glcabar_lca, "glcabar_lca", json);
    param_from_json(p.gcatbar_cat, "gcatbar_cat", json);
    param_from_json(p.gskbar_gskch, "gskbar_gskch", json);
    param_from_json(p.gkbar_cagk, "gkbar_cagk", json);
    param_from_json(p.gl_ichan2, "gl_ichan2", json);
    param_from_json(p.catau_ccanl, "catau_ccanl", json);
    param_from_json(p.caiinf_ccanl, "caiinf_ccanl", json);

    param_from_json(p.ra, "ra", json);
    param_from_json(p.cm, "cm", json);
    param_from_json(p.ra_mult, "ra_mult", json);
    param_from_json(p.cm_mult, "cm_mult", json);

    param_from_json(p.enat, "enat", json);
    param_from_json(p.ekf, "ekf", json);
    param_from_json(p.eks, "eks", json);
    param_from_json(p.ek, "ek", json);
    param_from_json(p.elca, "elca", json);
    param_from_json(p.etca, "etca", json);
    param_from_json(p.esk, "esk", json);
    param_from_json(p.el_ichan2, "el_ichan2", json);
    param_from_json(p.cao, "cao", json);

    param_from_json(p.syn_layer, "syn_layer", json);
    param_from_json(p.syn_id, "syn_id", json);
    param_from_json(p.seg_res, "seg_res", json);
    param_from_json(p.weight, "weight", json);
    param_from_json(p.morph_file, "morph_file", json);


    for (auto it=json.begin(); it!=json.end(); ++it) {
        std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
    }
    std::cout << "\n";

    return p;

}

/*---------------------------------------------------------------------------------------*/

template <typename T>
struct layer_map {
    std::map<std::string, T> map =
            { {"granuleCellLayer", T()},
              {"innerThird",   T()},
              {"middleThird",  T()},
              {"outerThird",   T()},
              {"soma",    T()}
            };
};

struct synapse_id {
    unsigned segment;
    double pos;
};

using layer_to_bool         = layer_map<bool>;
using layer_to_vec_double   = layer_map<std::vector<double>>;
using layer_to_vec_unsigned = layer_map<std::vector<unsigned>>;
using layer_to_pair_double  = layer_map<std::pair<double, double>>;
using layer_to_vec_syn      = layer_map<std::vector<synapse_id>>;

struct cell_layers {
    layer_to_vec_syn      synapses; // map: layer -> vector of (branch_id, pos)
    layer_to_vec_unsigned segments; // map: layer -> vector of branch_ids belonging in the layer
};

// Returns a vector containing a, b and num_pts-2 equidistant points between a and b.
// Returns empty vector if num_pts == 0
// Returns vector containing a if num_pts == 1
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

// x, y, z: vectors containing coordinates
// Calculate the distance between 2 consecutive points in x, y, z
//     dist = sqrt[(x[i-1]-x[i])^2 + (y[i-1]-y[i])^2) + (z[i-1]-z0[i])^2)]
// Normalize by length
//     norm = dist/length
// Return vector of norm: doubles in [0, 1] that represent locations along a branch
//
// |---|---|----|---|-----------|----------|----------------|
// 0   0.1 0.2  0.3 0.4         0.6        0.75             1
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

// x, y, z: vectors containing coordinates
// Return vector of distance between every point in (x,y,z) and center
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

// Returns layer information:
// cell_layers::synapses is a map from every layer to a vector of synapses (defined as [branch_id, pos] pairs)
// cell_layers::segments is a map from every layer to a vector of branch_ids that lie fully in the layer
cell_layers get_layer_info(std::string filename, double res) {
    cell_layers layer_info;

    std::ifstream f(filename);
    if (!f) throw std::runtime_error("unable to open file");

    // Construct a temporary cell, will NOT be used in the simulation
    auto morph = arb::swc_as_morphology(arb::parse_swc_file(f));
    arb::cable_cell cell = arb::make_cable_cell(morph);

    // Find the furthest swc point from the soma (projected on the z axis)
    double max_extent = 0.0;

    auto soma_loc = cell.soma()->center();
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

    // branch_layer_dist will contain information for every branch of the cell
    // Every branch has a map: layer -> locations on the branch that lie in this layer
    std::vector<layer_to_vec_double> branch_layer_dist;

    // ncomp will contain the number of compartments in every branch of the cell
    std::vector<unsigned> ncomp;

    // Set the layer extents
    layer_to_pair_double layer_extents;
    layer_extents.map["soma"] = {0,0};
    layer_extents.map["granuleCellLayer"] = {0,0.1*max_extent};
    layer_extents.map["innerThird"] = {0.1*max_extent,0.3*max_extent};
    layer_extents.map["middleThird"] =  {0.3*max_extent,0.6*max_extent};
    layer_extents.map["outerThird"] = {0.6*max_extent,max_extent};


    unsigned seg_id = 0;
    for (auto& seg: cell.segments()) {
        if (!seg->is_soma()) {
            std::vector<double> x_interp, y_interp, z_interp;

            // Find and save num compartment in this branch
            auto length = seg->as_cable()->length();
            auto num_comp = (unsigned)std::ceil(length/res);

            ncomp.push_back(num_comp);

            // Set compartments of branch remove
            seg->as_cable()->set_compartments(num_comp);

            // Get all the points that define the branch (including first and last point)
            auto locs = seg->as_cable()->locations();

            // Generate vectors of interpolated points along the branch
            // that are equidistant between 2 consecutive locations in `locs`
            for (unsigned i = 0; i < locs.size()- 1; i++) {
                auto x0 = locs[i].x;
                auto y0 = locs[i].y;
                auto z0 = locs[i].z;
                auto x1 = locs[i + 1].x;
                auto y1 = locs[i + 1].y;
                auto z1 = locs[i + 1].z;
                auto x_t = linspace(x0, x1, num_comp+2);
                auto y_t = linspace(y0, y1, num_comp+2);
                auto z_t = linspace(z0, z1, num_comp+2);
                x_interp.insert(x_interp.end(), x_t.begin(), x_t.end()-1);
                y_interp.insert(y_interp.end(), y_t.begin(), y_t.end()-1);
                z_interp.insert(z_interp.end(), z_t.begin(), z_t.end()-1);
            }
            x_interp.push_back(locs.back().x);
            y_interp.push_back(locs.back().y);
            z_interp.push_back(locs.back().z);

            // Normalize the interpolated points to the length of the cable
            auto norm = normalize(x_interp, y_interp, z_interp, seg->as_cable()->length());
            // Generate the distance from the inetrpolated points to the soma
            auto ext =     extent(x_interp, y_interp, z_interp, cell.soma()->center());


            layer_to_bool       branch_in_layer;           // layer -> bool that indicates whether or not the branch is attributed to the layer
            layer_to_vec_double branch_layer_to_positions; // layer -> subsection of `norm` that is within the boundaries of the layer

            for (unsigned i = 0; i< ext.size(); i++) {
                for (auto e: layer_extents.map) {
                    auto layer_name = e.first;
                    auto layer_extents = e.second;

                    if (ext[i] >= layer_extents.first && ext[i] < layer_extents.second) {
                        branch_layer_to_positions.map[layer_name].push_back(norm[i]);

                        // if any of the points of the branch lie in a layer,
                        // then the whole branch is attributed to that layer
                        branch_in_layer.map[e.first] = true;
                    }
                }
            }

            // Add the branch info to the vector of all branch infos
            branch_layer_dist.push_back(std::move(branch_layer_to_positions));

            // Attribute the branch id to the correct layer(s)
            for (auto l: branch_in_layer.map) {
                if(l.second) {
                    layer_info.segments.map[l.first].push_back(seg_id);
                }
            }

        }
        seg_id++;
    }

    // At this point branch_layer_dist contains all the possible synapse locations per layer for every branch
    // However, we only want one synapse per compartment per branch.
    // So, we reduce branch_layer_dist and save the final results in layer_info::synapses

    // For every branch
    for (unsigned i = 0; i < branch_layer_dist.size(); i++) {

        // n is the number of compartments in a branch, so it is also the number of synapses we want
        unsigned n = ncomp[i];

        // Generate n reference points on the branch
        auto branch_pos = linspace(0.5/n, 1-0.5/n, n);

        std::vector<double> abs_diff;

        for (auto layer: branch_layer_dist[i].map) {
            auto layer_name = layer.first;
            auto branch_possible_syanpses = layer.second;

            auto branch_selected_synapses = layer_info.synapses.map[layer_name];
            auto last_synapse_index = branch_selected_synapses.end() - branch_selected_synapses.begin();

            for (auto s : branch_possible_syanpses) {
                // for every possible synapse, find the distance from the reference points, and store in abs_diff
                abs_diff = branch_pos;
                for (auto& a: abs_diff) {
                    a = std::abs(a - s);
                }

                // Get the closest point from the reference points to the synapse
                unsigned idx = std::min_element(abs_diff.begin(), abs_diff.end()) - abs_diff.begin();

                // Use that point as a synapse on the layer
                layer_info.synapses.map[layer_name].push_back({i+1, branch_pos[idx]});
            }

            // Remove duplicate synapses
            auto last = std::unique(last_synapse_index + layer_info.synapses.map[layer_name].begin(),
                                    layer_info.synapses.map[layer_name].end(),
                                    [](synapse_id a, synapse_id b) { return a.segment == b.segment && a.pos == b.pos;}
                                    );
            layer_info.synapses.map[layer_name].erase(last, layer_info.synapses.map[layer_name].end());
        }
    }

    layer_info.synapses.map["soma"].push_back({0, 0.5});


    // Sort the vector of branch_ids for every layer (for binary search later)
    for (auto& s: layer_info.segments.map) {
        std::sort(s.second.begin(), s.second.end());
    }

    return layer_info;
}

