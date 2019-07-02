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
    std::vector<double> spikes = {
            25.269724183039855, 29.37076391451496, 58.472477010286546, 93.80268485203328,
            112.71090127018375, 142.6472406502223/*, 293.3318516075217, 456.96763867081177,
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
            9956.054906300105*/
    };
    double temp, v_init;
    double tau1_reg, tau2_reg, e_reg;
    double tau1_syn, tau2_syn, e_syn;

    double hh_gnabar, hh_gkbar, hh_gl, hh_ena, hh_ek;

    double gnatbar_ichan2, gkfbar_ichan2, gksbar_ichan2, gkabar_borgka, gncabar_nca, glcabar_lca;
    double gcatbar_cat, gskbar_gskch, gkbar_cagk, gl_ichan2, catau_ccanl, caiinf_ccanl;

    double enat, ekf, eks, ek, elca, etca, esk, el_ichan2, cao;

    double ra_mult, cm_mult;
    double ra, cm;

    std::string syn_layer;
    unsigned syn_id;
    double seg_res, dt, weight;


};

granule_params read_params(int argc, char** argv) {
    granule_params p;
    p.morph_file = "/home/abiakarn/git/usc_neuron/morphologies/output2_updated.swc";

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
    param_from_json(p.tau1_reg, "tau1_reg", json);
    param_from_json(p.tau2_reg, "tau2_reg", json);
    param_from_json(p.e_reg, "e_reg", json);
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);

    param_from_json(p.hh_gnabar, "hh_gnabar", json);
    param_from_json(p.hh_gkbar, "hh_gkbar", json);
    param_from_json(p.hh_gl, "hh_gl", json);
    param_from_json(p.hh_ena, "hh_ena", json);
    param_from_json(p.hh_ek, "hh_ek", json);

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
struct layers {
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

using double_layers = layers<std::vector<double>>;
using unsigned_layers = layers<std::vector<unsigned>>;
using pair_layers = layers<std::pair<double, double>>;
using synapse_layers = layers<std::vector<synapse_id>>;

struct cell_layers {
    synapse_layers synapses;
    unsigned_layers segments;
};

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

cell_layers get_layer_info(std::string filename, double res) {
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

    std::vector<double_layers> segment_layer_dist;
    std::vector<unsigned> nseg;

    pair_layers layer_extents;
    layer_extents.map["soma"] = {0,0};
    layer_extents.map["granuleCellLayer"] = {0,0.1*max_extent};
    layer_extents.map["innerThird"] = {0.1*max_extent,0.3*max_extent};
    layer_extents.map["middleThird"] =  {0.3*max_extent,0.6*max_extent};
    layer_extents.map["outerThird"] = {0.6*max_extent,max_extent};

    cell_layers ret;

    unsigned seg_id = 0;
    for (auto& seg: cell.segments()) {
        if (!seg->is_soma()) {
            std::vector<double> x_interp, y_interp, z_interp;
            auto length = seg->as_cable()->length();
            auto n = (unsigned)std::ceil(length/res);
            nseg.push_back(n);
            seg->as_cable()->set_compartments(n);

            auto locs = seg->as_cable()->locations();
            for (unsigned i = 0; i < locs.size()- 1; i++) {
                auto x0 = locs[i].x;
                auto y0 = locs[i].y;
                auto z0 = locs[i].z;
                auto x1 = locs[i + 1].x;
                auto y1 = locs[i + 1].y;
                auto z1 = locs[i + 1].z;
                auto x_t = linspace(x0, x1, n+2);
                auto y_t = linspace(y0, y1, n+2);
                auto z_t = linspace(z0, z1, n+2);
                x_interp.insert(x_interp.end(), x_t.begin(), x_t.end()-1);
                y_interp.insert(y_interp.end(), y_t.begin(), y_t.end()-1);
                z_interp.insert(z_interp.end(), z_t.begin(), z_t.end()-1);
            }
            x_interp.push_back(locs.back().x);
            y_interp.push_back(locs.back().y);
            z_interp.push_back(locs.back().z);

            auto norm = normalize(x_interp, y_interp, z_interp, seg->as_cable()->length());
            auto ext = extent(x_interp, y_interp, z_interp, cell.soma()->center());

            double_layers l;
            for (unsigned i = 0; i< ext.size(); i++) {
                for (auto e: layer_extents.map) {
                    if (ext[i] >= e.second.first && ext[i] < e.second.second) {
                        l.map[e.first].push_back(norm[i]);
                        ret.segments.map[e.first].push_back(seg_id);
                    }
                }
            }
            segment_layer_dist.push_back(std::move(l));
        }
        seg_id++;
    }
    for (auto s: ret.segments.map) {
        std::sort(s.second.begin(), s.second.end());
        auto last = std::unique(s.second.begin(), s.second.end());
        s.second.erase(last, s.second.end());
    }

    for (unsigned i = 0; i < segment_layer_dist.size(); i++) {
        unsigned n = nseg[i];
        auto branch_pos = linspace(0.5/n, 1-0.5/n, n);
        std::vector<double> abs_diff;

        for (auto layer: segment_layer_dist[i].map) {
            auto prev_end_idx = ret.synapses.map[layer.first].end() - ret.synapses.map[layer.first].begin();

            for (auto v : layer.second) {
                abs_diff = branch_pos;
                for (auto& a: abs_diff) {
                    a = std::abs(a - v);
                }
                unsigned idx = std::min_element(abs_diff.begin(), abs_diff.end() ) - abs_diff.begin();
                ret.synapses.map[layer.first].push_back({i+1, branch_pos[idx]});
            }
            auto last = std::unique(prev_end_idx + ret.synapses.map[layer.first].begin(), ret.synapses.map[layer.first].end(),
                    [](synapse_id a, synapse_id b) { return a.segment == b.segment && a.pos == b.pos; });
            ret.synapses.map[layer.first].erase(last, ret.synapses.map[layer.first].end());
        }
    }

    ret.synapses.map["soma"].push_back({0, 0.5});

    return ret;
}

