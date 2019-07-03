#include <cmath>
#include <memory>

#include <arbor/recipe.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>

// Inspecting FVM shared state internals.
#include "backends/multicore/fvm.hpp"
#include "fvm_lowered_cell_impl.hpp"
#include "execution_context.hpp"

#include "../simple_recipes.hpp"
#include "../gtest.h"
#include "usc_catalogue.hpp"

using recipe_ptr = std::unique_ptr<arb::recipe>;
using fvm_cell = arb::fvm_lowered_cell_impl<arb::multicore::backend>;

using namespace arb;

static recipe_ptr make_test_recipe(double inca = 0, double ilca = 0, double itca = 0) {
    // Make a single compartment cell with nca, lca, tca ions
    // and ccanl, ccanlrev, lca mechanisms.
    //
    // We need at least one mechanism that reads a reversal potential
    // written by ccanlrev for us to be able to check the generated
    // value.

    mechanism_catalogue cat = make_usc_catalogue();
    cat.derive("constant/nca", "constant_ionic_current", {{"current", inca}}, {{"x", "nca"}});
    cat.derive("constant/lca", "constant_ionic_current", {{"current", ilca}}, {{"x", "lca"}});
    cat.derive("constant/tca", "constant_ionic_current", {{"current", itca}}, {{"x", "tca"}});

    cable_cell cell;

    auto soma = cell.add_soma(6.);
    soma->add_mechanism("ccanl");
    soma->add_mechanism("constant/nca");
    soma->add_mechanism("constant/lca");
    soma->add_mechanism("constant/tca");

    cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";

    std::unique_ptr<cable1d_recipe> rec(new cable1d_recipe(cell));
    rec->catalogue() = cat;

    rec->add_ion("nca", 2, 0, 0, 0);
    rec->add_ion("lca", 2, 0, 0, 0);
    rec->add_ion("tca", 2, 0, 0, 0);

    return recipe_ptr(rec.release());
}

// Private member access tricks, copied from unit tests:

namespace access_ {
    template <typename V, V& store, V value>
    struct bind {
        static struct binder {
            binder() { store = value; }
        } init;
    };

    template <typename V, V& store, V value>
    typename bind<V, store, value>::binder bind<V, store, value>::init;
} // namespace access

#define ACCESS_BIND(type, global, value)\
namespace { using global ## _type_ = type; global ## _type_ global; }\
template struct ::access_::bind<type, global, value>;

using shared_state = arb::multicore::shared_state;
ACCESS_BIND(std::unique_ptr<shared_state> fvm_cell::*, private_state_ptr, &fvm_cell::state_)

TEST(usc, ccanl_init) {
    auto rec = make_test_recipe();

    // Make an fvm lowered cell directly:

    std::vector<target_handle> targets;
    std::vector<fvm_index_type> cell_to_intdom;
    probe_association_map<probe_handle> probe_map;

    execution_context context;
    fvm_cell fvcell{context};
    fvcell.initialize({0}, *rec.get(), cell_to_intdom, targets, probe_map);

    auto& state = *(fvcell.*private_state_ptr).get();

    auto nca_ion = state.ion_data.at("nca");
    auto lca_ion = state.ion_data.at("lca");
    auto tca_ion = state.ion_data.at("tca");

    double expected_cai = 50.e-6/3;
    EXPECT_DOUBLE_EQ(expected_cai, nca_ion.Xi_[0]);
    EXPECT_DOUBLE_EQ(expected_cai, lca_ion.Xi_[0]);
    EXPECT_DOUBLE_EQ(expected_cai, tca_ion.Xi_[0]);

    double expected_erev = 1000*8.3134*(6.3+273.15)/(2*96520)*std::log(2.0/50.e-6);
    EXPECT_DOUBLE_EQ(expected_erev, nca_ion.eX_[0]);
    EXPECT_DOUBLE_EQ(expected_erev, lca_ion.eX_[0]);
    EXPECT_DOUBLE_EQ(expected_erev, tca_ion.eX_[0]);
}

#define EXPECT_NEAR_RELATIVE(expected, value, reltol)\
EXPECT_NEAR((expected), (value), std::abs(expected)*(reltol))

TEST(usc, ccanl_dynamic) {
    const double inca = 1.2;
    const double ilca = 3.4;
    const double itca = 5.6;

    const double depth = 200;
    const double catau = 9;
    const double xi0 = 50.e-6/3;

    auto rec = make_test_recipe(inca, ilca, itca);

    // Each of the ion concentrations is governed by
    //    c' = -a·i + (c₀-c)/τ
    // which with constant current i gives
    //    c(t) = c₀ - a·i·τ·(1 - exp(-t/t)).
    //
    // a is equal to 1e6/(depth*96520).

    auto analytic = [=](double t, double ix) {
        return xi0 - 1e6/(depth*96520)*ix*catau*(1-std::exp(-t/catau));
    };

    // Make an fvm lowered cell directly:

    std::vector<target_handle> targets;
    std::vector<fvm_index_type> cell_to_intdom;
    probe_association_map<probe_handle> probe_map;

    execution_context context;
    fvm_cell fvcell{context};
    fvcell.initialize({0}, *rec.get(), cell_to_intdom, targets, probe_map);

    auto& state = *(fvcell.*private_state_ptr).get();

    auto& nca_ion = state.ion_data.at("nca");
    auto& lca_ion = state.ion_data.at("lca");
    auto& tca_ion = state.ion_data.at("tca");

    (void)fvcell.integrate(1.0, 0.001, {}, {});

    double expected_ncai = analytic(1., inca);
    double expected_lcai = analytic(1., ilca);
    double expected_tcai = analytic(1., itca);

    EXPECT_NEAR_RELATIVE(expected_ncai, nca_ion.Xi_[0], 1e-6);
    EXPECT_NEAR_RELATIVE(expected_lcai, lca_ion.Xi_[0], 1e-6);
    EXPECT_NEAR_RELATIVE(expected_tcai, tca_ion.Xi_[0], 1e-6);
}
