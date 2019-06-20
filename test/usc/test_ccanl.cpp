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

static recipe_ptr make_test_recipe() {
    // Make a single compartment cell with nca, lca, tca ions
    // and ccanl, ccanlrev, lca mechanisms.
    //
    // We need at least one mechanism that reads a reversal potential
    // written by ccanlrev for us to be able to check the generated
    // value.

    cable_cell cell;

    auto soma = cell.add_soma(6.);
    soma->add_mechanism("ccanl");
    soma->add_mechanism("lca");

    cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";

    std::unique_ptr<cable1d_recipe> rec(new cable1d_recipe(cell));
    rec->catalogue() = make_usc_catalogue();

    rec->add_ion("nca", 2, 0, 2.0/3, 0);
    rec->add_ion("lca", 2, 0, 2.0/3, 0);
    rec->add_ion("tca", 2, 0, 2.0/3, 0);

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

TEST(usc, ccanl) {
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
