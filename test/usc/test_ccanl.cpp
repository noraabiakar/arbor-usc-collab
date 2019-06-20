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
    // and ccanl, ccanlrev mechanisms.

    cable_cell cell;

    auto soma = cell.add_soma(6.);
    soma->add_mechanism("ccanl");

    cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
    cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";

    std::unique_ptr<cable1d_recipe> rec(new cable1d_recipe(cell));
    rec->catalogue() = make_usc_catalogue();

    // (Interior concentration and reversal potentials are provided by mechanisms.)
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

    // Expect initial concentrations to be set to values from ccanl,
    // reversal potentials from ccanlrev.

    auto& state = *(fvcell.*private_state_ptr).get();

    auto nca_ion = state.ion_data.at("nca");
    auto lca_ion = state.ion_data.at("lca");
    auto tca_ion = state.ion_data.at("tca");

    double expected_xcai = 50.e-6/3;
    EXPECT_DOUBLE_EQ(expected_xcai, nca_ion.Xi_[0]);
    EXPECT_DOUBLE_EQ(expected_xcai, lca_ion.Xi_[0]);
    EXPECT_DOUBLE_EQ(expected_xcai, tca_ion.Xi_[0]);
}
