#include <arbor/mechcat.hpp>
#include <arbor/version.hpp>

#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif
#include "backends/multicore/fvm.hpp"

#include "usc_catalogue.hpp"

#include "mechanisms/borgka.hpp"
#include "mechanisms/cagk.hpp"
#include "mechanisms/cat.hpp"
#include "mechanisms/ccanl.hpp"
#include "mechanisms/ccanlrev.hpp"
#include "mechanisms/gskch.hpp"
#include "mechanisms/ichan2.hpp"
#include "mechanisms/lca.hpp"
#include "mechanisms/nca.hpp"

#include "../gtest.h"

#ifndef ARB_GPU_ENABLED
#define ADD_MECH(c, x)\
c.add(#x, usc::mechanism_##x##_info());\
c.register_implementation(#x, usc::make_mechanism_##x<multicore::backend>());
#else
#define ADD_MECH(c, x)\
c.add(#x, usc::mechanism_##x##_info());\
c.register_implementation(#x, usc::make_mechanism_##x<multicore::backend>());\
c.register_implementation(#x, usc::make_mechanism_##x<gpu::backend>());
#endif

using namespace arb;

mechanism_catalogue make_usc_catalogue() {
    mechanism_catalogue usc_cat(global_default_catalogue());

    ADD_MECH(usc_cat, borgka)
    ADD_MECH(usc_cat, cagk)
    ADD_MECH(usc_cat, cat)
    ADD_MECH(usc_cat, ccanl)
    ADD_MECH(usc_cat, ccanlrev)
    ADD_MECH(usc_cat, gskch)
    ADD_MECH(usc_cat, ichan2)
    ADD_MECH(usc_cat, lca)
    ADD_MECH(usc_cat, nca)

    return usc_cat;
}

