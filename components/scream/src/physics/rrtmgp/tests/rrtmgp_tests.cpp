#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/mo_garand_atmos_io.h"
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"

#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

#include "cpp/rrtmgp/mo_gas_concentrations.h"

#include "YAKL.h"
#include "ekat/util/ekat_test_utils.hpp"

#include <iostream>
#include <cmath>

using namespace scream;

void expect_another_arg (int i, int argc) {
  EKAT_REQUIRE_MSG(i != argc-1, "Expected another cmd-line arg.");
}

int main(int argc, char** argv) {
    // Parse command line arguments
    if (argc < 3) {
        std::cout << "argc: " << argc << std::endl;
        std::cout << argv[0] << " -i <inputfile> -b <baseline file> [options]\n";
        return 1;
    }
    std::string inputfile, baseline, device;

    for (int i = 1; i < argc-1; ++i) {
      if (ekat::argv_matches(argv[i], "-b", "--baseline-file")) {
        expect_another_arg(i, argc);
        ++i;
        baseline = argv[i];
      }
      if (ekat::argv_matches(argv[i], "-i", "--input-file")) {
        expect_another_arg(i, argc);
        ++i;
        inputfile = argv[i];
      }
      // RRTMGP baselines tests to not use kokoks. Swallow the arg, but ignore it
      if (std::string(argv[i])=="--ekat-kokkos-device") {
        expect_another_arg(i, argc);
        ++i;
        device = argv[i];
      }
    }

    // Check to see that inputfiles exist
    if (!rrtmgpTest::file_exists(inputfile.c_str())) {
        std::cout << "Inputfile " << inputfile << " does not exist." << std::endl;
        return -1;
    }
    if (!rrtmgpTest::file_exists(baseline.c_str())) {
        std::cout << "Baseline " << baseline << " does not exist." << std::endl;
        return -1;
    }

    // Initialize yakl
    std::cout << "Initialize yakl...\n";
    yakl::init();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    std::cout << "Read fluxes...\n";
    real2d sw_flux_up_ref;
    real2d sw_flux_dn_ref;
    real2d sw_flux_dir_ref;
    real2d lw_flux_up_ref;
    real2d lw_flux_dn_ref;
    rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

    // Get dimension sizes
    int ncol = sw_flux_up_ref.dimension[0];
    int nlev = sw_flux_up_ref.dimension[1];
    int nlay = nlev - 1;

    // Read in dummy Garand atmosphere; if this were an actual model simulation, 
    // these would be passed as inputs to the driver
    // NOTE: set ncol to size of col_flx dimension in the input file. This is so
    // that we can compare to the reference data provided in that file. Note that
    // this will copy the first column of the input data (the first profile) ncol
    // times. We will then fill some fraction of these columns with clouds for
    // the test problem.
    std::cout << "Read dummy atmos...\n";
    real2d p_lay("p_lay", ncol, nlay);
    real2d t_lay("t_lay", ncol, nlay);
    real2d p_lev("p_lev", ncol, nlay+1);
    real2d t_lev("t_lev", ncol, nlay+1);
    GasConcs gas_concs;
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, ncol);

    // Initialize absorption coefficients
    std::cout << "Initialize RRTMGP...\n";
    scream::rrtmgp::rrtmgp_initialize(gas_concs);

    // Setup our dummy atmosphere based on the input data we read in
    std::cout << "Setup dummy atmos...\n";
    real1d sfc_alb_dir_vis("sfc_alb_dir_vis", ncol);
    real1d sfc_alb_dir_nir("sfc_alb_dir_nir", ncol);
    real1d sfc_alb_dif_vis("sfc_alb_dif_vis", ncol);
    real1d sfc_alb_dif_nir("sfc_alb_dif_nir", ncol);
    real1d mu0("mu0", ncol);
    real2d lwp("lwp", ncol, nlay);
    real2d iwp("iwp", ncol, nlay);
    real2d rel("rel", ncol, nlay);
    real2d rei("rei", ncol, nlay);
    real2d cld("cld", ncol, nlay);
    rrtmgpTest::dummy_atmos(
            inputfile, ncol, p_lay, t_lay,
            sfc_alb_dir_vis, sfc_alb_dir_nir,
            sfc_alb_dif_vis, sfc_alb_dif_nir,
            mu0,
            lwp, iwp, rel, rei, cld
        );

    // Setup flux outputs; In a real model run, the fluxes would be
    // input/outputs into the driver (persisting between calls), and
    // we would just have to setup the pointers to them in the
    // FluxesBroadband object
    std::cout << "Setup fluxes...\n";
    const auto nswbands = scream::rrtmgp::k_dist_sw.get_nband();
    const auto nlwbands = scream::rrtmgp::k_dist_lw.get_nband();
    real2d sw_flux_up ("sw_flux_up" , ncol, nlay+1);
    real2d sw_flux_dn ("sw_flux_dn" , ncol, nlay+1);
    real2d sw_flux_dir("sw_flux_dir", ncol, nlay+1);
    real2d lw_flux_up ("lw_flux_up" , ncol, nlay+1);
    real2d lw_flux_dn ("lw_flux_dn" , ncol, nlay+1);
    real3d sw_bnd_flux_up ("sw_bnd_flux_up" , ncol, nlay+1, nswbands);
    real3d sw_bnd_flux_dn ("sw_bnd_flux_dn" , ncol, nlay+1, nswbands);
    real3d sw_bnd_flux_dir("sw_bnd_flux_dir", ncol, nlay+1, nswbands);
    real3d lw_bnd_flux_up ("lw_bnd_flux_up" , ncol, nlay+1, nlwbands);
    real3d lw_bnd_flux_dn ("lw_bnd_flux_dn" , ncol, nlay+1, nlwbands);

    // Compute band-by-band surface_albedos.
    real2d sfc_alb_dir("sfc_alb_dir", ncol, nswbands);
    real2d sfc_alb_dif("sfc_alb_dif", ncol, nswbands);
    rrtmgp::compute_band_by_band_surface_albedos(
      ncol, nswbands,
      sfc_alb_dir_vis, sfc_alb_dir_nir,
      sfc_alb_dif_vis, sfc_alb_dif_nir,
      sfc_alb_dir, sfc_alb_dif);

    // Run RRTMGP code on dummy atmosphere
    std::cout << "Run RRTMGP...\n";
    scream::rrtmgp::rrtmgp_main(
            ncol, nlay,
            p_lay, t_lay, p_lev, t_lev, gas_concs,
            sfc_alb_dir, sfc_alb_dif, mu0,
            lwp, iwp, rel, rei,
            sw_flux_up, sw_flux_dn, sw_flux_dir,
            lw_flux_up, lw_flux_dn,
            sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
            lw_bnd_flux_up, lw_bnd_flux_dn);

    // Check values against baseline
    std::cout << "Check values...\n";
    rrtmgpTest::read_fluxes(
        baseline, 
        sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dir_ref,
        lw_flux_up_ref, lw_flux_dn_ref
    );
    int nerr = 0;
    if (!rrtmgpTest::all_close(sw_flux_up_ref , sw_flux_up , 0.001)) nerr++;
    if (!rrtmgpTest::all_close(sw_flux_dn_ref , sw_flux_dn , 0.001)) nerr++;
    if (!rrtmgpTest::all_close(sw_flux_dir_ref, sw_flux_dir, 0.001)) nerr++;
    if (!rrtmgpTest::all_close(lw_flux_up_ref , lw_flux_up , 0.001)) nerr++;
    if (!rrtmgpTest::all_close(lw_flux_dn_ref , lw_flux_dn , 0.001)) nerr++;

    // Clean up or else YAKL will throw errors
    scream::rrtmgp::rrtmgp_finalize();
    sw_flux_up_ref.deallocate();
    sw_flux_dn_ref.deallocate();
    sw_flux_dir_ref.deallocate();
    lw_flux_up_ref.deallocate();
    lw_flux_dn_ref.deallocate();
    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();
    sw_bnd_flux_up.deallocate();
    sw_bnd_flux_dn.deallocate();
    sw_bnd_flux_dir.deallocate();
    lw_bnd_flux_up.deallocate();
    lw_bnd_flux_dn.deallocate();
    p_lay.deallocate();
    t_lay.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();
    gas_concs.reset();
    sfc_alb_dir_vis.deallocate();
    sfc_alb_dir_nir.deallocate();
    sfc_alb_dif_vis.deallocate();
    sfc_alb_dif_nir.deallocate();
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();
    mu0.deallocate();
    lwp.deallocate();
    iwp.deallocate();
    rel.deallocate();
    rei.deallocate();
    cld.deallocate();
    yakl::finalize();

    return nerr != 0 ? 1 : 0;

}  // end of main driver code
