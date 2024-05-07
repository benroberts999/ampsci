#include "DiracOperator.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <utility>

TEST_CASE("DiracOperator", "[DiracOperator][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracOperator\n";

  // for dumb reason, generateOperator takes in wf...
  Wavefunction wf({2000, 1.0e-6, 300.0, 50.0, "loglinear"},
                  {"Cs", 133, "Fermi"});

  auto &orbs = wf.valence();
  for (int ik = 0; ik <= Angular::indexFromKappa(2); ik++) {
    int kappa = Angular::kappaFromIndex(ik);
    int n = Angular::l_k(kappa) + 1;
    orbs.push_back(DiracSpinor::exactHlike(n, kappa, wf.grid_sptr(), 1.0));
  }

  //--------------------------------------------------------------------
  SECTION("E1") {
    std::cout << "E1\n";
    const auto e1data = std::vector{std::tuple{"1s+", "2p-", -1.0534784551e+00},
                                    {"1s+", "2p+", 1.4898420883e+00},
                                    {"2p-", "1s+", -1.0534784551e+00},
                                    {"2p-", "3d-", 5.4823578421e+00},
                                    {"2p+", "1s+", -1.4898420883e+00},
                                    {"2p+", "3d-", -2.4518347077e+00},
                                    {"3d-", "2p-", -5.4823578421e+00},
                                    {"3d-", "2p+", -2.4518347077e+00}};

    const IO::InputBlock options{""};
    auto h = DiracOperator::generate("E1", options, wf);

    for (auto &[a, b, me] : e1data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hab = h->reducedME(Fa, Fb);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
    }
    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 1);
    REQUIRE(h->parity() == -1);
    REQUIRE(h->name() == "E1");
    REQUIRE(!h->units().empty());

    // v-form
    auto hv = DiracOperator::generate("E1v", {}, wf);
    REQUIRE(hv->get_d_order() == 0);
    REQUIRE(hv->imaginaryQ() == false);
    REQUIRE(hv->rank() == 1);
    REQUIRE(hv->parity() == -1);
    REQUIRE(hv->name() == "E1v");
    REQUIRE(!hv->units().empty());
    for (auto &[a, b, me] : e1data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      hv->updateFrequency(Fa.en() - Fb.en());
      const auto hvab = hv->reducedME(Fa, Fb);
      REQUIRE(std::abs(hvab - me) < 1.0e-4);
    }

    auto hk1 = DiracOperator::generate("Ek", {"", "k=1;"}, wf);
    for (auto &[a, b, me] : e1data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hkab = hk1->reducedME(Fa, Fb);
      REQUIRE(std::abs(hkab - me) < 1.0e-4);
    }
  }

  //--------------------------------------------------------------------
  SECTION("E2") {
    std::cout << "Ek\n";
    const auto data = std::vector{std::tuple{"1s+", "3d-", 1.5500083192e+00},
                                  {"2p-", "2p+", 2.6832065561e+01},
                                  {"2p+", "2p-", -2.6832065561e+01},
                                  {"2p+", "2p+", 2.6832553768e+01},
                                  {"3d-", "1s+", -1.5500083192e+00},
                                  {"3d-", "3d-", 1.1269636147e+02}};

    const IO::InputBlock options{"", "k=2;"};
    auto h = DiracOperator::generate("Ek", options, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 2);
    REQUIRE(h->parity() == 1);
    REQUIRE(h->name() == "E2");
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hab = h->reducedME(Fa, Fb);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
    }
  }

  //--------------------------------------------------------------------
  SECTION("M1") {
    std::cout << "M1\n";
    const auto data = std::vector{std::tuple{"1s+", "1s+", 2.4494462627e+00},
                                  {"1s+", "3d-", 1.0524536312e-06},
                                  {"2p-", "2p-", 8.1648571087e-01},
                                  {"2p-", "2p+", -1.1546966950e+00},
                                  {"2p+", "2p-", 1.1546966950e+00},
                                  {"2p+", "2p+", 5.1639502961e+00},
                                  {"3d-", "1s+", -1.0524536312e-06},
                                  {"3d-", "3d-", 3.0983744553e+00}};

    const IO::InputBlock options{""};
    auto h = DiracOperator::generate("M1", options, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 1);
    REQUIRE(h->parity() == 1);
    REQUIRE(h->name() == "M1");
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      h->updateFrequency(std::abs(Fa.en() - Fb.en()));
      const auto hab = h->reducedME(Fa, Fb);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
    }
  }

  //--------------------------------------------------------------------
  SECTION("hfs") {
    // test data generated with "old" mu = 2.582025
    std::cout << "hfs\n";
    const auto data0 = std::vector{std::tuple{"1s+", "1s+", 2.2989921474e+02},
                                   {"1s+", "3d-", -1.2365855864e+00},
                                   {"2p-", "2p-", 9.5793869234e+00},
                                   {"2p-", "2p+", 1.6932864919e+00},
                                   {"2p+", "2p-", -1.6932864919e+00},
                                   {"2p+", "2p+", 6.0579985492e+00},
                                   {"3d-", "1s+", 1.2365855864e+00},
                                   {"3d-", "3d-", 1.0769837860e+00}};

    const auto dataB = std::vector{std::tuple{"1s+", "1s+", 2.2985924212e+02},
                                   {"1s+", "3d-", -1.2365855713e+00},
                                   {"2p-", "2p-", 9.5793868143e+00},
                                   {"2p-", "2p+", 1.6932864709e+00},
                                   {"2p+", "2p-", -1.6932864709e+00},
                                   {"2p+", "2p+", 6.0579985492e+00},
                                   {"3d-", "1s+", 1.2365855713e+00},
                                   {"3d-", "3d-", 1.0769837860e+00}};

    const auto dataSP = std::vector{std::tuple{"1s+", "1s+", 2.2987811905e+02},
                                    {"1s+", "3d-", -1.2365855826e+00},
                                    {"2p-", "2p-", 9.5793868821e+00},
                                    {"2p-", "2p+", 1.6932864866e+00},
                                    {"2p+", "2p-", -1.6932864866e+00},
                                    {"2p+", "2p+", 6.0579985492e+00},
                                    {"3d-", "1s+", 1.2365855826e+00},
                                    {"3d-", "3d-", 1.0769837860e+00}};

    const auto dataBW_SP =
        std::vector{std::tuple{"1s+", "1s+", -9.1760577658e-05},
                    {"1s+", "3d-", -3.0880838988e-09},
                    {"2p-", "2p-", -4.3123067284e-09},
                    {"2p-", "2p+", -3.1267423824e-09},
                    {"2p+", "2p-", -3.1267423824e-09},
                    {"2p+", "2p+", -8.7967510637e-15},
                    {"3d-", "1s+", -3.0880838988e-09},
                    {"3d-", "3d-", -0.0000000000e+00}};

    // test data generated with "old" mu = 2.582025
    const IO::InputBlock options{""};
    auto h0 = DiracOperator::generate(
        "hfs", {"hfs", "F(r)=pointlike; mu=2.582025;"}, wf);
    auto hB =
        DiracOperator::generate("hfs", {"hfs", "F(r)=Ball; mu=2.582025;"}, wf);
    auto hS = DiracOperator::generate(
        "hfs", {"hfs", "F(r)=SingleParticle; mu=2.582025;"}, wf);
    auto h0_au = DiracOperator::generate(
        "hfs", {"hfs", "F(r)=pointlike; units=au; mu=2.582025;"}, wf);

    REQUIRE(h0->get_d_order() == 0);
    REQUIRE(h0->imaginaryQ() == false);
    REQUIRE(h0->rank() == 1);
    REQUIRE(h0->parity() == 1);
    REQUIRE(h0->name() == "hfs1");
    REQUIRE(!h0->units().empty());

    for (std::size_t i = 0; i < data0.size(); ++i) {
      const auto [a, b, a0] = data0.at(i);
      const auto ab = std::get<2>(dataB.at(i));
      const auto as = std::get<2>(dataSP.at(i));
      const auto bw = std::get<2>(dataBW_SP.at(i));
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto A0 = h0->reducedME(Fa, Fb);
      const auto AB = hB->reducedME(Fa, Fb);
      const auto AS = hS->reducedME(Fa, Fb);
      const auto A0_2 = h0_au->reducedME(Fa, Fb) * PhysConst::muN_CGS_MHz;
      const auto BW = (AS - A0) / A0;
      REQUIRE(std::abs(A0 - a0) < 1.0e-6);
      REQUIRE(std::abs(AB - ab) < 1.0e-6);
      REQUIRE(std::abs(AS - as) < 1.0e-6);
      REQUIRE(std::abs(BW - bw) < 1.0e-8);
      // REQUIRE(std::abs(A0 - A0_2) < 1.0e-14);
      // Above fails on M1 mac? Different floating pt? It was very tight anyway
      REQUIRE(std::abs(A0 - A0_2) < 1.0e-12);
    }
  }

  //--------------------------------------------------------------------
  SECTION("hfs(2)") {
    std::cout << "hfs(2)\n";
    const auto data = std::vector{std::tuple{"1s+", "3d-", 2.1316106092e+00},
                                  {"2p-", "2p+", 8.7571740357e+00},
                                  {"2p+", "2p-", -8.7571740357e+00},
                                  {"2p+", "2p+", 8.7568728570e+00},
                                  {"3d-", "1s+", -2.1316106092e+00},
                                  {"3d-", "3d-", 5.1893316793e-01}};

    const IO::InputBlock options{""};
    auto h =
        DiracOperator::generate("hfs", {"hfs", "k=2; F(r)=pointlike;"}, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 2);
    REQUIRE(h->parity() == 1);
    REQUIRE(h->name() == "hfs2");
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hab = h->reducedME(Fa, Fb);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
    }
  }

  //--------------------------------------------------------------------
  SECTION("r") {
    std::cout << "r\n";
    const auto data = std::vector{std::tuple{"1s+", "1s+", 2.1212826887e+00},
                                  {"2p-", "2p-", 7.0709030723e+00},
                                  {"2p+", "2p+", 9.9999467487e+00},
                                  {"3d-", "3d-", 2.0999852080e+01}};

    auto h = DiracOperator::generate("r", {}, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 0);
    REQUIRE(h->parity() == 1);
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hab = h->reducedME(Fa, Fb);
      const auto Rab1 = h->radialIntegral(Fa, Fb);
      const auto Rab2 = Fa * (wf.grid().r() * Fb);
      const auto Rab3 = hab * h->rme3js(Fa.twoj(), Fb.twoj(), 1, 0);
      const auto Rab4 = hab * h->rme3js(Fa.twoj(), Fb.twoj(),
                                        std::min(Fa.twoj(), Fb.twoj()), 0);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
      REQUIRE(std::abs(Rab1 - Rab2) < 1.0e-10);
      REQUIRE(std::abs(Rab1 - Rab3) < 1.0e-10);
      REQUIRE(std::abs(Rab1 - Rab4) < 1.0e-10);
    }
  }

  //--------------------------------------------------------------------
  SECTION("pnc") {
    std::cout << "pnc\n";
    const auto data = std::vector{std::tuple{"1s+", "2p-", 3.9534392752e-07},
                                  {"2p-", "1s+", -3.9534392752e-07},
                                  {"2p+", "3d-", 3.4613851942e-17},
                                  {"3d-", "2p+", -3.4613851942e-17}};

    const IO::InputBlock options{""};
    auto h =
        DiracOperator::generate("pnc", {"pnc", "c=5.67073; t=2.3; N=-1;"}, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == true);
    REQUIRE(h->rank() == 0);
    REQUIRE(h->parity() == -1);
    REQUIRE(h->name() == "pnc-nsi");
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto &Fb = *wf.getState(b);
      const auto hab = h->reducedME(Fa, Fb);
      REQUIRE(std::abs((hab - me) / me) < 1.0e-5);
    }
  }

  //--------------------------------------------------------------------
  SECTION("dr") {
    std::cout << "dr\n";

    auto h = DiracOperator::generate("dr", {}, wf);

    REQUIRE(h->get_d_order() == 1);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 0);
    REQUIRE(h->parity() == 1);
    REQUIRE(!h->units().empty());

    for (const auto &Fa : wf.valence()) {
      const auto hab = h->reducedME(Fa, Fa);
      REQUIRE(std::abs(hab) < 1.0e-10);
    }
  }

  //--------------------------------------------------------------------
  SECTION("p") {
    std::cout << "p\n";

    std::cout << "Test p operator (Li)\n";
    Wavefunction wf2({2000, 1.0e-5, 80.0, 20.0, "loglinear"},
                     {"Li", 0, "pointlike"});
    wf2.solve_core("HartreeFock", 0.0, "[He]");
    wf2.solve_valence("2sp");

    const auto test_data = std::vector{
        std::pair{"2s+", 0.0}, {"2p-", -0.04163}, {"2p+", -0.04162}};

    const auto p = DiracOperator::generate("p", {}, wf);
    for (const auto &[state, expected] : test_data) {
      const auto v = wf2.getState(state);
      REQUIRE(v != nullptr);
      double t = 0.0;
      for (const auto &c : wf2.core()) {
        auto pvc = p->reducedME(*v, c);
        t += -1.0 * pvc * pvc;
      }
      t /= v->twojp1();
      std::cout << *v << " " << t << " [" << expected << "]\n";
      REQUIRE(t == Approx(expected).margin(0.000005));
    }

    //
  }

  //--------------------------------------------------------------------
  SECTION("QED") {
    std::cout << "QED (Vrad)\n";
    const auto data = std::vector{std::tuple{"1s+", "1s+", 4.4297913070e-05},
                                  {"2p-", "2p-", -1.4054241678e-07},
                                  {"2p+", "2p+", 6.3439374927e-08},
                                  {"3d-", "3d-", -6.9277166579e-09}};

    // Uehling only (fast)
    auto h = DiracOperator::generate(
        "Vrad",
        {"Vrad", "readwrite=false; Ueh=1.0; SE_h=1.0; SE_l=1.0; SE_m=1.0"}, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 0);
    REQUIRE(h->parity() == 1);
    REQUIRE(h->name() == "Vrad");
    REQUIRE(!h->units().empty());

    for (auto &[a, b, me] : data) {
      const auto &Fa = *wf.getState(a);
      const auto hab = h->radialIntegral(Fa, Fa);
      REQUIRE(std::abs(hab - me) < 1.0e-6);
    }
  }

  //--------------------------------------------------------------------
  SECTION("jL") {
    std::cout << "jL\n";

    const auto qgrid = Grid(0.01, 100.0, 10, GridType::logarithmic);

    std::size_t in_max_l = 3;
    auto jl = DiracOperator::jL(wf.grid(), qgrid, in_max_l);
    auto g0jl = DiracOperator::g0jL(wf.grid(), qgrid, in_max_l);
    auto ig5jL = DiracOperator::ig5jL(wf.grid(), qgrid, in_max_l);
    // auto ig0g5jL = DiracOperator::ig0g5jL(wf.grid(), qgrid, in_max_l);
    // constructing from existing operator: copies table
    // (should only see "Filling jL lookup table:" 3 times)
    auto ig0g5jL = DiracOperator::ig0g5jL(jl);

    REQUIRE(&jl.q_grid() == &qgrid);
    REQUIRE(&jl.r_grid() == &wf.grid());

    // We never call set_L_q on these ones!
    const DiracOperator::jL jl_0 = jl;
    const DiracOperator::g0jL g0jl_0 = jl;
    const DiracOperator::ig5jL ig5jL_0 = jl;
    const DiracOperator::ig0g5jL ig0g5jL_0 = jl;

    REQUIRE(jl.get_d_order() == 0);
    REQUIRE(jl.imaginaryQ() == false);
    REQUIRE(jl.max_L() == in_max_l);

    // test copy constructor
    auto jl_2 = jl;

    auto C = [](int ka, int kb, std::size_t l) {
      auto c = Angular::Ck_kk(int(l), ka, kb);
      return double(2 * l + 1) * c * c;
    };
    auto D = [C](int ka, int kb, std::size_t l) { return C(ka, -kb, l); };

    // For each L and q, test against directly calculating ME
    for (auto l = 0ul; l <= jl.max_L(); ++l) {
      for (auto &q : qgrid.r()) {
        jl.set_L_q(l, q);
        jl_2.set_L_q(l, q);
        // std::cout << l << " " << q << "\n";
        REQUIRE(jl.rank() == int(l));
        REQUIRE(jl.L() == l);
        const auto jl2 =
            SphericalBessel::fillBesselVec_kr(int(l), q, wf.grid().r());

        for (const auto &a : wf.valence()) {
          for (const auto &b : wf.valence()) {
            REQUIRE(jl.isZero(a, b) == jl.is_zero(a, b, l));
            if (jl.isZero(a, b))
              continue;
            //test radial integral vs. diract calculation
            auto ri1 = jl.radialIntegral(a, b);

            const auto ri2 = a * (jl2 * b);
            REQUIRE(ri1 == Approx(ri2));

            // test square of RME (also test copy of jl_2)
            auto rme1 = double(2 * l + 1) * std::pow(jl_2.reducedME(a, b), 2);
            auto rme2 = ri2 * ri2 * C(a.kappa(), b.kappa(), l);
            REQUIRE(rme1 == Approx(rme2));

            // test direct
            REQUIRE(jl.reducedME(a, b) == Approx(jl_0.rme(a, b, l, q)));
          }
        }
      }
    }

    // test copy of derived:
    auto g0jl_2 = g0jl;

    // For each L and q, test against directly calculating ME
    // For the g0, g5, and g0g5 versions:
    for (auto l = 0ul; l <= jl.max_L(); ++l) {
      for (auto &q : qgrid.r()) {
        g0jl.set_L_q(l, q);
        g0jl_2.set_L_q(l, q);
        ig5jL.set_L_q(l, q);
        ig0g5jL.set_L_q(l, q);
        REQUIRE(g0jl.L() == l);
        REQUIRE(g0jl_2.L() == l);
        REQUIRE(ig5jL.L() == l);
        REQUIRE(ig0g5jL.L() == l);

        const auto J_L =
            SphericalBessel::fillBesselVec_kr(int(l), q, wf.grid().r());

        for (const auto &a : wf.valence()) {
          for (const auto &b : wf.valence()) {
            REQUIRE(g0jl.isZero(a, b) == g0jl_0.is_zero(a, b, l));
            REQUIRE(ig0g5jL.isZero(a, b) == ig0g5jL_0.is_zero(a, b, l));
            REQUIRE(ig5jL.isZero(a, b) == ig5jL_0.is_zero(a, b, l));

            // test copy:
            REQUIRE(g0jl.reducedME(a, b) == Approx(g0jl_2.reducedME(a, b)));

            const auto &gr = wf.grid();
            const auto Rf =
                NumCalc::integrate(gr.du(), 0, 0, a.f(), b.f(), J_L, gr.drdu());
            const auto Rg =
                NumCalc::integrate(gr.du(), 0, 0, a.g(), b.g(), J_L, gr.drdu());
            const auto Rfg =
                NumCalc::integrate(gr.du(), 0, 0, a.f(), b.g(), J_L, gr.drdu());
            const auto Rgf =
                NumCalc::integrate(gr.du(), 0, 0, a.g(), b.f(), J_L, gr.drdu());

            const auto rme2_g0 =
                C(a.kappa(), b.kappa(), l) * (Rf * Rf - 2 * Rf * Rg + Rg * Rg);
            const auto rme2_g0g5 = D(a.kappa(), b.kappa(), l) *
                                   (Rfg * Rfg + 2 * Rfg * Rgf + Rgf * Rgf);
            const auto rme2_g5 = D(a.kappa(), b.kappa(), l) *
                                 (Rfg * Rfg - 2 * Rfg * Rgf + Rgf * Rgf);

            const auto rme1_g0 =
                double(2 * l + 1) * std::pow(g0jl.reducedME(a, b), 2);
            const auto rme1_g0g5 =
                double(2 * l + 1) * std::pow(ig0g5jL.reducedME(a, b), 2);
            const auto rme1_g5 =
                double(2 * l + 1) * std::pow(ig5jL.reducedME(a, b), 2);

            // use 'margin' here, since some MEs are zero (don't skip, since some have different selection rules)

            REQUIRE(rme2_g0 == Approx(rme1_g0).margin(1.0e-12));
            REQUIRE(rme2_g0g5 == Approx(rme1_g0g5).margin(1.0e-12));
            REQUIRE(rme2_g5 == Approx(rme1_g5).margin(1.0e-12));

            // test direct
            REQUIRE(g0jl.reducedME(a, b) ==
                    Approx(g0jl_0.rme(a, b, l, q)).margin(1.0e-12));
            REQUIRE(ig0g5jL.reducedME(a, b) ==
                    Approx(ig0g5jL_0.rme(a, b, l, q)).margin(1.0e-12));
            REQUIRE(ig5jL.reducedME(a, b) ==
                    Approx(ig5jL_0.rme(a, b, l, q)).margin(1.0e-12));
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------
  SECTION("SpinorMatrix") {
    std::cout << "SpinorMatrix\n";
    DiracOperator::SpinorMatrix a;

    for (const auto &Fv : wf.valence()) {
      double zero = Fv * (a * Fv);
      REQUIRE(zero == Approx(0.0));
    }

    DiracOperator::SpinorMatrix b{1.0, 0.0, 0.0, 1.0};
    a = b;
    for (const auto &Fv : wf.valence()) {
      double lhs = Fv * (a * Fv);
      double rhs = Fv * (Fv);
      REQUIRE(lhs == Approx(rhs));
    }

    auto c = 2.0 * b;
    c *= 1.0;
    auto d = a + b;
    for (const auto &Fv : wf.valence()) {
      double lhs = Fv * (c * Fv);
      double rhs = Fv * (d * Fv);
      REQUIRE(lhs == Approx(rhs));
    }

    DiracOperator::SpinorMatrix e{0.0, 1.0, 1.0, 0.0};
    for (const auto &Fv : wf.valence()) {
      double lhs = Fv * (e * Fv);

      auto mFv = Fv;
      mFv.f() = Fv.g();
      mFv.g() = Fv.f();
      double rhs = Fv * (mFv);

      REQUIRE(lhs == Approx(rhs));
    }
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// Used to generate test data (with 8000 points)
// for (auto &a : wf.valence()) {
//   for (auto &b : wf.valence()) {
//     if (h->isZero(a, b))
//       continue;
//     auto me = h->reducedME(a, b);
//     std::cout << "{\"" << a << "\", \"" << b << "\", ";
//     printf("%.10e},\n", me);
//   }
// }
