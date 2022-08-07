#include "DiracOperator.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <utility>

// {"E1", &generate_E1},
// {"Ek", &generate_Ek},
// {"M1", &generate_M1},
// {"hfs", &generate_hfsA},
// {"hfsK", &generate_hfsK},
// {"r", &generate_r},
// {"pnc", &generate_pnc},
// {"Hrad", &generate_Hrad}

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
    auto hv = DiracOperator::generate("E1", {"", "gauge=vform;"}, wf);
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
    std::cout << "hfs\n";
    const auto data0 = std::vector{std::tuple{"1s+", "1s+", -2.2989921474e+02},
                                   {"1s+", "3d-", 1.2365855864e+00},
                                   {"2p-", "2p-", -9.5793869234e+00},
                                   {"2p-", "2p+", -1.6932864919e+00},
                                   {"2p+", "2p-", 1.6932864919e+00},
                                   {"2p+", "2p+", -6.0579985492e+00},
                                   {"3d-", "1s+", -1.2365855864e+00},
                                   {"3d-", "3d-", -1.0769837860e+00}};

    const auto dataB = std::vector{std::tuple{"1s+", "1s+", -2.2985924212e+02},
                                   {"1s+", "3d-", 1.2365855713e+00},
                                   {"2p-", "2p-", -9.5793868143e+00},
                                   {"2p-", "2p+", -1.6932864709e+00},
                                   {"2p+", "2p-", 1.6932864709e+00},
                                   {"2p+", "2p+", -6.0579985492e+00},
                                   {"3d-", "1s+", -1.2365855713e+00},
                                   {"3d-", "3d-", -1.0769837860e+00}};

    const auto dataSP = std::vector{std::tuple{"1s+", "1s+", -2.2987811905e+02},
                                    {"1s+", "3d-", 1.2365855826e+00},
                                    {"2p-", "2p-", -9.5793868821e+00},
                                    {"2p-", "2p+", -1.6932864866e+00},
                                    {"2p+", "2p-", 1.6932864866e+00},
                                    {"2p+", "2p+", -6.0579985492e+00},
                                    {"3d-", "1s+", -1.2365855826e+00},
                                    {"3d-", "3d-", -1.0769837860e+00}};

    const auto dataBW_SP =
        std::vector{std::tuple{"1s+", "1s+", -9.1760577658e-05},
                    {"1s+", "3d-", -3.0880838988e-09},
                    {"2p-", "2p-", -4.3123067284e-09},
                    {"2p-", "2p+", -3.1267423824e-09},
                    {"2p+", "2p-", -3.1267423824e-09},
                    {"2p+", "2p+", -8.7967510637e-15},
                    {"3d-", "1s+", -3.0880838988e-09},
                    {"3d-", "3d-", -0.0000000000e+00}};

    const IO::InputBlock options{""};
    auto h0 = DiracOperator::generate("hfs", {"hfs", "F(r)=pointlike;"}, wf);
    auto hB = DiracOperator::generate("hfs", {"hfs", "F(r)=Ball;"}, wf);
    auto hS = DiracOperator::generate("hfs", {"hfs", "F(r)=VolotkaBW;"}, wf);

    REQUIRE(h0->get_d_order() == 0);
    REQUIRE(h0->imaginaryQ() == false);
    REQUIRE(h0->rank() == 1);
    REQUIRE(h0->parity() == 1);
    REQUIRE(h0->name() == "hfs");
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
      const auto BW = (AS - A0) / A0;
      REQUIRE(std::abs(A0 - a0) < 1.0e-6);
      REQUIRE(std::abs(AB - ab) < 1.0e-6);
      REQUIRE(std::abs(AS - as) < 1.0e-6);
      REQUIRE(std::abs(BW - bw) < 1.0e-8);
    }
  }

  //--------------------------------------------------------------------
  SECTION("hfsK(2)") {
    std::cout << "hfsK(2)\n";
    const auto data = std::vector{std::tuple{"1s+", "3d-", 2.1316106092e+00},
                                  {"2p-", "2p+", 8.7571740357e+00},
                                  {"2p+", "2p-", -8.7571740357e+00},
                                  {"2p+", "2p+", 8.7568728570e+00},
                                  {"3d-", "1s+", -2.1316106092e+00},
                                  {"3d-", "3d-", 5.1893316793e-01}};

    const IO::InputBlock options{""};
    auto h =
        DiracOperator::generate("hfsK", {"hfsK", "k=2; F(r)=pointlike;"}, wf);

    REQUIRE(h->get_d_order() == 0);
    REQUIRE(h->imaginaryQ() == false);
    REQUIRE(h->rank() == 2);
    REQUIRE(h->parity() == 1);
    REQUIRE(h->name() == "hfs(2)");
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
    auto h = DiracOperator::generate("pnc", {"pnc", "c=5.67073; t=2.3;"}, wf);

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
}

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
