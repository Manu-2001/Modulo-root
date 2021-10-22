#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <iostream>
#include <memory>

#include "doctest.h"
#include "impulse.hpp"
#include "particle.hpp"
#include "particletype.hpp"
#include "resonancetype.hpp"

std::unique_ptr<ParticleType> Create(bool);

TEST_CASE("Impulse") {
  {
    Impulse<float> P1{5.F, 9.F, 7.F};
    Impulse<float> P2{2.F, 58.F, 31.5F};
    Impulse<float> P3 = -P1;

    CHECK(P3 == Impulse<float>{-5.F, -9.F, -7.F});
    CHECK(P2 + P1 == Impulse<float>{7.F, 67.F, 38.5F});
    CHECK(P2 + P1 == P1 + P2);
    CHECK(P2 - P1 == Impulse<float>{-3.F, 49.F, 24.5F});
    CHECK(P2 * 5.F == Impulse<float>{10.F, 290.F, 157.5F});
    CHECK(P2 * 5. == Impulse<float>{10.F, 290.F, 157.5F});
    CHECK(P1 + P2 + P3 == P2);
    CHECK(P1.Norm() == (float)sqrt(155));
    CHECK(P1.Norm() == P3.Norm());

    P2 *= 5.;
    CHECK(P2 == Impulse<float>{10.F, 290.F, 157.5F});
  }
}

TEST_CASE("Particle") {
  {
    // test 1
    ParticleType p("Proton", 1.6, 1);
    ResonanceType z("Res", 1.854, -1, .25);
    std::vector<std::unique_ptr<ParticleType>> ParticleVector{};

    ParticleVector.push_back(Create(1));
    ParticleVector.push_back(Create(0));
    ParticleVector.push_back(Create(1));
    ParticleVector.push_back(Create(0));

    std::cout << "*** TESTING ISTANCE ***\n";

    p.Print();
    z.Print();

    std::cout << "\n*** TESTING DYNAMIC POLIMORFISM ***\n";

    for (auto const& Particle : ParticleVector) {
      Particle->Print();
    }
  }

  {
    // test 2
    CHECK_THROWS(ParticleType{"Particle", -5, -1});
    CHECK_THROWS(ResonanceType{"Res", -58, +1, 10E-15});
    CHECK_THROWS(ResonanceType{"p+", 1.6, 0, -10E-15});
    CHECK_THROWS(ResonanceType{"Res", -58, +1, -10E-15});
  }

  {
    // test 3
    std::cout << "\n*** TESTING PARTICLE ***";
    Particle::AddParticleType("electron", 9.11, -1);
    Particle::AddParticleType("proton", 1.67, +1);
    Particle::AddParticleType("pione+", 2.14, +1, 3.2 * 10E-23);
    Particle::AddParticleType("pione-", 2.14, +1, 3.2 * 10E-23);

    std::cout<<"\nTest costructor error:";
    Particle try_("ndef");
    Particle electron1("electron", Impulse<double>{1.2, -2.5, 8.75});
    Particle proton1("proton");
    Particle pion1("pione-", Impulse<double>{1.2, -2.5, 455.2});

    Particle rdef("electron", Impulse<double>{1.2, -2.5, 8.75});

    CHECK(rdef.GetTypeIndex() == 0);
    rdef.SetTypeIndex("pione-");
    CHECK(rdef.GetTypeIndex() == 3);

    CHECK_THROWS(rdef.SetTypeIndex("neutron"));
    CHECK_THROWS(rdef.SetTypeIndex("notdef"));
    CHECK_THROWS(rdef.SetTypeIndex(15));
    CHECK_THROWS(rdef.SetTypeIndex(9));

    rdef.SetMomentum(Impulse<double>{});
    CHECK(rdef.GetMomentum() == Impulse<double>{});

    CHECK(electron1.GetTypeIndex() == 0);
    CHECK(proton1.GetTypeIndex() == 1);
    CHECK(pion1.GetTypeIndex() == 3);

    Particle::PrintParticleType();

    std::cout << '\n';
    electron1.Print();
    proton1.Print();
    pion1.Print();
  }
}

std::unique_ptr<ParticleType> Create(bool const x) {
  if (x) {
    return std::make_unique<ParticleType>("Proton", 1.6, 1);
  }
  return std::make_unique<ResonanceType>("Res", 1.854, -1, .25);
}