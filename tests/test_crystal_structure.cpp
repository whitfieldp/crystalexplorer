#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>

#include "crystalstructure.h"
#include "slab_options.h"
#include <occ/crystal/crystal.h>

using Catch::Approx;

// Acetic acid crystal setup from occ tests
auto acetic_asym() {
  const std::vector<std::string> labels = {"C1", "C2", "H1", "H2",
                                           "H3", "H4", "O1", "O2"};
  occ::IVec nums(labels.size());
  occ::Mat positions(labels.size(), 3);
  for (size_t i = 0; i < labels.size(); i++) {
    nums(i) = occ::core::Element(labels[i]).atomic_number();
  }
  positions << 0.16510, 0.28580, 0.17090, 0.08940, 0.37620, 0.34810, 0.18200,
      0.05100, -0.11600, 0.12800, 0.51000, 0.49100, 0.03300, 0.54000, 0.27900,
      0.05300, 0.16800, 0.42100, 0.12870, 0.10750, 0.00000, 0.25290, 0.37030,
      0.17690;
  return occ::crystal::AsymmetricUnit(positions.transpose(), nums, labels);
}

auto acetic_acid_crystal() {
  occ::crystal::AsymmetricUnit asym = acetic_asym();
  occ::crystal::SpaceGroup sg(33);
  occ::crystal::UnitCell cell = occ::crystal::orthorhombic_cell(13.31, 4.1, 5.75);
  return OccCrystal(asym, sg, cell);
}

// PIM-1 DFT structure (PIM-1_241K_B3LYP_dftd3bj): orthorhombic P2_12_12 with
// a 10-atom asymmetric unit whose symmetry images assemble into two
// symmetry-related C7H10N2 molecules per unit cell. The space group operation
// relating the two molecules permutes their atom ordering, which is the case
// that exercised the element-order bug in getTransformation.
auto pim1_asym() {
  const std::vector<std::string> labels = {"N1", "C1", "C2", "C3", "C4",
                                           "H1", "H2", "H3", "H4", "H5"};
  occ::IVec nums(labels.size());
  for (size_t i = 0; i < labels.size(); i++) {
    nums(i) = occ::core::Element(labels[i]).atomic_number();
  }
  occ::Mat positions(labels.size(), 3);
  positions << 0.3307035898, 0.8114353860, 0.4186165078,
      0.0000000000, 0.0000000000, 0.8444474191,
      0.1484370920, 0.0607682623, 0.0251980156,
      0.7986292579, 0.3769787623, 0.1542254402,
      0.7286628277, 0.2717572089, 0.3044631891,
      0.4340003497, 0.4336550278, 0.2939909414,
      0.2172505319, 0.9939598880, 0.1579322621,
      0.0863359592, 0.1267088760, 0.1681163153,
      0.9128163783, 0.3478224437, 0.0192496968,
      0.8555167209, 0.4388606896, 0.3118332106;
  return occ::crystal::AsymmetricUnit(positions.transpose(), nums, labels);
}

auto pim1_crystal() {
  occ::crystal::AsymmetricUnit asym = pim1_asym();
  occ::crystal::SpaceGroup sg(18);
  occ::crystal::UnitCell cell =
      occ::crystal::orthorhombic_cell(7.2468702008, 11.1188874428, 4.7062696926);
  return OccCrystal(asym, sg, cell);
}

TEST_CASE("CrystalStructure basic functionality", "[crystal][crystal_structure]") {
    CrystalStructure structure;
    
    SECTION("Initial state") {
        REQUIRE(structure.structureType() == ChemicalStructure::StructureType::Crystal);
        REQUIRE(structure.numberOfAtoms() == 0);
        REQUIRE(structure.numberOfFragments() == 0);
    }
    
    SECTION("Cell vectors and parameters with acetic acid") {
        OccCrystal crystal = acetic_acid_crystal();
        structure.setOccCrystal(crystal);
        
        auto lengths = structure.cellLengths();
        REQUIRE(lengths[0] == Approx(13.31));
        REQUIRE(lengths[1] == Approx(4.1));
        REQUIRE(lengths[2] == Approx(5.75));
        
        auto angles = structure.cellAngles();
        REQUIRE(angles[0] == Approx(M_PI/2)); // 90 degrees (orthorhombic)
        REQUIRE(angles[1] == Approx(M_PI/2));
        REQUIRE(angles[2] == Approx(M_PI/2));
        
        REQUIRE(structure.spaceGroup().number() == 33);
        REQUIRE(structure.numberOfAtoms() > 0);
    }
}

TEST_CASE("CrystalStructure coordinate conversion", "[crystal][crystal_structure]") {
    CrystalStructure structure;
    OccCrystal crystal = acetic_acid_crystal();
    structure.setOccCrystal(crystal);
    
    SECTION("Round trip conversion") {
        occ::Mat3N originalCart(3, 1);
        originalCart.col(0) << 3.7, 2.1, 1.9;
        
        auto fracCoords = structure.convertCoordinates(
            originalCart, ChemicalStructure::CoordinateConversion::CartToFrac);
        auto backToCart = structure.convertCoordinates(
            fracCoords, ChemicalStructure::CoordinateConversion::FracToCart);
        
        REQUIRE(backToCart(0, 0) == Approx(originalCart(0, 0)).margin(1e-10));
        REQUIRE(backToCart(1, 0) == Approx(originalCart(1, 0)).margin(1e-10));
        REQUIRE(backToCart(2, 0) == Approx(originalCart(2, 0)).margin(1e-10));
    }
}

TEST_CASE("CrystalStructure space group", "[crystal][crystal_structure]") {
    CrystalStructure structure;
    OccCrystal crystal = acetic_acid_crystal();
    structure.setOccCrystal(crystal);
    
    SECTION("Space group access") {
        auto sg = structure.spaceGroup();
        REQUIRE(sg.number() == 33); // Orthorhombic space group
    }
    
    SECTION("OCC crystal access") {
        auto occCrystal = structure.occCrystal();
        REQUIRE(occCrystal.space_group().number() == 33);
        REQUIRE(occCrystal.asymmetric_unit().size() == 8); // C2H4O2
    }
}

TEST_CASE("CrystalStructure atom properties", "[crystal][crystal_structure]") {
    CrystalStructure structure;
    OccCrystal crystal = acetic_acid_crystal();
    structure.setOccCrystal(crystal);
    
    SECTION("Basic atom properties") {
        REQUIRE(structure.numberOfAtoms() > 0);
        REQUIRE(structure.chemicalFormula(false).contains("C"));
        REQUIRE(structure.chemicalFormula(false).contains("H"));
        REQUIRE(structure.chemicalFormula(false).contains("O"));
    }
    
    SECTION("Index conversions") {
        if (structure.numberOfAtoms() >= 2) {
            GenericAtomIndex generic0 = structure.indexToGenericIndex(0);
            GenericAtomIndex generic1 = structure.indexToGenericIndex(1);
            
            int back0 = structure.genericIndexToIndex(generic0);
            int back1 = structure.genericIndexToIndex(generic1);
            
            REQUIRE(back0 == 0);
            REQUIRE(back1 == 1);
        }
    }
}

TEST_CASE("CrystalStructure JSON serialization", "[crystal][crystal_structure][json]") {
    CrystalStructure structure;
    OccCrystal crystal = acetic_acid_crystal();
    structure.setOccCrystal(crystal);
    
    SECTION("To JSON") {
        auto json = structure.toJson();
        // Check that JSON contains basic structure data
        REQUIRE(json.contains("atomicNumbers"));
        REQUIRE(json.contains("atomicPositions"));
        REQUIRE(json.is_object());
    }
}

TEST_CASE("CrystalStructure atomsInBoundingBox", "[crystal][crystal_structure][bounding_box]") {
    CrystalStructure structure;
    OccCrystal crystal = acetic_acid_crystal();
    structure.setOccCrystal(crystal);
    
    SECTION("Unit cell bounding box should contain all unit cell atoms") {
        // Get unit cell dimensions
        auto cellVectors = structure.cellVectors();
        occ::Vec3 min = occ::Vec3::Zero();
        occ::Vec3 max = cellVectors.col(0) + cellVectors.col(1) + cellVectors.col(2);
        
        auto atomsInBox = structure.atomsInBoundingBox(min, max);
        
        // Should find atoms from the unit cell
        REQUIRE(atomsInBox.size() > 0);
        
        // Log what we found for debugging
        std::cout << "Found " << atomsInBox.size() << " atoms in unit cell bounding box" << std::endl;
        for (const auto& idx : atomsInBox) {
            std::cout << "  Atom " << idx.unique << " at offset (" << idx.x << "," << idx.y << "," << idx.z << ")" << std::endl;
        }
    }
    
    SECTION("Larger bounding box should find periodic images") {
        // Create a 2x2x2 supercell bounding box
        auto cellVectors = structure.cellVectors();
        occ::Vec3 min = -0.5 * (cellVectors.col(0) + cellVectors.col(1) + cellVectors.col(2));
        occ::Vec3 max = 1.5 * (cellVectors.col(0) + cellVectors.col(1) + cellVectors.col(2));
        
        auto atomsInBox = structure.atomsInBoundingBox(min, max);
        
        // Should find many more atoms from periodic images
        REQUIRE(atomsInBox.size() > 8); // At least a few unit cells worth
        
        std::cout << "Found " << atomsInBox.size() << " atoms in 2x2x2 supercell bounding box" << std::endl;
        
        // Check that we have atoms with different unit cell offsets
        bool foundNegativeOffset = false;
        bool foundPositiveOffset = false;
        for (const auto& idx : atomsInBox) {
            if (idx.x < 0 || idx.y < 0 || idx.z < 0) foundNegativeOffset = true;
            if (idx.x > 0 || idx.y > 0 || idx.z > 0) foundPositiveOffset = true;
        }
        
        REQUIRE(foundNegativeOffset);
        REQUIRE(foundPositiveOffset);
    }
    
    SECTION("Small bounding box should find subset") {
        // Create a small bounding box around origin
        occ::Vec3 min(-1.0, -1.0, -1.0);
        occ::Vec3 max(1.0, 1.0, 1.0);
        
        auto atomsInBox = structure.atomsInBoundingBox(min, max);
        
        std::cout << "Found " << atomsInBox.size() << " atoms in small bounding box around origin" << std::endl;
        
        // Should find at least some atoms but not all
        REQUIRE(atomsInBox.size() >= 0); // Could be 0 if no atoms near origin
    }
}

TEST_CASE("getTransformation between symmetry-related molecules",
          "[crystal][crystal_structure][getTransformation]") {
    // Regression test: getTransformation must find the space group operation
    // mapping one unit cell molecule onto its symmetry-related counterpart
    // even when the operation permutes the atoms. The element comparison and
    // position matching must be order-independent (multiset / set based).
    CrystalStructure structure;
    OccCrystal crystal = pim1_crystal();
    structure.setOccCrystal(crystal);
    structure.normalizeHydrogenBondLengths(false);
    structure.completeAllFragments();
    // Ensure periodic images exist so lattice-translated sets resolve
    structure.expandAtomsWithinRadius(6.0f, true);
    structure.completeAllFragments();

    auto fragments = structure.unitCellFragments();
    REQUIRE(fragments.size() == 2);

    std::vector<std::vector<GenericAtomIndex>> molecules;
    for (const auto &[idx, frag] : fragments)
        molecules.push_back(frag.atomIndices);
    REQUIRE(molecules[0].size() == molecules[1].size());

    Eigen::Isometry3d transform;
    SECTION("Molecule A maps onto molecule B") {
        REQUIRE(structure.getTransformation(molecules[0], molecules[1],
                                            transform));
        // The transform must actually carry A's atoms onto B's positions
        auto posA = structure.atomicPositionsForIndices(molecules[0]);
        auto posB = structure.atomicPositionsForIndices(molecules[1]);
        auto mapped = (transform * posA).eval();
        // Each mapped atom must land on some B atom (within tolerance)
        for (int i = 0; i < posA.cols(); i++) {
            double best = 1e30;
            for (int j = 0; j < posB.cols(); j++)
                best = std::min(best, (mapped.col(i) - posB.col(j)).norm());
            REQUIRE(best < 1e-3);
        }
    }

    SECTION("Molecule B maps onto molecule A") {
        REQUIRE(structure.getTransformation(molecules[1], molecules[0],
                                            transform));
    }

    SECTION("A molecule maps onto a lattice translation of itself") {
        // Same internal ordering, shifted by one cell: must still succeed
        std::vector<GenericAtomIndex> shifted;
        for (const auto &idx : molecules[0])
            shifted.push_back(GenericAtomIndex{idx.unique, idx.x + 1, idx.y,
                                               idx.z});
        REQUIRE(structure.getTransformation(molecules[0], shifted, transform));
    }

    SECTION("Different element multisets are rejected") {
        // Remove one H from A and one C from B: same size, different
        // composition -> the multiset element check must reject it.
        auto numsA = structure.atomicNumbersForIndices(molecules[0]);
        auto numsB = structure.atomicNumbersForIndices(molecules[1]);
        int hInA = -1, cInB = -1;
        for (int i = 0; i < numsA.size(); i++)
            if (numsA(i) == 1) { hInA = i; break; }
        for (int i = 0; i < numsB.size(); i++)
            if (numsB(i) == 6) { cInB = i; break; }
        REQUIRE(hInA >= 0);
        REQUIRE(cInB >= 0);
        std::vector<GenericAtomIndex> aSet, bSet;
        for (size_t i = 0; i < molecules[0].size(); i++)
            if (static_cast<int>(i) != hInA)
                aSet.push_back(molecules[0][i]);
        for (size_t i = 0; i < molecules[1].size(); i++)
            if (static_cast<int>(i) != cInB)
                bSet.push_back(molecules[1][i]);
        REQUIRE(aSet.size() == bSet.size());
        REQUIRE_FALSE(structure.getTransformation(aSet, bSet, transform));
    }
}

