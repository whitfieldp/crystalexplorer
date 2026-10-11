// Headless reproduction of the pair-energy workflow for PIM-1_241K:
// import CIF -> complete fragments -> select molecule -> create b3lyp/def2-tzvpd
// wavefunction -> expand to 6 A and complete fragments -> set up pair
// interaction calculation using the previously created wavefunction.
//
// This exercises the same code paths as Crystalx::calculatePairEnergies and
// Crystalx::calculatePairEnergiesWithExistingWavefunctions without launching
// any OCC processes, and reports why individual fragment pairs fail to find
// a matching wavefunction.
#include <QCoreApplication>
#include <QDebug>
#include <cmath>

#include "ciffile.h"
#include "crystalstructure.h"
#include "molecular_wavefunction.h"
#include "pair_energy_parameters.h"
#include "settings.h"

static void printFragments(const ChemicalStructure *s, const QString &label) {
  const auto &frags = s->getFragments();
  qInfo().noquote()
      << QString("---- %1: %2 fragments, %3 atoms total, incomplete=%4 ----")
             .arg(label)
             .arg(frags.size())
             .arg(s->numberOfAtoms())
             .arg(s->hasIncompleteFragments());
  for (const auto &[idx, frag] : frags) {
    qInfo().noquote()
        << QString("  frag[%1,%2,%3,%4] size=%5 asym=%6")
               .arg(idx.u)
               .arg(idx.h)
               .arg(idx.k)
               .arg(idx.l)
               .arg(frag.atomIndices.size())
               .arg(frag.asymmetricFragmentIndex.u);
  }
}

static QString composition(const ChemicalStructure *s,
                           const std::vector<GenericAtomIndex> &idxs) {
  auto nums = s->atomicNumbersForIndices(idxs);
  QMap<int, int> counts;
  for (int i = 0; i < nums.size(); i++)
    counts[nums(i)]++;
  QStringList parts;
  for (auto it = counts.constBegin(); it != counts.constEnd(); ++it)
    parts << QString("%1x%2").arg(it.key()).arg(it.value());
  return parts.join(",");
}

// Diagnose why getTransformation(wfnAtoms, fragAtoms) fails
static QString transformFailureReason(const ChemicalStructure *s,
                                      const std::vector<GenericAtomIndex> &a,
                                      const std::vector<GenericAtomIndex> &b) {
  if (a.size() != b.size())
    return QString("size mismatch (%1 vs %2)").arg(a.size()).arg(b.size());
  auto numsA = s->atomicNumbersForIndices(a);
  auto numsB = s->atomicNumbersForIndices(b);
  auto sa = a;
  auto sb = b;
  std::sort(sa.begin(), sa.end());
  std::sort(sb.begin(), sb.end());
  auto na = s->atomicNumbersForIndices(sa);
  auto nb = s->atomicNumbersForIndices(sb);
  if (!(na.array() == nb.array()).all())
    return "element mismatch";
  return "geometry RMSD > 1e-3";
}

int main(int argc, char **argv) {
  QCoreApplication app(argc, argv);
  qSetMessagePattern("[%{type}] %{message}");

  QString cif = argc > 1
                    ? QString::fromLocal8Bit(argv[1])
                    : "/home/whitfieldp/AI_projects/crystalexplorer/"
                      "example_cif/PIM-1_241K_B3LYP_dftd3bj.cif";

  CifFile cifReader;
  if (!cifReader.readFromFile(cif)) {
    qCritical() << "Failed to read CIF:" << cif;
    return 1;
  }

  const auto &crystal = cifReader.getCrystalStructure(0);
  qInfo() << "OCC crystal: uc atoms =" << crystal.unit_cell_atoms().size()
          << ", asym atoms =" << crystal.asymmetric_unit().size();
  qInfo() << "unit cell molecules:" << crystal.unit_cell_molecules().size()
          << ", symmetry unique molecules:"
          << crystal.symmetry_unique_molecules().size();
  for (const auto &mol : crystal.unit_cell_molecules()) {
    qInfo().noquote()
        << QString("  uc mol %1: asym_mol_idx=%2 size=%3 asym_idx=[%4]")
               .arg(mol.unit_cell_molecule_idx())
               .arg(mol.asymmetric_molecule_idx())
               .arg(mol.size())
               .arg(mol.asymmetric_unit_idx().size());
    QStringList parts;
    for (int i = 0; i < mol.asymmetric_unit_idx().size(); i++)
      parts << QString::number(mol.asymmetric_unit_idx()(i));
    qInfo().noquote() << QString("    asym_unit_idx: %1").arg(parts.join(","));
  }
  const auto &uc_atoms = crystal.unit_cell_atoms();
  QStringList asymMap;
  for (int i = 0; i < uc_atoms.size(); i++)
    asymMap << QString("%1->%2").arg(i).arg(uc_atoms.asym_idx(i));
  qInfo().noquote() << QString("uc atom -> asym atom: %1").arg(asymMap.join(" "));

  // Diagnostic: does any space group operation map UC molecule 0 onto UC
  // molecule 1? Reports, per symop, the maximum (over the reference
  // molecule's atoms) of the minimum fractional distance to the target
  // molecule's atoms (minimum-image over the 27 nearest cells).
  {
    const auto &sg = crystal.space_group();
    qInfo() << "Space group:" << sg.number() << sg.symbol().c_str();
    QStringList ops;
    for (const auto &op : sg.symmetry_operations())
      ops << QString::fromStdString(op.to_string());
    qInfo().noquote() << QString("symops: %1").arg(ops.join(" | "));

    const auto &ucMols = crystal.unit_cell_molecules();
    if (ucMols.size() >= 2) {
      auto fracOf = [&](const occ::core::Molecule &m) {
        occ::Mat3N f(3, static_cast<int>(m.size()));
        for (size_t j = 0; j < m.size(); j++)
          f.col(j) = uc_atoms.frac_pos.col(m.unit_cell_idx()(j)) +
                     m.unit_cell_atom_shift().col(j).cast<double>();
        return f;
      };
      auto wrap = [](const occ::Mat3N &f) {
        return (f.array() - f.array().floor()).matrix();
      };
      occ::Mat3N fracA = fracOf(ucMols[0]);
      occ::Mat3N wB = wrap(fracOf(ucMols[1]));
      for (const auto &op : sg.symmetry_operations()) {
        occ::Mat3N wT = wrap(op.apply(fracA));
        double maxMin = 0.0;
        for (int j = 0; j < wT.cols(); j++) {
          double dmin = 1e9;
          for (int k = 0; k < wB.cols(); k++) {
            occ::Vec3 d = wT.col(j) - wB.col(k);
            for (int ix = -1; ix <= 1; ix++)
              for (int iy = -1; iy <= 1; iy++)
                for (int iz = -1; iz <= 1; iz++)
                  dmin = std::min(
                      dmin, (d - occ::Vec3(ix, iy, iz)).norm());
          }
          maxMin = std::max(maxMin, dmin);
        }
        qInfo().noquote()
            << QString("symop '%1': max-min frac dist A->B = %2")
                   .arg(QString::fromStdString(op.to_string()))
                   .arg(maxMin, 0, 'f', 5);
      }
    }
  }

  auto *structure = new CrystalStructure();
  structure->setOccCrystal(crystal);
  structure->normalizeHydrogenBondLengths(false);
  qInfo() << "=== Step 1: import ===";
  printFragments(structure, "after import");

  qInfo() << "=== Step 2: complete fragments ===";
  structure->completeAllFragments();
  printFragments(structure, "after complete #1");

  qInfo() << "unit cell fragments after complete #1:"
          << structure->unitCellFragments().size();
  for (const auto &[idx, frag] : structure->unitCellFragments()) {
    qInfo().noquote()
        << QString("  uc frag u=%1 size=%2 composition=[%3] asymIdx=%4")
               .arg(idx.u)
               .arg(frag.atomIndices.size())
               .arg(composition(structure, frag.atomIndices))
               .arg(frag.asymmetricFragmentIndex.u);
    auto sorted = frag.atomIndices;
    std::sort(sorted.begin(), sorted.end());
    QStringList parts;
    for (const auto &a : sorted)
      parts << QString("(%1,%2,%3,%4)").arg(a.unique).arg(a.x).arg(a.y).arg(a.z);
    qInfo().noquote() << QString("    atoms: %1").arg(parts.join(" "));
  }

  qInfo() << "=== Step 3: select molecule, create b3lyp/def2-tzvpd wfn ===";
  const auto &fragments0 = structure->getFragments();
  const Fragment &keyFragment0 = fragments0.begin()->second;
  qInfo().noquote()
      << QString("key fragment: %1 atoms, asym index u=%2")
             .arg(keyFragment0.atomIndices.size())
             .arg(keyFragment0.asymmetricFragmentIndex.u);

  std::vector<GenericAtomIndex> selected = keyFragment0.atomIndices;
  structure->setFlagForAtoms(selected, AtomFlag::Selected);

  auto *wfn = new MolecularWavefunction();
  wfn::Parameters wp;
  wp.structure = structure;
  wp.atoms = selected;
  wp.charge = 0;
  wp.multiplicity = 1;
  wp.method = "b3lyp";
  wp.basis = "def2-tzvpd";
  wp.accepted = true;
  wp.program = wfn::Program::Occ;
  wfn->setParameters(wp);
  wfn->setAtomIndices(selected);
  wfn->setObjectName("b3lyp def2-tzvpd [repro]");
  wfn->setParent(structure);
  qInfo().noquote() << QString("wfn composition: [%1]")
                         .arg(composition(structure, selected));

  qInfo() << "=== Step 4: expand to 6 A radius and complete fragments ===";
  structure->expandAtomsWithinRadius(6.0f, true);
  structure->completeAllFragments();
  printFragments(structure, "after expand+complete #2");
  qInfo().noquote()
      << QString("wfn atoms still valid composition: [%1] size=%2")
             .arg(composition(structure, wfn->atomIndices()))
             .arg(wfn->atomIndices().size());

  // Compare unit-cell fragment (unique,hkl) sets with periodic instances
  auto dumpIdxs = [&](const std::vector<GenericAtomIndex> &idxs) {
    auto sorted = idxs;
    std::sort(sorted.begin(), sorted.end());
    QStringList parts;
    for (const auto &a : sorted)
      parts << QString("(%1,%2,%3,%4)").arg(a.unique).arg(a.x).arg(a.y).arg(a.z);
    return parts.join(" ");
  };
  qInfo() << "unit cell fragments after step 4:"
          << structure->unitCellFragments().size();
  for (const auto &[idx, frag] : structure->unitCellFragments()) {
    qInfo().noquote() << QString("  uc u=%1 size=%2: %3").arg(idx.u).arg(frag.atomIndices.size()).arg(dumpIdxs(frag.atomIndices));
  }
  int shown = 0;
  for (const auto &[idx, frag] : structure->getFragments()) {
    if (shown++ >= 4)
      break;
    qInfo().noquote()
        << QString("  periodic frag[%1,%2,%3,%4] asym=%5 size=%6: %7")
               .arg(idx.u)
               .arg(idx.h)
               .arg(idx.k)
               .arg(idx.l)
               .arg(frag.asymmetricFragmentIndex.u)
               .arg(frag.atomIndices.size())
               .arg(dumpIdxs(frag.atomIndices));
  }

  qInfo() << "=== Step 5: energy dialog - find fragment pairs ===";
  structure->ensureDimerMappingTableForCurrentExtent();

  auto selectedFragments = structure->selectedFragments();
  qInfo() << "selected fragments:" << selectedFragments.size();
  if (selectedFragments.size() == 0) {
    qWarning() << "No selected fragments after expansion - aborting repro";
    return 1;
  }
  const auto keyFragmentIndex = selectedFragments[0];
  FragmentPairSettings pairSettings;
  pairSettings.keyFragment = keyFragmentIndex;
  auto fragmentPairs = structure->findFragmentPairs(pairSettings);
  qInfo() << "unique pairs:" << fragmentPairs.uniquePairs.size();

  const auto &fragments = structure->getFragments();
  const auto &keyFragment = fragments.at(keyFragmentIndex);
  FragmentIndexSet wavefunctionsNeeded;
  wavefunctionsNeeded.insert(keyFragment.asymmetricFragmentIndex);
  for (const auto &pair : fragmentPairs.uniquePairs) {
    wavefunctionsNeeded.insert(pair.a.asymmetricFragmentIndex);
    wavefunctionsNeeded.insert(pair.b.asymmetricFragmentIndex);
  }

  // required wavefunctions as built by EnergyCalculationDialog
  std::vector<wfn::Parameters> requiredWavefunctions;
  const auto &uniqueFragments = structure->symmetryUniqueFragments();
  for (const auto &uniqueIndex : wavefunctionsNeeded) {
    const auto &uniqueFrag = uniqueFragments.at(uniqueIndex);
    wfn::Parameters params;
    params.charge = uniqueFrag.state.charge;
    params.multiplicity = uniqueFrag.state.multiplicity;
    // user selected the existing b3lyp/def2-tzvpd wavefunction, so use its
    // method/basis for the request
    params.method = "b3lyp";
    params.basis = "def2-tzvpd";
    params.structure = structure;
    params.atoms = uniqueFrag.atomIndices;
    params.accepted = true;
    requiredWavefunctions.push_back(params);
    qInfo().noquote()
        << QString("required wfn for asym u=%1: %2 atoms")
               .arg(uniqueIndex.u)
               .arg(params.atoms.size());
  }

  qInfo() << "=== Step 6: calculatePairEnergies equivalent check ===";
  for (auto &wfnParams : requiredWavefunctions) {
    auto candidates =
        structure->wavefunctionsAndTransformsForAtoms(wfnParams.atoms);
    qInfo() << "candidates for required wfn (" << wfnParams.atoms.size()
            << "atoms):" << candidates.size();
    bool found = false;
    for (auto &candidate : candidates) {
      qInfo().noquote()
          << QString("  candidate: %1 method=%2 basis=%3 atoms=%4")
                 .arg(candidate.wavefunction->objectName())
                 .arg(candidate.wavefunction->parameters().method)
                 .arg(candidate.wavefunction->parameters().basis)
                 .arg(candidate.wavefunction->atomIndices().size());
      if (wfnParams.hasEquivalentMethodTo(
              candidate.wavefunction->parameters())) {
        found = true;
        break;
      }
    }
    if (!found) {
      qWarning() << "  -> NO equivalent wavefunction: would be recalculated";
    }
  }

  qInfo() << "=== Step 7: pair matching (calculatePairEnergiesWithExistingWavefunctions) ===";
  std::vector<MolecularWavefunction *> wavefunctions;
  for (const auto &wfnParams : requiredWavefunctions) {
    auto candidates =
        structure->wavefunctionsAndTransformsForAtoms(wfnParams.atoms);
    MolecularWavefunction *selectedWfn = nullptr;
    for (auto &candidate : candidates) {
      if (wfnParams.hasEquivalentMethodTo(
              candidate.wavefunction->parameters())) {
        selectedWfn = candidate.wavefunction;
        break;
      }
    }
    if (!selectedWfn && !candidates.empty()) {
      qInfo() << "  fallback to most recent candidate for"
              << wfnParams.atoms.size() << "atoms";
      selectedWfn = candidates.back().wavefunction;
    }
    if (selectedWfn) {
      wavefunctions.push_back(selectedWfn);
    } else {
      qWarning().noquote()
          << QString("  *** required wfn (%1 atoms) has NO candidate at all")
                 .arg(wfnParams.atoms.size());
    }
  }
  qInfo() << "wavefunctions available for pairs:" << wavefunctions.size();

  int computed = 0, skipped = 0, aborted = 0;
  for (const auto &pair : fragmentPairs.uniquePairs) {
    qInfo().noquote()
        << QString("  pair: A=[%1] B=[%2]")
               .arg(composition(structure, pair.a.atomIndices))
               .arg(composition(structure, pair.b.atomIndices));
    pair_energy::Parameters p;
    p.fragmentDimer = pair;
    p.structure = structure;
    p.atomsA = pair.a.atomIndices;
    p.atomsB = pair.b.atomIndices;
    p.model = "ce-1p";

    bool foundA = false, foundB = false;
    for (auto *wf : wavefunctions) {
      if (foundA && foundB)
        break;
      if (!foundA) {
        foundA = structure->getTransformation(wf->atomIndices(), p.atomsA,
                                              p.transformA);
        if (foundA)
          p.wfnA = wf;
      }
      if (!foundB) {
        foundB = structure->getTransformation(wf->atomIndices(), p.atomsB,
                                              p.transformB);
        if (foundB)
          p.wfnB = wf;
      }
    }
    if (!foundA && foundB) {
      aborted++;
      qWarning().noquote()
          << QString("  pair A[%1 atoms] B[%2 atoms]: foundB only -> "
                     "WHOLE CALCULATION ABORTS (return)")
                 .arg(p.atomsA.size())
                 .arg(p.atomsB.size());
      continue;
    }
    if (!foundA || !foundB) {
      skipped++;
      qWarning().noquote()
          << QString("  SKIPPED pair: A composition=[%1] B composition=[%2] "
                     "wfn composition=[%3]")
                 .arg(composition(structure, p.atomsA))
                 .arg(composition(structure, p.atomsB))
                 .arg(composition(structure, wfn->atomIndices()));
      if (skipped == 1) {
        auto dump = [&](const char *name,
                        const std::vector<GenericAtomIndex> &idxs) {
          auto sorted = idxs;
          std::sort(sorted.begin(), sorted.end());
          auto nums = structure->atomicNumbersForIndices(sorted);
          QStringList parts;
          for (int i = 0; i < sorted.size(); i++)
            parts << QString("(%1,%2,%3,%4)->%5")
                         .arg(sorted[i].unique)
                         .arg(sorted[i].x)
                         .arg(sorted[i].y)
                         .arg(sorted[i].z)
                         .arg(nums(i));
          qWarning().noquote()
              << QString("  %1: %2").arg(name).arg(parts.join(" "));
        };
        dump("wfn", wfn->atomIndices());
        dump("A  ", p.atomsA);
        dump("B  ", p.atomsB);
      }
      for (auto *wf : wavefunctions) {
        if (!foundA)
          qWarning().noquote()
              << QString("  pair A[%1 atoms] B[%2 atoms]: no wfn for A "
                         "(wfn %3 has %4 atoms): %5")
                     .arg(p.atomsA.size())
                     .arg(p.atomsB.size())
                     .arg(wf->objectName())
                     .arg(wf->atomIndices().size())
                     .arg(transformFailureReason(structure,
                                                  wf->atomIndices(), p.atomsA));
        if (!foundB)
          qWarning().noquote()
              << QString("  pair A[%1 atoms] B[%2 atoms]: no wfn for B "
                         "(wfn %3 has %4 atoms): %5")
                     .arg(p.atomsA.size())
                     .arg(p.atomsB.size())
                     .arg(wf->objectName())
                     .arg(wf->atomIndices().size())
                     .arg(transformFailureReason(structure,
                                                  wf->atomIndices(), p.atomsB));
      }
      continue;
    }
    computed++;
  }
  qInfo().noquote()
      << QString("RESULT: %1 pairs computed, %2 pairs skipped (missing "
                 "wavefunction), %3 pairs triggered whole-calculation abort")
             .arg(computed)
             .arg(skipped)
             .arg(aborted);
  return 0;
}
