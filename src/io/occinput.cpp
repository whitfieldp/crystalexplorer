#include "occinput.h"
#include "json.h"
#include "chemicalstructure.h"

namespace io {

QString getOccWavefunctionJson(const wfn::Parameters &params) {
  nlohmann::json root;
  nlohmann::json topology;

  auto nums = params.structure->atomicNumbersForIndices(params.atoms);
  auto pos = params.structure->atomicPositionsForIndices(params.atoms);

  std::vector<std::string> symbols;
  std::vector<double> positions;

  // should be in angstroms
  for (int i = 0; i < nums.rows(); i++) {
    symbols.push_back(occ::core::Element(nums(i)).symbol());
    positions.push_back(pos(0, i));
    positions.push_back(pos(1, i));
    positions.push_back(pos(2, i));
  }
  root["schema_name"] = "qcschema_input";
  root["scema_version"] = 1;
  root["return_output"] = true;

  root["molecule"]["geometry"] = positions;
  root["molecule"]["symbols"] = symbols;
  root["molecule"]["molecular_charge"] = params.charge;
  root["molecule"]["molecular_multiplicity"] = params.multiplicity;
  root["driver"] = "energy";

  root["model"]["method"] = params.method;
  root["model"]["basis"] = params.basis;

  return QString::fromStdString(root.dump(2));
}

bool updateWfnParametersFromOccJson(const QString &json,
                                    wfn::Parameters &params) {
  nlohmann::json doc =
      nlohmann::json::parse(json.toStdString(), nullptr, false);
  if (doc.is_discarded() || !doc.contains("model") ||
      !doc["model"].is_object()) {
    qWarning()
        << "Edited OCC input has no parseable model section; keeping "
           "original method and basis set";
    return false;
  }
  const auto &model = doc["model"];
  if (model.contains("method") && model["method"].is_string()) {
    params.method = QString::fromStdString(model["method"].get<std::string>());
  }
  if (model.contains("basis") && model["basis"].is_string()) {
    params.basis = QString::fromStdString(model["basis"].get<std::string>());
  }
  return true;
}

} // namespace io
