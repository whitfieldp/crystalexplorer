#include "load_pair_energy_json.h"
#include "json.h"
#include <QFile>

PairInteraction *load_pair_energy_json(const QString &filename) {
  constexpr double hartree_to_kj_per_mol{2625.5};
  PairInteraction *result{nullptr};

  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly)) {
    qWarning() << "Couldn't open pair energy json file:" << filename;
    return result;
  }
  QByteArray data = file.readAll();

  nlohmann::json doc;
  try {
    doc = nlohmann::json::parse(data.constData());
  } catch (const nlohmann::json::exception &e) {
    qWarning() << "JSON parse error in pair energy file:" << filename << e.what();
    return result;
  }

  if (!doc.contains("interaction_model") || !doc.contains("interaction_energy") ||
      !doc["interaction_energy"].is_object()) {
    qWarning() << "Pair energy JSON missing/invalid required fields:" << filename
               << "type:" << doc.type_name();
    return result;
  }

  QString interactionModel;
  try {
    interactionModel =
        doc.at("interaction_model").at("name").get<QString>();
  } catch (const nlohmann::json::exception &e) {
    qWarning() << "Could not read interaction_model.name:" << e.what();
    interactionModel = "unknown";
  }

  result = new PairInteraction(interactionModel);
  for (const auto &item : doc["interaction_energy"].items()) {
    // Only accept numeric components; skip any non-numeric entry (e.g. a
    // nested breakdown or an error placeholder) rather than throwing.
    if (!item.value().is_number())
      continue;
    double value = item.value().get<double>() * hartree_to_kj_per_mol;
    result->addComponent(QString::fromStdString(item.key()), value);
  }

  return result;
}
