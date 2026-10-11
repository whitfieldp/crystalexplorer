#pragma once
#include "wavefunction_parameters.h"

namespace io {

QString getOccWavefunctionJson(const wfn::Parameters &params);

// Read the model section of (possibly user edited) OCC wavefunction input
// and update params.method / params.basis if present.
// Returns false if the JSON could not be parsed or has no usable model
// section, in which case params is left unchanged.
bool updateWfnParametersFromOccJson(const QString &json,
                                    wfn::Parameters &params);

}
