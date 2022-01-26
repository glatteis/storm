#ifndef STORM_EQUALPARAMATERREDUCER_H
#define STORM_EQUALPARAMATERREDUCER_H

#include "adapters/RationalFunctionAdapter.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "storm-pars/utility/parametric.h"

namespace storm {
namespace derivative {

class EqualParameterReducer {
    
   public:
        EqualParameterReducer() {}
        models::sparse::Dtmc<RationalFunction> minimizeEqualParameters(models::sparse::Dtmc<RationalFunction> model, modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask);
};

}  // namespace derivative
}  // namespace storm
    
#endif