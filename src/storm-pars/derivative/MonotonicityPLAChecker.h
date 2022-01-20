#ifndef STORM_MONOTONICITYPLACHECKER_H
#define STORM_MONOTONICITYPLACHECKER_H

#include <cstdint>
#include <memory>
#include "DerivativeBoundFinder.h"
#include "environment/Environment.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "storm-pars/modelchecker/region/SparseDtmcParameterLiftingModelChecker.h"
namespace storm {
namespace derivative {

template<typename FunctionType, typename ConstantType>
class MonotonicityPLAChecker {
   public:
    MonotonicityPLAChecker<FunctionType, ConstantType>(storm::models::sparse::Dtmc<FunctionType> const model) : model(model), boundFinder(model) {}

    void specifyFormula(Environment const& env, modelchecker::CheckTask<logic::Formula, FunctionType> const& checkTask) {
        boundFinder.specifyFormula(env, checkTask);
        auto currentFormula = checkTask.getFormula().asSharedPointer();
        currentCheckTask = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(checkTask.substituteFormula(*currentFormula));
    }

    void performMonotonicityPLA(Environment const& env, typename utility::parametric::VariableType<FunctionType>::type wrt, ConstantType terminateArea);

    std::pair<std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>,
              std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>>
    getDerivativeBound(Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region,
                       typename utility::parametric::VariableType<FunctionType>::type parameter);
    uint_fast64_t getInitialState();

    const models::sparse::Dtmc<FunctionType> model;
    derivative::DerivativeBoundFinder<FunctionType, ConstantType> boundFinder;
    modelchecker::SparseDtmcParameterLiftingModelChecker<models::sparse::Dtmc<FunctionType>, ConstantType> modelChecker;
    std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, FunctionType>> currentCheckTask;
};
}  // namespace derivative
}  // namespace storm

#endif