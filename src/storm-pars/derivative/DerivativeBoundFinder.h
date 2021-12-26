#ifndef STORM_DERIVATIVEBOUNDFINDER_H
#define STORM_DERIVATIVEBOUNDFINDER_H

#include <memory>
#include <queue>
#include "builder/TerminalStatesGetter.h"
#include "environment/Environment.h"
#include "logic/Formula.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "storm-pars/analysis/MonotonicityResult.h"
#include "storm-pars/transformer/ParameterLifter.h"
#include "storm-pars/utility/parametric.h"
#include "utility/macros.h"
namespace storm {
namespace derivative {

template<typename FunctionType, typename ConstantType>
class DerivativeBoundFinder {
   public:
    DerivativeBoundFinder<FunctionType, ConstantType>(storm::models::sparse::Dtmc<FunctionType> const model) : model(model) {}

    void specifyFormula(Environment const& env, modelchecker::CheckTask<logic::Formula, FunctionType> const& checkTask) {
        this->currentFormula = checkTask.getFormula().asSharedPointer();
        this->currentCheckTask = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(
            checkTask.substituteFormula(*currentFormula).template convertValueType<FunctionType>());

        if (!checkTask.getFormula().isRewardOperatorFormula()) {
            this->currentFormulaNoBound = std::make_shared<storm::logic::ProbabilityOperatorFormula>(
                checkTask.getFormula().asProbabilityOperatorFormula().getSubformula().asSharedPointer(),
                storm::logic::OperatorInformation(boost::none, boost::none));
            /* auto subformulaConstructor =
             * checkTask.getFormula().asProbabilityOperatorFormula().getSubformula().asEventuallyFormula().getSubformula().asSharedPointer(); */
            /* this->subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none); */
            this->formulaOperatorInformation = checkTask.getFormula().asProbabilityOperatorFormula().getOperatorInformation();
        } else {
            // No worries, this works as intended, the API is just weird.
            this->currentFormulaNoBound =
                std::make_shared<storm::logic::RewardOperatorFormula>(checkTask.getFormula().asRewardOperatorFormula().getSubformula().asSharedPointer());
            /* auto subformulaConstructor =
             * checkTask.getFormula().asRewardOperatorFormula().getSubformula().asEventuallyFormula().getSubformula().asSharedPointer(); */
            /* this->subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none); */
            this->formulaOperatorInformation = checkTask.getFormula().asRewardOperatorFormula().getOperatorInformation();
        }
        auto terminalStates = storm::builder::getTerminalStatesFromFormula(*this->currentFormulaNoBound);
        STORM_LOG_ASSERT(terminalStates.terminalExpressions.size() == 1, "Model needs one terminal label!");
        this->terminalExpression = terminalStates.terminalExpressions[0];
        
        this->currentCheckTaskNoBound = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*currentFormulaNoBound);
        this->currentCheckTaskNoBoundConstantType =
            std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, ConstantType>>(*currentFormulaNoBound);

        this->parameters = storm::models::sparse::getProbabilityParameters(model);
        if (checkTask.getFormula().isRewardOperatorFormula()) {
            for (auto const& rewardParameter : storm::models::sparse::getRewardParameters(model)) {
                this->parameters.insert(rewardParameter);
            }
        }
    }

    std::pair<models::sparse::Dtmc<FunctionType>, std::pair<std::shared_ptr<storm::logic::Formula>,
                                                            std::shared_ptr<storm::logic::Formula>>>
    computeMonotonicityTasks(
        Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region, std::vector<ConstantType> minValues,
        std::vector<ConstantType> maxValues,
        std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> localMonotonicityResult,
        // std::shared_ptr<storm::analysis::MonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> globalMonotonicityResult,
        typename utility::parametric::VariableType<FunctionType>::type const& parameter
        );

    void
    updateMonotonicityResult(
        std::vector<ConstantType> derivativeMinValues, std::vector<ConstantType> derivativeMaxValues,
        std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> localMonotonicityResult,
        // std::shared_ptr<storm::analysis::MonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> globalMonotonicityResult,
        typename utility::parametric::VariableType<FunctionType>::type const& parameter,
        uint_fast64_t initialState
    );

    // void derivativePLASketch(Environment const& env, typename utility::parametric::VariableType<FunctionType>::type parameter, ConstantType terminateArea);

    // std::pair<std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>,
    //           std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>>
    // getDerivativeBound(Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region,
    //                    typename utility::parametric::VariableType<FunctionType>::type parameter);

   private:
    std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, FunctionType>> currentCheckTask;
    std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, FunctionType>> currentCheckTaskNoBound;
    std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, ConstantType>> currentCheckTaskNoBoundConstantType;
    std::shared_ptr<storm::logic::Formula const> currentFormula;
    std::shared_ptr<storm::logic::Formula const> currentFormulaNoBound;
    std::unique_ptr<storm::logic::Formula const> terminalExpression;

    /* std::shared_ptr<storm::logic::EventuallyFormula const> subformula; */
    logic::OperatorInformation formulaOperatorInformation;

    const models::sparse::Dtmc<FunctionType> model;
    std::set<typename utility::parametric::VariableType<FunctionType>::type> parameters;
};

}  // namespace derivative
}  // namespace storm

#endif
