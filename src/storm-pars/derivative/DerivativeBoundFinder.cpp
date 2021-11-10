#include "DerivativeBoundFinder.h"
#include "storm-pars/modelchecker/instantiation/SparseDtmcInstantiationModelChecker.h"
#include "storm-pars/utility/parametric.h"
#include "storm-parsers/parser/FormulaParser.h"

namespace storm {
namespace derivative {

template<typename FunctionType>
using VariableType = typename utility::parametric::VariableType<FunctionType>::type;
template<typename FunctionType>
using CoefficientType = typename utility::parametric::CoefficientType<FunctionType>::type;

template<typename FunctionType, typename ConstantType>
std::pair<std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>,
          std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>>
DerivativeBoundFinder<FunctionType, ConstantType>::getDerivativeBound(Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region,
                                                                      VariableType<FunctionType> parameter) {
    this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(model), *this->currentCheckTaskNoBound);
    std::vector<ConstantType> min = liftingModelChecker->getBound(env, region, OptimizationDirection::Minimize, nullptr)
                                        ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                        .getValueVector();
    std::vector<ConstantType> max = liftingModelChecker->getBound(env, region, OptimizationDirection::Maximize, nullptr)
                                        ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                        .getValueVector();

    const storage::SparseMatrix<FunctionType> transitionMatrix = model.getTransitionMatrix();
    models::sparse::Dtmc<FunctionType> modelCopy = model;

    std::vector<FunctionType> stateRewardsMax(transitionMatrix.getRowCount());
    std::vector<FunctionType> stateRewardsMin(transitionMatrix.getRowCount());

    for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
        stateRewardsMax[i] = model.getUniqueRewardModel().getStateRewardVector()[i];
        stateRewardsMin[i] = model.getUniqueRewardModel().getStateRewardVector()[i];

        for (auto const& entry : transitionMatrix.getRow(i)) {
            ConstantType derivative = utility::convertNumber<ConstantType>(entry.getValue().derivative(parameter));
            uint_fast64_t toState = entry.getColumn();
            /* std::cout << entry << " , " << parameter << " , " << derivative << std::endl; */
            if (derivative != utility::zero<ConstantType>()) {
                ConstantType extremalValueMin = 0;
                ConstantType extremalValueMax = 0;
                if (derivative < utility::zero<ConstantType>()) {
                    extremalValueMax = min[toState];
                    extremalValueMin = max[toState];
                } else if (derivative > utility::zero<ConstantType>()) {
                    extremalValueMax = max[toState];
                    extremalValueMin = min[toState];
                }

                stateRewardsMax[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivative * extremalValueMax)));
                stateRewardsMin[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivative * extremalValueMin)));
            }
        }
    }

    models::sparse::StandardRewardModel<FunctionType> rewardModelMax(std::move(stateRewardsMax));
    models::sparse::StandardRewardModel<FunctionType> rewardModelMin(std::move(stateRewardsMin));

    modelCopy.addRewardModel("derivative-max", rewardModelMax);
    modelCopy.addRewardModel("derivative-min", rewardModelMin);

    storage::BitVector target = modelCopy.getStates("target");
    storm::storage::BitVector probZero =
        storm::utility::graph::performProbGreater0(modelCopy.getBackwardTransitions(), storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);
    probZero.complement();
    storm::storage::BitVector probOne =
        storm::utility::graph::performProb1(modelCopy.getBackwardTransitions(), storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);

    storm::storage::BitVector newTarget(probZero.size());
    newTarget |= probZero;
    newTarget |= probOne;

    modelCopy.getStateLabeling().addLabel("derivative-target");
    modelCopy.getStateLabeling().setStates("derivative-target", newTarget);

    auto subformulaConstructor = std::make_shared<logic::AtomicLabelFormula>("derivative-target");
    auto subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none);

    /* std::shared_ptr<const storm::logic::Formula> formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula,
     * std::string("derivative-max"), this->formulaOperatorInformation,
     * logic::RewardMeasureType::Expectation)->asRewardOperatorFormula().asSharedPointer(); */
    auto formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-max"), this->formulaOperatorInformation);
    auto formulaMin = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-min"), this->formulaOperatorInformation);
    auto checkTaskMin = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);
    auto checkTaskMax = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);

    this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMax);
    auto derivativeMax = liftingModelChecker->getBound(env, region, OptimizationDirection::Maximize)
                             ->template asExplicitQuantitativeCheckResult<ConstantType>()
                             .getValueVector();

    this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMin);
    auto derivativeMin = liftingModelChecker->getBound(env, region, OptimizationDirection::Minimize)
                             ->template asExplicitQuantitativeCheckResult<ConstantType>()
                             .getValueVector();

    uint_fast64_t initialState;
    const storm::storage::BitVector initialVector = model.getInitialStates();
    for (uint_fast64_t x : initialVector) {
        initialState = x;
        break;
    }

    std::cout << "max at init: " << derivativeMax[initialState] << std::endl;
    std::cout << "min at init: " << derivativeMin[initialState] << std::endl;

    auto resultMax = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMax);
    auto resultMin = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMin);
}

template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
template class DerivativeBoundFinder<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm
