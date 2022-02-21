#include "DerivativeBoundFinder.h"
#include <cstdint>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include "adapters/RationalFunctionAdapter.h"
#include "api/bisimulation.h"
#include "logic/FormulaContext.h"
#include "logic/FormulasForwardDeclarations.h"
#include "logic/ProbabilityOperatorFormula.h"
#include "logic/UntilFormula.h"
#include "modelchecker/prctl/SparseDtmcPrctlModelChecker.h"
#include "modelchecker/propositional/SparsePropositionalModelChecker.h"
#include "modelchecker/results/CheckResult.h"
#include "modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "models/sparse/Dtmc.h"
#include "models/sparse/ItemLabeling.h"
#include "models/sparse/Model.h"
#include "models/sparse/StateLabeling.h"
#include "settings/modules/GeneralSettings.h"
#include "storage/BitVector.h"
#include "storage/SparseMatrix.h"
#include "storage/bisimulation/BisimulationType.h"
#include "storm-pars/analysis/MonotonicityResult.h"
#include "storm-pars/modelchecker/instantiation/SparseDtmcInstantiationModelChecker.h"
#include "storm-pars/transformer/SparseParametricDtmcSimplifier.h"
#include "storm-pars/utility/parametric.h"
#include "storm-parsers/parser/FormulaParser.h"
#include "utility/constants.h"
#include "utility/graph.h"
#include "utility/logging.h"
#include "utility/macros.h"
#include "modelchecker/results/ExplicitQuantitativeCheckResult.h"

namespace storm {
namespace derivative {

template<typename FunctionType>
using VariableType = typename utility::parametric::VariableType<FunctionType>::type;
template<typename FunctionType>
using CoefficientType = typename utility::parametric::CoefficientType<FunctionType>::type;

// template<typename FunctionType, typename ConstantType>
// std::pair<std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>,
//           std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>>
// DerivativeBoundFinder<FunctionType, ConstantType>::getDerivativeBound(Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region,
//                                                                       VariableType<FunctionType> parameter) {
//     auto const sharedModelPointer = std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(model);
//     this->liftingModelChecker->specify_internal(env, sharedModelPointer, *this->currentCheckTaskNoBound, false, true);

//     auto modelCopyInitial = model;
//     std::vector<ConstantType> max = liftingModelChecker->check(env, region, OptimizationDirection::Maximize, nullptr)
//                                         ->template asExplicitQuantitativeCheckResult<ConstantType>()
//                                         .getValueVector();
//     /* std::cout << storm::utility::vector::toString(liftingModelChecker->check(env, region, OptimizationDirection::Minimize, nullptr) */
//     /*                                     ->template asExplicitQuantitativeCheckResult<ConstantType>() */
//     /*                                     .getValueVector()) << std::endl; */
//     std::vector<ConstantType> min = liftingModelChecker->check(env, region, OptimizationDirection::Minimize, nullptr)
//                                         ->template asExplicitQuantitativeCheckResult<ConstantType>()
//                                         .getValueVector();

//     const storage::SparseMatrix<FunctionType> transitionMatrix = model.getTransitionMatrix();
//     models::sparse::Dtmc<FunctionType> modelCopy = model;

//     std::vector<FunctionType> stateRewardsMax(transitionMatrix.getRowCount());
//     std::vector<FunctionType> stateRewardsMin(transitionMatrix.getRowCount());

//     /* for (uint_fast64_t i = 0; i < modelCopy.getNumberOfStates(); i++) { */
//     /*     std::cout << "max" << i << ":" << max[i] << std::endl; */
//     /*     std::cout << "min" << i << ":" << min[i] << std::endl; */
//     /* } */

//     for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
//         if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
//             if (currentCheckTaskNoBound->isRewardModelSet()) {
//                 /* std::cout << "reward for " << i << ":" << model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i] <<
//                  * std::endl; */
//                 stateRewardsMax[i] = model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i].derivative(parameter);
//                 stateRewardsMin[i] = model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i].derivative(parameter);
//             } else {
//                 /* std::cout << "reward for " << i << ":" << model.getRewardModel("").getStateRewardVector()[i] << std::endl; */
//                 stateRewardsMax[i] = model.getRewardModel("").getStateRewardVector()[i].derivative(parameter);
//                 stateRewardsMin[i] = model.getRewardModel("").getStateRewardVector()[i].derivative(parameter);
//             }
//         } else {
//             stateRewardsMax[i] = utility::zero<FunctionType>();
//             stateRewardsMin[i] = utility::zero<FunctionType>();
//         }

//         for (auto const& entry : transitionMatrix.getRow(i)) {
//             if (!entry.getValue().derivative(parameter).isConstant()) {
//                 STORM_LOG_ERROR("Non-simple models not supported. This model has a transition with derivative " << entry.getValue().derivative(parameter));
//             }
//             ConstantType derivative = utility::convertNumber<ConstantType>(entry.getValue().derivative(parameter));
//             uint_fast64_t toState = entry.getColumn();
//             if (derivative != utility::zero<ConstantType>()) {
//                 ConstantType extremalValueMin = 0;
//                 ConstantType extremalValueMax = 0;
//                 if (derivative < utility::zero<ConstantType>()) {
//                     extremalValueMax = min[toState];
//                     extremalValueMin = max[toState];
//                 } else if (derivative > utility::zero<ConstantType>()) {
//                     extremalValueMax = max[toState];
//                     extremalValueMin = min[toState];
//                 }

//                 stateRewardsMax[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivative * extremalValueMax)));
//                 stateRewardsMin[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivative * extremalValueMin)));
//             }
//         }
//     }

//     /* for (uint_fast64_t i = 0; i < modelCopy.getNumberOfStates(); i++) { */
//     /*     std::cout << "for " << i << ": " << utility::convertNumber<double>(stateRewardsMin[i]) << "<=" << utility::convertNumber<double>(stateRewardsMax[i])
//      * << std::endl; */
//     /* } */

//     models::sparse::StandardRewardModel<FunctionType> rewardModelMax(std::move(stateRewardsMax));
//     models::sparse::StandardRewardModel<FunctionType> rewardModelMin(std::move(stateRewardsMin));

//     modelCopy.addRewardModel("derivative-max", rewardModelMax);
//     modelCopy.addRewardModel("derivative-min", rewardModelMin);

//     storage::BitVector target = modelCopy.getStates("target");
//     storm::storage::BitVector newTarget(target.size());

//     if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
//         newTarget = target;
//     } else {
//         storm::storage::BitVector probZero = storm::utility::graph::performProbGreater0(modelCopy.getBackwardTransitions(),
//                                                                                         storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);
//         probZero.complement();
//         storm::storage::BitVector probOne =
//             storm::utility::graph::performProb1(modelCopy.getBackwardTransitions(), storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);

//         newTarget |= probZero;
//         newTarget |= probOne;
//     }

//     modelCopy.getStateLabeling().addLabel("derivative-target");
//     modelCopy.getStateLabeling().setStates("derivative-target", newTarget);

//     /* modelCopy.writeDotToStream(std::cout); */
//     /* for (uint_fast64_t i = 0; i < modelCopy.getNumberOfStates(); i++) { */
//     /*     std::cout << "R" << i << ":" << rewardModelMax.getStateReward(i) << std::endl; */
//     /* } */

//     auto subformulaConstructor = std::make_shared<logic::AtomicLabelFormula>("derivative-target");
//     auto subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none);

//     /* std::shared_ptr<const storm::logic::Formula> formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula,
//      * std::string("derivative-max"), this->formulaOperatorInformation,
//      * logic::RewardMeasureType::Expectation)->asRewardOperatorFormula().asSharedPointer(); */
//     auto formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-max"), this->formulaOperatorInformation);
//     auto formulaMin = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-min"), this->formulaOperatorInformation);
//     auto checkTaskMin = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMin);
//     auto checkTaskMax = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);

//     this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMax);
//     auto derivativeMax = liftingModelChecker->getBound(env, region, OptimizationDirection::Maximize)
//                              ->template asExplicitQuantitativeCheckResult<ConstantType>()
//                              .getValueVector();

//     this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMin);
//     auto derivativeMin = liftingModelChecker->getBound(env, region, OptimizationDirection::Minimize)
//                              ->template asExplicitQuantitativeCheckResult<ConstantType>()
//                              .getValueVector();

//     uint_fast64_t initialState;
//     const storm::storage::BitVector initialVector = model.getInitialStates();
//     for (uint_fast64_t x : initialVector) {
//         initialState = x;
//         break;
//     }

//     /* std::cout << "max at init: " << derivativeMax[initialState] << std::endl; */
//     /* std::cout << "min at init: " << derivativeMin[initialState] << std::endl; */

//     auto resultMax = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMax);
//     auto resultMin = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMin);

//     return std::make_pair(std::move(resultMax), std::move(resultMin));
// }


template<typename FunctionType, typename ConstantType>
storm::models::sparse::Dtmc<FunctionType> const DerivativeBoundFinder<FunctionType, ConstantType>::getInternalModel() {
    return model;
}

template<typename FunctionType, typename ConstantType>
std::pair<std::pair<models::sparse::Dtmc<FunctionType>, models::sparse::Dtmc<FunctionType>>,
          std::pair<std::shared_ptr<storm::logic::Formula>, std::shared_ptr<storm::logic::Formula>>>
DerivativeBoundFinder<FunctionType, ConstantType>::computeMonotonicityTasks(
    Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region, std::vector<ConstantType> minValues,
    std::vector<ConstantType> maxValues, std::shared_ptr<storm::analysis::LocalMonotonicityResult<VariableType<FunctionType>>> localMonotonicityResult,
    VariableType<FunctionType> const& parameter) {
    auto const sharedModelPointer = std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(model);
    // this->liftingModelChecker->specify_internal(env, sharedModelPointer, *this->currentCheckTaskNoBound, false, true);

    auto modelCopyInitial = model;

    const storage::SparseMatrix<FunctionType> transitionMatrix = model.getTransitionMatrix();
    models::sparse::Dtmc<FunctionType> modelCopy = model;
    
    std::vector<FunctionType> stateRewardsMax(transitionMatrix.getRowCount());
    std::vector<FunctionType> stateRewardsMin(transitionMatrix.getRowCount());
    
    for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
        if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
            if (currentCheckTaskNoBound->isRewardModelSet()) {
                /* std::cout << "reward for " << i << ":" << model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i] <<
                 * std::endl; */
                stateRewardsMax[i] = model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i].derivative(parameter);
                stateRewardsMin[i] = model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector()[i].derivative(parameter);
            } else {
                /* std::cout << "reward for " << i << ":" << model.getRewardModel("").getStateRewardVector()[i] << std::endl; */
                stateRewardsMax[i] = model.getRewardModel("").getStateRewardVector()[i].derivative(parameter);
                stateRewardsMin[i] = model.getRewardModel("").getStateRewardVector()[i].derivative(parameter);
            }
        } else {
            stateRewardsMax[i] = utility::zero<FunctionType>();
            stateRewardsMin[i] = utility::zero<FunctionType>();
        }
        
        ConstantType derivativeMin = 0.0;
        ConstantType derivativeMax = 0.0;

        for (auto const& entry : transitionMatrix.getRow(i)) {
            if (!entry.getValue().derivative(parameter).isConstant()) {
                STORM_LOG_ERROR("Non-simple models not supported. This model has a transition with derivative " << entry.getValue().derivative(parameter));
            }
            ConstantType derivative = utility::convertNumber<ConstantType>(entry.getValue().derivative(parameter));
            uint_fast64_t toState = entry.getColumn();
            if (derivative != utility::zero<ConstantType>()) {
                ConstantType extremalValueMin = 0;
                ConstantType extremalValueMax = 0;
                if (derivative < utility::zero<ConstantType>()) {
                    extremalValueMax = minValues[toState];
                    extremalValueMin = maxValues[toState];
                } else if (derivative > utility::zero<ConstantType>()) {
                    extremalValueMax = maxValues[toState];
                    extremalValueMin = minValues[toState];
                }
                
                derivativeMax += derivative * extremalValueMax;
                derivativeMin += derivative * extremalValueMin;
            }
        }
        
        if (localMonotonicityResult != nullptr) {
            if (derivativeMax < 0) {
                localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Decr);
            } else if (derivativeMin > 0) {
                localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Incr);
            } else  if (derivativeMin == 0 && derivativeMax == 0) {
                localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Constant);
            } else {
                localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Unknown);
            }
        }

        stateRewardsMax[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivativeMax)));
        stateRewardsMin[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivativeMin)));
    }

    models::sparse::StandardRewardModel<FunctionType> rewardModelMax(std::move(stateRewardsMax));
    models::sparse::StandardRewardModel<FunctionType> rewardModelMin(std::move(stateRewardsMin));

    modelCopy.addRewardModel("derivative-max", rewardModelMax);
    modelCopy.addRewardModel("derivative-min", rewardModelMin);

    storm::modelchecker::SparsePropositionalModelChecker<models::sparse::Dtmc<FunctionType>> propositionalChecker(modelCopy);
    storage::BitVector phiStates;
    storage::BitVector psiStates;
    if (this->currentFormula->isProbabilityOperatorFormula()) {
        if (this->currentFormula->asProbabilityOperatorFormula().getSubformula().isUntilFormula()) {
            phiStates = propositionalChecker.check(this->currentFormula->asProbabilityOperatorFormula().getSubformula().asUntilFormula().getLeftSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
            psiStates = propositionalChecker.check(this->currentFormula->asProbabilityOperatorFormula().getSubformula().asUntilFormula().getRightSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
        } else {
            STORM_LOG_ASSERT(this->currentFormula->asProbabilityOperatorFormula().getSubformula().isEventuallyFormula(),
                             "Expecting formula to be until or eventually formula");
            phiStates = storage::BitVector(modelCopy.getNumberOfStates(), true);
            psiStates = propositionalChecker.check(this->currentFormula->asProbabilityOperatorFormula().getSubformula().asEventuallyFormula().getSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
        }
    } else if (this->currentFormula->isRewardOperatorFormula()) {
        if (this->currentFormula->asRewardOperatorFormula().getSubformula().isUntilFormula()) {
            phiStates = propositionalChecker.check(this->currentFormula->asRewardOperatorFormula().getSubformula().asUntilFormula().getLeftSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
            psiStates = propositionalChecker.check(this->currentFormula->asRewardOperatorFormula().getSubformula().asUntilFormula().getRightSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
        } else {
            STORM_LOG_ASSERT(this->currentFormula->asRewardOperatorFormula().getSubformula().isEventuallyFormula(),
                             "Expecting formula to be until or eventually formula");
            phiStates = storage::BitVector(modelCopy.getNumberOfStates(), true);
            psiStates = propositionalChecker.check(this->currentFormula->asRewardOperatorFormula().getSubformula().asEventuallyFormula().getSubformula())
                            ->asExplicitQualitativeCheckResult()
                            .getTruthValuesVector();
        }
    }
    std::pair<storage::BitVector, storage::BitVector> statesWithProbability01 =
        utility::graph::performProb01(modelCopy.getBackwardTransitions(), phiStates, psiStates);
    storage::BitVector target = statesWithProbability01.second;
    storage::BitVector bottomStates = statesWithProbability01.first;

    storm::storage::BitVector newPsiStates(target.size());

    if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
        newPsiStates = psiStates;
    } else {
        storm::storage::BitVector probZero = storm::utility::graph::performProbGreater0(modelCopy.getBackwardTransitions(),
                                                                                        storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);
        probZero.complement();
        storm::storage::BitVector probOne =
            storm::utility::graph::performProb1(modelCopy.getBackwardTransitions(), storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);

        // We are turning a probability operator formula into a reward operator formula.
        // The actual probability of reaching the target should now be one to collect all rewards.
        newPsiStates |= probZero;
        newPsiStates |= probOne;
    }

    modelCopy.getStateLabeling().addLabel("derivative-safe");
    modelCopy.getStateLabeling().setStates("derivative-safe", phiStates);
    modelCopy.getStateLabeling().addLabel("derivative-target");
    modelCopy.getStateLabeling().setStates("derivative-target", newPsiStates);

    /* modelCopy.writeDotToStream(std::cout); */
    /* for (uint_fast64_t i = 0; i < modelCopy.getNumberOfStates(); i++) { */
    /*     std::cout << "R" << i << ":" << rewardModelMax.getStateReward(i) << std::endl; */
    /* } */
    /* std::shared_ptr<const storm::logic::Formula> formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula,
     * std::string("derivative-max"), this->formulaOperatorInformation,
     * logic::RewardMeasureType::Expectation)->asRewardOperatorFormula().asSharedPointer(); */


    std::shared_ptr<storm::logic::Formula> subformula;
    if ((this->currentFormula->isRewardOperatorFormula() && this->currentFormula->asRewardOperatorFormula().getSubformula().isEventuallyFormula()) ||
        (this->currentFormula->isProbabilityOperatorFormula() && this->currentFormula->asProbabilityOperatorFormula().getSubformula().isEventuallyFormula())) {
        auto subformulaConstructor = std::make_shared<logic::AtomicLabelFormula>("derivative-target");
        subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none);
    } else {
        auto subformulaConstructor1 = std::make_shared<logic::AtomicLabelFormula>("derivative-safe");
        auto subformulaConstructor2 = std::make_shared<logic::AtomicLabelFormula>("derivative-target");
        subformula = std::make_shared<logic::UntilFormula>(subformulaConstructor1, subformulaConstructor2);
    }

    auto formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-max"), this->formulaOperatorInformation);
    auto formulaMin = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-min"), this->formulaOperatorInformation);
    
    // std::cout << *formulaMax << std::endl;
    // std::cout << *formulaMin << std::endl;

    storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<FunctionType>> simplifier(modelCopy);

    // STORM_LOG_ASSERT(simplifier.simplify(*formulaMax), "Could not simplify derivative model.");
    // auto modelMax = simplifier.getSimplifiedModel();

    // STORM_LOG_ASSERT(simplifier.simplify(*formulaMin), "Could not simplify derivative model.");
    // auto modelMin = simplifier.getSimplifiedModel();
    
    // modelMax->reduceToStateBasedRewards();
    // modelMin->reduceToStateBasedRewards();
    
    auto modelMax = model;
    auto modelMin = model;

    // std::cout << parameter << std::endl;
    // std::cout << "Model copy:" << std::endl;
    // modelCopy.writeDotToStream(std::cout, 30, true, nullptr, &modelCopy.getRewardModel("derivative-min").getStateRewardVector(),
    //                            &modelCopy.getRewardModel("derivative-max").getStateRewardVector());
    // std::cout << "Model max:" << std::endl;
    // modelMax.writeDotToStream(std::cout, 30, true, nullptr, &modelMax.getUniqueRewardModel().getStateRewardVector());
    // std::cout << "Model min:" << std::endl;
    // modelMin.writeDotToStream(std::cout, 30, true, nullptr, &modelMin.getUniqueRewardModel().getStateRewardVector());

    // modelMax = storm::api::performDeterministicSparseBisimulationMinimization(modelMax, {formulaMax}, storage::BisimulationType::Strong);
    // modelMin = storm::api::performDeterministicSparseBisimulationMinimization(modelMin, {formulaMin}, storage::BisimulationType::Strong);
    
    return std::make_pair(std::make_pair(modelCopy, modelCopy), std::make_pair(formulaMin, formulaMax));
}

template<typename FunctionType, typename ConstantType>
void DerivativeBoundFinder<FunctionType, ConstantType>::updateMonotonicityResult(
    std::vector<ConstantType> derivativeMinValues, std::vector<ConstantType> derivativeMaxValues,
    std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> localMonotonicityResult,
    VariableType<FunctionType> const& parameter, uint_fast64_t initialState) {
    // std::cout << derivativeMinValues[initialState] << " <= d" << parameter << " <= " << derivativeMaxValues[initialState] << std::endl;
    boost::optional<typename analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity> finalResult;
    if (derivativeMaxValues[initialState] < 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Decr;
    } else if (derivativeMinValues[initialState] > 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Incr;
    } else if (derivativeMaxValues[initialState] == 0 && derivativeMinValues[initialState] == 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Constant;
    } else {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Unknown;
    }

    if (finalResult) {
        localMonotonicityResult->getGlobalMonotonicityResult()->addMonotonicityResult(parameter, *finalResult);
        for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
            localMonotonicityResult->setMonotonicity(i, parameter, *finalResult);
        }
    }
}

template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
template class DerivativeBoundFinder<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm
