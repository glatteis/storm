#include "DerivativeBoundFinder.h"
#include <_types/_uint64_t.h>
#include <cstdint>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include "adapters/RationalFunctionAdapter.h"
#include "api/bisimulation.h"
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

    storm::storage::BitVector newTarget(target.size());

    if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
        newTarget = target;
    } else {
        storm::storage::BitVector probZero = storm::utility::graph::performProbGreater0(modelCopy.getBackwardTransitions(),
                                                                                        storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);
        probZero.complement();
        storm::storage::BitVector probOne =
            storm::utility::graph::performProb1(modelCopy.getBackwardTransitions(), storm::storage::BitVector(modelCopy.getNumberOfStates(), true), target);

        newTarget |= probZero;
        newTarget |= probOne;
    }

    modelCopy.getStateLabeling().addLabel("derivative-target");
    modelCopy.getStateLabeling().setStates("derivative-target", newTarget);

    /* modelCopy.writeDotToStream(std::cout); */
    /* for (uint_fast64_t i = 0; i < modelCopy.getNumberOfStates(); i++) { */
    /*     std::cout << "R" << i << ":" << rewardModelMax.getStateReward(i) << std::endl; */
    /* } */
    /* std::shared_ptr<const storm::logic::Formula> formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula,
     * std::string("derivative-max"), this->formulaOperatorInformation,
     * logic::RewardMeasureType::Expectation)->asRewardOperatorFormula().asSharedPointer(); */

    auto subformulaConstructor = std::make_shared<logic::AtomicLabelFormula>("derivative-target");
    auto subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Reward, boost::none);

    auto formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-max"), this->formulaOperatorInformation);
    auto formulaMin = std::make_shared<storm::logic::RewardOperatorFormula>(subformula, std::string("derivative-min"), this->formulaOperatorInformation);

    storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<FunctionType>> simplifier(modelCopy);

    STORM_LOG_ASSERT(simplifier.simplify(*formulaMax), "Could not simplify derivative model.");
    auto modelMax = simplifier.getSimplifiedModel();

    STORM_LOG_ASSERT(simplifier.simplify(*formulaMin), "Could not simplify derivative model.");
    auto modelMin = simplifier.getSimplifiedModel();

    modelMax = storm::api::performDeterministicSparseBisimulationMinimization(modelMax, {formulaMax}, storage::BisimulationType::Strong);
    modelMin = storm::api::performDeterministicSparseBisimulationMinimization(modelMin, {formulaMin}, storage::BisimulationType::Strong);

    // std::cout << parameter << std::endl;
    // std::cout << "Model copy:" << std::endl;
    // modelCopy.writeDotToStream(std::cout, 30, true, nullptr, &modelCopy.getRewardModel("derivative-min").getStateRewardVector(),
    //                            &modelCopy.getRewardModel("derivative-max").getStateRewardVector());
    // std::cout << "Model max:" << std::endl;
    // modelMax->writeDotToStream(std::cout, 30, true, nullptr, &modelMax->getUniqueRewardModel().getStateRewardVector());
    // std::cout << "Model min:" << std::endl;
    // modelMin->writeDotToStream(std::cout, 30, true, nullptr, &modelMin->getUniqueRewardModel().getStateRewardVector());

    return std::make_pair(std::make_pair(*modelMax, *modelMin), std::make_pair(formulaMin, formulaMax));
}

template<typename FunctionType, typename ConstantType>
void DerivativeBoundFinder<FunctionType, ConstantType>::updateMonotonicityResult(
    std::vector<ConstantType> derivativeMinValues, std::vector<ConstantType> derivativeMaxValues,
    std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> localMonotonicityResult,
    VariableType<FunctionType> const& parameter, uint_fast64_t initialState) {
    std::cout << derivativeMinValues[initialState] << " <= d" << parameter << " <= " << derivativeMaxValues[initialState] << std::endl;
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

template<typename FunctionType, typename ConstantType>
models::sparse::Dtmc<FunctionType> DerivativeBoundFinder<FunctionType, ConstantType>::minimizeParameterCountInDTMC(models::sparse::Dtmc<FunctionType> dtmc) {
    storage::SparseMatrix<FunctionType> transitionMatrix = model.getTransitionMatrix();
    
    // Repeat this algorithm until we can't mimimize anything anymore
    bool somethingChanged = true;
    auto allParameters = storm::models::sparse::getAllParameters(dtmc);
    
    models::sparse::StateLabeling newLabels(dtmc.getStateLabeling());

    // Check the reward model - do not touch states with rewards
    boost::optional<std::vector<FunctionType>> stateRewardVector;
    if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
        if (currentCheckTaskNoBound->isRewardModelSet()) {
            stateRewardVector =  model.getRewardModel(currentCheckTaskNoBound->getRewardModel()).getStateRewardVector();
        } else {
            stateRewardVector =  model.getRewardModel("").getStateRewardVector();
        }
    }
    
    while (somethingChanged) {
        somethingChanged = false;
        // Search backwards from parameters to where they join their paths
        // Tally up the probabilities
        // Multiple backwards transitions not supported yet
        //
        // Maps are from state that is the last one before the parameter transition occurs
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::vector<uint_fast64_t>>> alreadyVisitedStates;
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::vector<RationalNumber>>> talliedUpProbabilities;
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, bool>> doneSearching;
        // States visited from first with p to second
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, uint_fast64_t>> increasingSuccessors;
        // States visited from first with 1-p to second
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, uint_fast64_t>> decreasingSuccessors;
        
        // We only join states that have all the same labels.
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::set<std::string>>> joiningLabels;
        
        // As I can't create my own, save the ps from the matrix to re-use them
        std::map<storm::RationalFunctionVariable, FunctionType> pFunctions;

        uint_fast64_t numberOfSearchingTransitions = 0;
        for (auto const& parameter : allParameters) {
            alreadyVisitedStates[parameter] = std::map<uint_fast64_t, std::vector<uint_fast64_t>>();
            talliedUpProbabilities[parameter] = std::map<uint_fast64_t, std::vector<RationalNumber>>();
            increasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            decreasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            joiningLabels[parameter] = std::map<uint_fast64_t, std::set<std::string>>();
        }
        // Search for all occurences of parameters and fill up the maps
        for (uint_fast64_t row = 0; row < transitionMatrix.getRowCount(); row++) {
            for (auto const& entry : transitionMatrix.getRow(row)) {
                if (entry.getValue().isConstant()) {
                    continue;
                }
                STORM_LOG_ASSERT(entry.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
                auto parameter = *entry.getValue().gatherVariables().begin();
                auto cache = entry.getValue().nominatorAsPolynomial().pCache();
                FunctionType parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
                if (entry.getValue() == parameterAsFunction) {
                    pFunctions[parameter] = entry.getValue();
                    numberOfSearchingTransitions++;
                    auto parameter = entry.getValue().nominatorAsPolynomial().getSingleVariable();
                    STORM_LOG_ASSERT(alreadyVisitedStates[parameter].count(row) == 0, "Flip minimization only supports simple pMCs.");
                    alreadyVisitedStates[parameter][row] = std::vector<uint_fast64_t>();
                    alreadyVisitedStates[parameter][row].push_back(row);
                    talliedUpProbabilities[parameter][row] = std::vector<RationalNumber>();
                    talliedUpProbabilities[parameter][row].push_back(utility::one<RationalNumber>());
                    increasingSuccessors[parameter][row] = entry.getColumn();
                    doneSearching[parameter][row] = false;
                    joiningLabels[parameter][row] = newLabels.getLabelsOfState(row);
                } else if (utility::one<FunctionType>() - entry.getValue() == parameterAsFunction) {
                    decreasingSuccessors[parameter][row] = entry.getColumn();
                } else {
                    STORM_LOG_ASSERT(false, "Flip minimization only supports simple pMCs.");
                }
            }
        }
        // Step by step, do a backwards search
        uint_fast64_t transitionsDoneSearching = 0;
        auto backwardsTransitions = transitionMatrix.transpose(true);
        // Final results of the search: Sets of sets of states that we can transform together
        // The maps are <state where the paths join> -> <set of states that leads there, with this probability>
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::map<uint_fast64_t, RationalNumber>>> joinedFirstStates;
        for (auto const& parameter : allParameters) {
            joinedFirstStates[parameter] = std::map<uint_fast64_t, std::map<uint_fast64_t, RationalNumber>>();
        }
        while (transitionsDoneSearching < numberOfSearchingTransitions) {
            for (auto const& pair : alreadyVisitedStates) {
                auto parameter = pair.first;
                auto visitedStateMap = pair.second;
                // Step one step back for parameter and firstState
                for (auto const& pair2 : visitedStateMap) {
                    auto firstState = pair2.first;
                    auto visitedStates = pair2.second;
                    if (doneSearching[parameter][firstState]) {
                        continue;
                    }
                    auto row = backwardsTransitions.getRow(visitedStates.back());
                    
                    storage::MatrixEntry<uint_fast64_t, FunctionType> entry;
                    bool entryFound = false;
                    for (auto const& loopEntry : row) {
                        if (loopEntry.getValue().isZero()) {
                            continue;
                        }
                        // We only support one backlink per state for now, so if we already found an entry, abort
                        if (entryFound) {
                            entryFound = false;
                            break;
                        }
                        // We only support all involved states to have the same label
                        if (joiningLabels[parameter][visitedStates.back()] != dtmc.getLabelsOfState(entry.getColumn())) {
                            break;
                        }
                        // We only support all involved states to have no rewards
                        if (stateRewardVector && !(*stateRewardVector)[visitedStates.back()].isZero()) {
                            break;
                        }

                        entry = loopEntry;
                        entryFound = true;
                    }
                    if (!entryFound) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                        continue;
                    }

                    // Only search until another parametric transition occurs, then we can't continue anympore
                    if (!entry.getValue().isConstant()) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                        continue;
                    }
                    auto constantProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                    auto currentState = entry.getColumn();
                    
                    // The next probability is the previous one times the considered tranistion
                    talliedUpProbabilities[parameter][firstState].push_back(talliedUpProbabilities[parameter][firstState].back() * constantProbability);
                    alreadyVisitedStates[parameter][firstState].push_back(currentState);

                    // Check if any search of the same parameter has already gotten here
                    for (auto const& pair3 : alreadyVisitedStates[parameter]) {
                        auto otherFirstState = pair3.first;
                        // Don't join a state with itself
                        if (otherFirstState == firstState) {
                            continue;
                        }
                        // If the labels don't match, this is not a match
                        if (joiningLabels[parameter][otherFirstState] != joiningLabels[parameter][firstState]) {
                            continue;
                        }
                        auto otherVisitedStates = pair3.second;
                        // TODO maybe use a set if this is too slow
                        bool done = false;
                        for (uint_fast64_t i = 0; i < otherVisitedStates.size(); i++) {
                            auto state = otherVisitedStates[i];
                            if (state == currentState) {
                                // Wow!!! Another search of the same parameter has gotten here.
                                transitionsDoneSearching++;
                                doneSearching[parameter][firstState] = true;
                                transitionsDoneSearching++;
                                doneSearching[parameter][otherFirstState] = true;
                                
                                if (!joinedFirstStates[parameter].count(state)) {
                                    joinedFirstStates[parameter][state] = std::map<uint_fast64_t, RationalNumber>();
                                }
                                RationalNumber p1 = talliedUpProbabilities[parameter][firstState].back();
                                joinedFirstStates[parameter][state][firstState] = p1;
                                RationalNumber p2 = talliedUpProbabilities[parameter][otherFirstState][i];
                                joinedFirstStates[parameter][state][otherFirstState] = p2;

                                // We are all done. We can continue the search for the next, unrelated thing.
                                done = true;
                                break;
                            }
                        }
                        if (done) {
                            break;
                        }
                    }
                }
            }
        }
        
        // Check if we found any parameters that we can join into one transition, and if so, do it
        for (auto const& parameter : allParameters) {
            if (joinedFirstStates[parameter].size() > 0) {
                somethingChanged = true;
                break;
            }
        }
        if (!somethingChanged) {
            break;
        }

        storage::SparseMatrix<FunctionType> wipMatrix(transitionMatrix); 
        for (auto const& parameter : allParameters) {
            // We add three new states to the DTMC for this procedure for every state where parameters join
            if (joinedFirstStates[parameter].size() == 0) {
                continue;
            }
            std::cout << "joining " << joinedFirstStates[parameter].size() << " cases for " << parameter << std::endl;
            uint_fast64_t runningCounter = wipMatrix.getRowCount();

            std::map<uint_fast64_t, std::vector<uint_fast64_t>> newStates;
            for (auto const& pair : joinedFirstStates[parameter]) {
                // The state that will come before the parametric flip
                auto firstNewState = runningCounter;
                // The state after p 
                auto secondNewState = runningCounter + 1;
                // The state after 1-p
                auto thirdNewState = runningCounter + 2;
                newStates[pair.first] = {firstNewState, secondNewState, thirdNewState};
                runningCounter += 3;
            }
            
            // Extend labeling to more states (Found no better way to do it)
            models::sparse::StateLabeling nextNewLabels(runningCounter);
            for (auto const& label : newLabels.getLabels()) {
                nextNewLabels.addLabel(label);
            }
            for (uint_fast64_t state = 0; state < wipMatrix.getRowCount(); state++) {
                for (auto const& label : newLabels.getLabelsOfState(state)) {
                    nextNewLabels.addLabelToState(label, state);
                }
            }
            newLabels = nextNewLabels;

            storage::SparseMatrixBuilder<FunctionType> builder;
            for (uint_fast64_t row = 0; row < wipMatrix.getRowCount(); row++) {
                // If from row, we are starting with a transition we will flip...
                if (joinedFirstStates[parameter].count(row)) {
                    // We will join until these transitions
                    auto statesBeforeParams = joinedFirstStates[parameter][row];
                    RationalNumber summedProbability = 0;
                    for (auto const& entry2 : statesBeforeParams) {
                        summedProbability += entry2.second;
                    }
                    // Add the rest of the transitions before, because we need to preserve order
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        // If this state was (the second last) visited in any of the joinedFirstStates, do not add it
                        // Else add it
                        bool addState = true;
                        for (auto const& pair : statesBeforeParams) {
                            auto id = pair.first;
                            auto visitedStates = alreadyVisitedStates[parameter][id];
                            if (std::find(visitedStates.begin(), visitedStates.end(), entry.getColumn()) != visitedStates.end()) {
                                // This state was visited in id's path and we have thus joined it above, so don't add it
                                addState = false;
                            }
                        }
                        if (addState) {
                            builder.addNextValue(row, entry.getColumn(), entry.getValue());
                        }
                    }
                    builder.addNextValue(row, newStates[row][0], utility::convertNumber<FunctionType>(summedProbability));
                } else { // Else, just add the stuff back into the matrix
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        builder.addNextValue(row, entry.getColumn(), entry.getValue());
                    }
                }
            }
            
            // FunctionType parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
            FunctionType parameterAsFunction = pFunctions[parameter];

            for (auto const& pair : joinedFirstStates[parameter]) {
                auto row = pair.first;
                RationalNumber summedProbability = 0;
                for (auto const& entry2 : joinedFirstStates[parameter][row]) {
                    summedProbability += entry2.second;
                }
                // Add parametric transitions 
                builder.addNextValue(newStates[row][0], newStates[row][1], parameterAsFunction);
                builder.addNextValue(newStates[row][0], newStates[row][2], utility::one<FunctionType>() - parameterAsFunction);
                
                // Add transitions from p and 1-p states to actual successors
                std::priority_queue<std::pair<uint_fast64_t, FunctionType>, std::vector<std::pair<uint_fast64_t, FunctionType>>, std::greater<std::pair<uint_fast64_t, FunctionType>>> insertInOrder1;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto increasingSuccessor = increasingSuccessors[parameter][id];
                    insertInOrder1.emplace(std::make_pair(increasingSuccessor, utility::convertNumber<FunctionType>(probability)));
                }
                while (!insertInOrder1.empty()) {
                    auto pair = insertInOrder1.top();
                    insertInOrder1.pop();
                    builder.addNextValue(newStates[row][1], pair.first, pair.second);
                }
                
                std::priority_queue<std::pair<uint_fast64_t, FunctionType>, std::vector<std::pair<uint_fast64_t, FunctionType>>, std::greater<std::pair<uint_fast64_t, FunctionType>>> insertInOrder2;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto decreasingSuccessor = decreasingSuccessors[parameter][id];
                    insertInOrder2.emplace(std::make_pair(decreasingSuccessor, utility::convertNumber<FunctionType>(probability)));
                }
                while (!insertInOrder2.empty()) {
                    auto pair = insertInOrder2.top();
                    insertInOrder2.pop();
                    builder.addNextValue(newStates[row][2], pair.first, pair.second);
                }
                
                for (auto const& label : joiningLabels[parameter][row]) {
                    if (!newLabels.getLabels().count(label)) {
                        newLabels.addLabel(label);
                    }
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        newLabels.addLabelToState(label, newStates[row][i]);
                    }
                }
            }
            wipMatrix = builder.build();
        }
        transitionMatrix = wipMatrix;
    }
    // models::sparse::StateLabeling stateLabeling(transitionMatrix.getRowCount());
    // for (auto const& label : dtmc.getStateLabeling().getLabels()) {
    //     stateLabeling.addLabel(label);
    // }
    // for (uint_fast64_t state = 0; state < dtmc.getNumberOfStates(); state++) {
    //     for (auto const& label : dtmc.getStateLabeling().getLabelsOfState(state)) {
    //         stateLabeling.addLabelToState(label, state);
    //     }
    // }
    // for (uint_fast64_t state = dtmc.getNumberOfStates(); state < transitionMatrix.getRowCount(); state++) {
    //     for (auto const& label : newLabels.getLabelsOfState(state)) {
    //         stateLabeling.addLabelToState(label, state);
    //     }
    // }
    models::sparse::Dtmc<FunctionType> newDTMC(transitionMatrix, newLabels);
    for (auto const& rewardModel : dtmc.getRewardModels()) {
        newDTMC.addRewardModel(rewardModel.first, rewardModel.second);
    }
    return newDTMC;

    // Do a DFS and search for occuring parameters, stop there
    // storage::BitVector currentStates = model.getInitialStates();
    // storage::BitVector visitedStates(currentStates.size());
    // std::vector<FunctionType> constantProbabilities(currentStates.size(), utility::one<FunctionType>());
    // while (!currentStates.empty()) {
    //     for (uint_fast64_t state : currentStates) {
    //         storage::BitVector additionalStates;
    //         for (storage::MatrixEntry<uint_fast64_t, FunctionType> const& entry : transitionMatrix.getRow(state)) {
    //             FunctionType value = entry.getValue();
    //             uint_fast64_t nextState = entry.getColumn();
    //             if (value.isConstant() && !visitedStates.get(nextState)) {
    //                 additionalStates.set(nextState, true);
    //                 constantProbabilities[nextState] *=
    //             }
    //         }
    //         currentStates.set(state, false);
    //         visitedStates.set(state, true);
    //         currentStates |= additionalStates;
    //     }
    // }
    // What have you found? Maybe some nice transitions that we can join?
    

    // Make a DTMC that is both reversed and missing the parametric transitions
    // auto const reverseTransitions = dtmc.getBackwardTransitions();
    // storage::SparseMatrixBuilder<RationalNumber> reverseTransitionsWithDeletedParametersBuilder;
    // for (uint_fast64_t row = 0; row < reverseTransitions.getRowCount(); row++) {
    //     for (storage::MatrixEntry<uint_fast64_t, FunctionType> const& entry : reverseTransitions.getRow(row)) {
    //         if (entry.getValue().isConstant()) {
    //             reverseTransitionsWithDeletedParametersBuilder.addNextValue(row, entry.getColumn(), utility::convertNumber<RationalNumber>(entry.getValue()));
    //         }
    //     }
    // }
    // auto reverseTransitionsWithDeletedParametersMatrix = reverseTransitionsWithDeletedParametersBuilder.build();
    // models::sparse::Dtmc<RationalNumber> reverseTransitionsWithDeletedParameters(reverseTransitionsWithDeletedParametersMatrix, dtmc.getStateLabeling());
    
    // auto subformulaConstructor = std::make_shared<logic::AtomicLabelFormula>("init");
    // auto subformula = std::make_shared<logic::EventuallyFormula>(subformulaConstructor, logic::FormulaContext::Probability, boost::none);
    // auto formula = std::make_shared<storm::logic::ProbabilityOperatorFormula>(subformula, this->formulaOperatorInformation);
    // auto checkTask = std::make_shared<storm::modelchecker::CheckTask<storm::logic::Formula, RationalNumber>>(*formula);
    
    // modelchecker::SparseDtmcPrctlModelChecker<models::sparse::Dtmc<RationalNumber>> modelChecker(reverseTransitionsWithDeletedParameters);
    // auto results = modelChecker.check(Environment(), *checkTask);
    // auto resultVector = results->template asExplicitQuantitativeCheckResult<RationalNumber>().getValueVector();
    // for (auto const& entry : resultVector) {
    //     std::cout << entry << std::endl;
    // }
}

template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
template class DerivativeBoundFinder<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm
