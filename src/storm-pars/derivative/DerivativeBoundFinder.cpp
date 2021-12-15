#include "DerivativeBoundFinder.h"
#include <_types/_uint64_t.h>
#include <cstdint>
#include <memory>
#include "logic/FormulasForwardDeclarations.h"
#include "settings/modules/GeneralSettings.h"
#include "storm-pars/analysis/MonotonicityResult.h"
#include "storm-pars/utility/parametric.h"
#include "storm-parsers/parser/FormulaParser.h"
#include "utility/graph.h"

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
std::pair<models::sparse::Dtmc<FunctionType>, std::pair<std::shared_ptr<storm::logic::Formula>,
                                                        std::shared_ptr<storm::logic::Formula>>>
DerivativeBoundFinder<FunctionType, ConstantType>::computeMonotonicityTasks(
    Environment const& env, storm::storage::ParameterRegion<FunctionType> const& region, std::vector<ConstantType> minValues,
    std::vector<ConstantType> maxValues,
    std::shared_ptr<storm::analysis::LocalMonotonicityResult<VariableType<FunctionType>>> localMonotonicityResult,
    VariableType<FunctionType> const& parameter
) {

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
        
        if (derivativeMax < 0) {
            localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Decr);
        } else if (derivativeMin > 0) {
            localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Incr);
        } else  if (derivativeMin == 0 && derivativeMax == 0) {
            localMonotonicityResult->setMonotonicity(i, parameter, analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Constant);
        }

        stateRewardsMax[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivativeMax)));
        stateRewardsMin[i] += FunctionType(utility::convertNumber<CoefficientType<FunctionType>>((ConstantType)(derivativeMin)));
    }

    models::sparse::StandardRewardModel<FunctionType> rewardModelMax(std::move(stateRewardsMax));
    models::sparse::StandardRewardModel<FunctionType> rewardModelMin(std::move(stateRewardsMin));

    modelCopy.addRewardModel("derivative-max", rewardModelMax);
    modelCopy.addRewardModel("derivative-min", rewardModelMin);

    storage::BitVector target = modelCopy.getStates("target");
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
    
    return std::make_pair(modelCopy, std::make_pair(formulaMin, formulaMax));
}
        
template<typename FunctionType, typename ConstantType>
void
DerivativeBoundFinder<FunctionType, ConstantType>::updateMonotonicityResult(
    std::vector<ConstantType> derivativeMinValues, std::vector<ConstantType> derivativeMaxValues,
    std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename utility::parametric::VariableType<FunctionType>::type>> localMonotonicityResult,
    VariableType<FunctionType> const& parameter,
    uint_fast64_t initialState
) {
    boost::optional<typename analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity> finalResult;
    if (derivativeMaxValues[initialState] < 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Decr;
    } else if (derivativeMinValues[initialState] > 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Incr;
    } else if (derivativeMaxValues[initialState] == 0 && derivativeMinValues[initialState] == 0) {
        finalResult = analysis::MonotonicityResult<VariableType<FunctionType>>::Monotonicity::Constant;
    }

    if (finalResult) {
        localMonotonicityResult->getGlobalMonotonicityResult()->addMonotonicityResult(parameter, *finalResult);
        for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
            localMonotonicityResult->setMonotonicity(i, parameter, *finalResult);
        }
    }
}

// template<typename FunctionType, typename ConstantType>
// void DerivativeBoundFinder<FunctionType, ConstantType>::derivativePLASketch(Environment const& env, VariableType<FunctionType> wrt,
//                                                                             ConstantType terminateArea) {
//     STORM_LOG_INFO("Demo derivative PLA w.r.t. " << wrt);
//     auto positivelyMonotoneArea = storm::utility::zero<ConstantType>();
//     auto negativelyMonotoneArea = storm::utility::zero<ConstantType>();
//     auto unknownArea = storm::utility::zero<ConstantType>();
//     std::vector<storage::ParameterRegion<FunctionType>> positivelyMonotoneRegions;
//     std::vector<storage::ParameterRegion<FunctionType>> negativelyMonotoneRegions;
//     std::vector<storage::ParameterRegion<FunctionType>> unknownRegions;
//     uint_fast64_t initialState;
//     const storm::storage::BitVector initialVector = model.getInitialStates();
//     for (uint_fast64_t x : initialVector) {
//         initialState = x;
//         break;
//     }
//     const ConstantType precision = utility::convertNumber<ConstantType>(storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision());
//     std::queue<storage::ParameterRegion<FunctionType>> regionQueue;
//     // Make big region
//     std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigLower;
//     std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigUpper;
//     for (auto const& parameter : storm::models::sparse::getAllParameters(model)) {
//         bigLower[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(0.01);
//         bigUpper[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(0.99);
//     }
//     storage::ParameterRegion<FunctionType> bigRegion(bigLower, bigUpper);

//     regionQueue.push(bigRegion);
//     uint_fast64_t regionsComputed = 0;
//     while (!regionQueue.empty() && (positivelyMonotoneArea + negativelyMonotoneArea) < 1 - terminateArea && regionsComputed < 1000) {
//         regionsComputed++;
//         storage::ParameterRegion<FunctionType> currRegion = regionQueue.front();
//         regionQueue.pop();
//         STORM_LOG_INFO("Looking at region: " << currRegion);
//         auto results = getDerivativeBound(env, currRegion, wrt);
//         auto resultMax = std::move(results.first)->template asExplicitQuantitativeCheckResult<ConstantType>().getValueVector();
//         auto resultMin = std::move(results.second->template asExplicitQuantitativeCheckResult<ConstantType>().getValueVector());
//         if (resultMin[initialState] > precision) {
//             std::cout << "+";
//             STORM_LOG_INFO("Found positively monotone region " << currRegion);
//             positivelyMonotoneArea += utility::convertNumber<ConstantType>(currRegion.area());
//             positivelyMonotoneRegions.push_back(currRegion);
//         } else if (resultMax[initialState] < -precision) {
//             std::cout << "-";
//             STORM_LOG_INFO("Found negatively monotone region " << currRegion);
//             negativelyMonotoneArea += utility::convertNumber<ConstantType>(currRegion.area());
//             negativelyMonotoneRegions.push_back(currRegion);
//         } else {
//             std::cout << "?";
//             STORM_LOG_INFO("Splitting region " << currRegion);
//             std::vector<storm::storage::ParameterRegion<FunctionType>> newRegions;
//             currRegion.split(currRegion.getCenterPoint(), newRegions);
//             for (auto const& region : newRegions) {
//                 regionQueue.emplace(region);
//             }
//         }
//         std::cout << std::flush;
//         STORM_LOG_INFO("Positively monotone area: " << positivelyMonotoneArea);
//         STORM_LOG_INFO("Negatively monotone area: " << negativelyMonotoneArea);
//     }
//     std::cout << std::endl;

//     while (!regionQueue.empty()) {
//         auto region = regionQueue.front();
//         regionQueue.pop();
//         unknownRegions.push_back(region);
//         unknownArea += utility::convertNumber<ConstantType>(region.area());
//     }

//     std::cout << "Positively monotone area: " << positivelyMonotoneArea << std::endl;
//     std::cout << "Negatively monotone area: " << negativelyMonotoneArea << std::endl;
//     std::cout << "Unknown area: " << unknownArea << std::endl;

//     std::vector<std::pair<storm::storage::ParameterRegion<FunctionType>, storm::modelchecker::RegionResult>> regionResults;

//     for (auto const& region : positivelyMonotoneRegions) {
//         regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::AllSat));
//     }
//     for (auto const& region : negativelyMonotoneRegions) {
//         regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::AllViolated));
//     }
//     for (auto const& region : unknownRegions) {
//         regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::Unknown));
//     }

//     std::cout << "Quick visualization hack, AllSat means positively monotone and AllViolated means negatively monotone :-)" << std::endl;
//     modelchecker::RegionRefinementCheckResult<FunctionType> regionRefinementCheckResult(regionResults, bigRegion);
//     regionRefinementCheckResult.writeToStream(std::cout);
//     regionRefinementCheckResult.writeIllustrationToStream(std::cout);
// }

template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
template class DerivativeBoundFinder<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm
