#include "DerivativeBoundFinder.h"
#include "settings/modules/GeneralSettings.h"
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
        if (currentCheckTaskNoBound->getFormula().isRewardOperatorFormula()) {
            auto rewardModel = currentCheckTaskNoBound->getFormula().asRewardOperatorFormula().getRewardModelName();
            stateRewardsMax[i] = model.getRewardModel(rewardModel).getStateRewardVector()[i];
            stateRewardsMin[i] = model.getRewardModel(rewardModel).getStateRewardVector()[i];
        } else {
            stateRewardsMax[i] = utility::zero<FunctionType>();
            stateRewardsMin[i] = utility::zero<FunctionType>();
        }

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

    /* std::cout << "max at init: " << derivativeMax[initialState] << std::endl; */
    /* std::cout << "min at init: " << derivativeMin[initialState] << std::endl; */

    auto resultMax = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMax);
    auto resultMin = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeMin);

    return std::make_pair(std::move(resultMax), std::move(resultMin));
}

template<typename FunctionType, typename ConstantType>
void DerivativeBoundFinder<FunctionType, ConstantType>::derivativePLASketch(Environment const& env, VariableType<FunctionType> wrt, ConstantType terminateArea) {
    STORM_LOG_INFO("Demo derivative PLA w.r.t. " << wrt);
    auto positivelyMonotoneArea = storm::utility::zero<ConstantType>();
    auto negativelyMonotoneArea = storm::utility::zero<ConstantType>();
    auto unknownArea = storm::utility::zero<ConstantType>();
    std::vector<storage::ParameterRegion<FunctionType>> positivelyMonotoneRegions;
    std::vector<storage::ParameterRegion<FunctionType>> negativelyMonotoneRegions;
    std::vector<storage::ParameterRegion<FunctionType>> unknownRegions;
    uint_fast64_t initialState;
    const storm::storage::BitVector initialVector = model.getInitialStates();
    for (uint_fast64_t x : initialVector) {
        initialState = x;
        break;
    }
    const ConstantType precision =
        utility::convertNumber<ConstantType>(storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision());
    std::queue<storage::ParameterRegion<FunctionType>> regionQueue;
    // Make big region 
    std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigLower;
    std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigUpper;
    for (auto const& parameter : storm::models::sparse::getAllParameters(model)) {
        bigLower[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(1e-6);
        bigUpper[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(1 - 1e-6);
    }
    storage::ParameterRegion<FunctionType> bigRegion(bigLower, bigUpper);

    regionQueue.push(bigRegion);
    uint_fast64_t regionsComputed = 0;
    while (!regionQueue.empty() && (positivelyMonotoneArea + negativelyMonotoneArea) < 1 - terminateArea && regionsComputed < 1000) {
        regionsComputed++;
        storage::ParameterRegion<FunctionType> currRegion = regionQueue.front();
        regionQueue.pop();
        STORM_LOG_INFO("Looking at region: " << currRegion);
        auto results = getDerivativeBound(env, currRegion, wrt);
        auto resultMax = std::move(results.first)->template asExplicitQuantitativeCheckResult<ConstantType>().getValueVector();
        auto resultMin = std::move(results.second->template asExplicitQuantitativeCheckResult<ConstantType>().getValueVector());
        if (resultMin[initialState] > precision) {
            std::cout << "+";
            STORM_LOG_INFO("Found positively monotone region " << currRegion);
            positivelyMonotoneArea += utility::convertNumber<ConstantType>(currRegion.area());
            positivelyMonotoneRegions.push_back(currRegion);
        } else if (resultMax[initialState] < -precision) {
            std::cout << "-";
            STORM_LOG_INFO("Found negatively monotone region " << currRegion);
            negativelyMonotoneArea += utility::convertNumber<ConstantType>(currRegion.area());
            negativelyMonotoneRegions.push_back(currRegion);
        } else {
            std::cout << "?";
            STORM_LOG_INFO("Splitting region " << currRegion);
            std::vector<storm::storage::ParameterRegion<FunctionType>> newRegions;
            currRegion.split(currRegion.getCenterPoint(), newRegions);
            for (auto const& region : newRegions) {
                regionQueue.emplace(region);
            }
        }
        std::cout << std::flush;
        STORM_LOG_INFO("Positively monotone area: " << positivelyMonotoneArea);
        STORM_LOG_INFO("Negatively monotone area: " << negativelyMonotoneArea);
    }
    std::cout << std::endl;

    while (!regionQueue.empty()) {
        auto region = regionQueue.front();
        regionQueue.pop();
        unknownRegions.push_back(region);
        unknownArea += utility::convertNumber<ConstantType>(region.area());
    }

    std::cout << "Positively monotone area: " << positivelyMonotoneArea << std::endl;
    std::cout << "Negatively monotone area: " << negativelyMonotoneArea << std::endl;
    std::cout << "Unknown area: " << unknownArea << std::endl;

    std::vector<std::pair<storm::storage::ParameterRegion<FunctionType>, storm::modelchecker::RegionResult>> regionResults;

    for (auto const& region : positivelyMonotoneRegions) {
        regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::AllSat));
    }
    for (auto const& region : negativelyMonotoneRegions) {
        regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::AllViolated));
    }
    for (auto const& region : unknownRegions) {
        regionResults.push_back(std::make_pair(region, modelchecker::RegionResult::Unknown));
    }

    std::cout << "Quick visualization hack, AllSat means positively monotone and AllViolated means negatively monotone :-)" << std::endl;
    modelchecker::RegionRefinementCheckResult<FunctionType> regionRefinementCheckResult(regionResults, bigRegion);
    regionRefinementCheckResult.writeToStream(std::cout);
    regionRefinementCheckResult.writeIllustrationToStream(std::cout);
}

template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
template class DerivativeBoundFinder<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm
