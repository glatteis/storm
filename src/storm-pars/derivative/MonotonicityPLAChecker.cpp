#include "MonotonicityPLAChecker.h"
#include <_types/_uint64_t.h>
#include <cstdint>
#include <memory>
#include "models/ModelBase.h"
#include "models/sparse/Dtmc.h"
#include "settings/modules/GeneralSettings.h"
#include "storm-pars/utility/parametric.h"

namespace storm {
namespace derivative {

template<typename FunctionType>
using VariableType = typename utility::parametric::VariableType<FunctionType>::type;
template<typename FunctionType>
using CoefficientType = typename utility::parametric::CoefficientType<FunctionType>::type;

template<typename FunctionType, typename ConstantType>
std::pair<std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>,
          std::unique_ptr<storm::modelchecker::QuantitativeCheckResult<ConstantType>>>
MonotonicityPLAChecker<FunctionType, ConstantType>::getDerivativeBound(Environment const& env, storm::storage::ParameterRegion<FunctionType> const& currRegion,
                                                                       VariableType<FunctionType> parameter) {
    auto initialState = getInitialState();
    // model.writeDotToStream(std::cout);
    modelChecker.specify(env, std::make_shared<models::sparse::Dtmc<FunctionType>>(model), *currentCheckTask, false, false);
    std::vector<ConstantType> minBound = modelChecker.getBound(env, currRegion, storm::solver::OptimizationDirection::Minimize, nullptr)
                                             ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                             .getValueVector();
    std::vector<ConstantType> maxBound = modelChecker.getBound(env, currRegion, storm::solver::OptimizationDirection::Maximize, nullptr)
                                             ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                             .getValueVector();
    auto derivativeCheckStuff = boundFinder.computeMonotonicityTasks(env, currRegion, minBound, maxBound, nullptr, parameter);
    auto modelMax = derivativeCheckStuff.first.first;
    auto modelMin = derivativeCheckStuff.first.second;
    auto formulaMin = derivativeCheckStuff.second.first;
    auto formulaMax = derivativeCheckStuff.second.second;

    // std::cout << "min bound: " << minBound[initialState] << ", max bound: " << maxBound[initialState] << std::endl;
    // std::cout << *formulaMin << std::endl;
    // std::cout << *formulaMax << std::endl;

    auto checkTaskMax = std::make_shared<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);
    auto checkTaskMin = std::make_shared<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMin);
    modelChecker.specify(env, std::make_shared<models::sparse::Dtmc<FunctionType>>(modelMax), *checkTaskMax, false, false);
    auto derivativeResultsMax = modelChecker.getBound(env, currRegion, OptimizationDirection::Maximize, nullptr)
                                    ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                    .getValueVector();
    modelChecker.specify(env, std::make_shared<models::sparse::Dtmc<FunctionType>>(modelMin), *checkTaskMin, false, false);
    auto derivativeResultsMin = modelChecker.getBound(env, currRegion, OptimizationDirection::Minimize, nullptr)
                                    ->template asExplicitQuantitativeCheckResult<ConstantType>()
                                    .getValueVector();
    // std::cout << "der. min bound: " << derivativeResultsMin[initialState] << ", der. max bound: " << derivativeResultsMax[initialState] << std::endl;
    STORM_LOG_INFO("Derivative monotonicity result computed for " << parameter);
    auto resultMax = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeResultsMax);
    auto resultMin = std::make_unique<modelchecker::ExplicitQuantitativeCheckResult<ConstantType>>(derivativeResultsMin);

    return std::make_pair(std::move(resultMax), std::move(resultMin));
}

template<typename FunctionType, typename ConstantType>
uint_fast64_t MonotonicityPLAChecker<FunctionType, ConstantType>::getInitialState() {
    auto internalModel = boundFinder.getInternalModel();
    uint_fast64_t initialState;
    const storm::storage::BitVector initialVector = internalModel.getInitialStates();
    for (uint_fast64_t x : initialVector) {
        initialState = x;
        break;
    }
    return initialState;
}

template<typename FunctionType, typename ConstantType>
void MonotonicityPLAChecker<FunctionType, ConstantType>::performMonotonicityPLA(Environment const& env, VariableType<FunctionType> wrt,
                                                                                ConstantType terminateArea) {
    STORM_LOG_INFO("Demo derivative PLA w.r.t. " << wrt);
    auto positivelyMonotoneArea = storm::utility::zero<ConstantType>();
    auto negativelyMonotoneArea = storm::utility::zero<ConstantType>();
    auto unknownArea = storm::utility::zero<ConstantType>();
    std::vector<storage::ParameterRegion<FunctionType>> positivelyMonotoneRegions;
    std::vector<storage::ParameterRegion<FunctionType>> negativelyMonotoneRegions;
    std::vector<storage::ParameterRegion<FunctionType>> unknownRegions;
    const ConstantType precision = utility::convertNumber<ConstantType>(storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision());
    std::queue<storage::ParameterRegion<FunctionType>> regionQueue;
    // Make big region
    std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigLower;
    std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> bigUpper;
    for (auto const& parameter : storm::models::sparse::getAllParameters(model)) {
        bigLower[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(0.01);
        bigUpper[parameter] = utility::convertNumber<CoefficientType<FunctionType>>(0.99);
    }
    storage::ParameterRegion<FunctionType> bigRegion(bigLower, bigUpper);
    
    auto initialState = getInitialState();
    
    const uint_fast64_t maxIterations = 10000;

    regionQueue.push(bigRegion);

    uint_fast64_t regionsComputed = 0;
    while (!regionQueue.empty() && (positivelyMonotoneArea + negativelyMonotoneArea) < 1 - terminateArea && maxIterations < 1000) {
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

template class MonotonicityPLAChecker<RationalFunction, RationalNumber>;
template class MonotonicityPLAChecker<RationalFunction, double>;
}  // namespace derivative
}  // namespace storm