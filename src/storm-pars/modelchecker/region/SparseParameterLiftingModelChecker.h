#pragma once

#include "storm-pars/modelchecker/region/RegionModelChecker.h"
#include "storm-pars/modelchecker/instantiation/SparseInstantiationModelChecker.h"
#include "storm-pars/storage/ParameterRegion.h"
#include "storm-pars/utility/parametric.h"

#include "storm/logic/Formulas.h"
#include "storm/modelchecker/CheckTask.h"
#include "storm/modelchecker/results/CheckResult.h"
#include "storm/solver/OptimizationDirection.h"

namespace storm {
    namespace modelchecker {
            
        /*!
         * Class to approximatively check a formula on a parametric model for all parameter valuations within a region
         * It is assumed that all considered valuations are graph-preserving and well defined, i.e.,
         * * all non-const transition probabilities evaluate to some non-zero value
         * * the sum of all outgoing transitions is one
         */
        template <typename SparseModelType, typename ConstantType>
        class SparseParameterLiftingModelChecker : public RegionModelChecker<typename SparseModelType::ValueType> {
        public:
            SparseParameterLiftingModelChecker();
            virtual ~SparseParameterLiftingModelChecker() = default;
            

            /*!
             * Analyzes the given region by means of parameter lifting.
             */
            virtual RegionResult analyzeRegion(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, RegionResultHypothesis const& hypothesis = RegionResultHypothesis::Unknown, RegionResult const& initialResult = RegionResult::Unknown, bool sampleVerticesOfRegion = false, storm::analysis::Order* reachabilityOrder = nullptr, std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename RegionModelChecker<typename SparseModelType::ValueType>::VariableType>> localMonotonicityResult = nullptr) override;

            /*!
             * Analyzes the 2^#parameters corner points of the given region.
             */
            RegionResult sampleVertices(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, RegionResult const& initialResult = RegionResult::Unknown);
            
            /*!
             * Checks the specified formula on the given region by applying parameter lifting (Parameter choices are lifted to nondeterministic choices)
             * This yields a (sound) approximative model checking result.

             * @param region the region on which parameter lifting is applied
             * @param dirForParameters  The optimization direction for the parameter choices. If this is, e.g., minimize, then the returned result will be a lower bound for all results induced by the parameter evaluations inside the region.
             */
            std::unique_ptr<CheckResult> check(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, storm::solver::OptimizationDirection const& dirForParameters, storm::analysis::Order* reachabilityOrder = nullptr, std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename RegionModelChecker<typename SparseModelType::ValueType>::VariableType>> localMonotonicityResult = nullptr);
            
            std::unique_ptr<QuantitativeCheckResult<ConstantType>> getBound(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, storm::solver::OptimizationDirection const& dirForParameters);
            virtual typename SparseModelType::ValueType getBoundAtInitState(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, storm::solver::OptimizationDirection const& dirForParameters) override;

            
            /*!
             * Finds the extremal value within the given region and with the given precision.
             * The returned value v corresponds to the value at the returned valuation.
             * The actual maximum (minimum) lies in the interval [v, v+precision] ([v-precision, v])
             */
            virtual std::pair<typename SparseModelType::ValueType, typename storm::storage::ParameterRegion<typename SparseModelType::ValueType>::Valuation> computeExtremalValue(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, storm::solver::OptimizationDirection const& dirForParameters, typename SparseModelType::ValueType const& precision) override;
            
            SparseModelType const& getConsideredParametricModel() const;
            CheckTask<storm::logic::Formula, ConstantType> const& getCurrentCheckTask() const;
            
        protected:
            void specifyFormula(Environment const& env, CheckTask<storm::logic::Formula, typename SparseModelType::ValueType> const& checkTask);
            
            // Resets all data that correspond to the currently defined property.
            virtual void reset() = 0;

            virtual void specifyBoundedUntilFormula(Environment const& env, CheckTask<storm::logic::BoundedUntilFormula, ConstantType> const& checkTask);
            virtual void specifyUntilFormula(Environment const& env, CheckTask<storm::logic::UntilFormula, ConstantType> const& checkTask);
            virtual void specifyReachabilityProbabilityFormula(Environment const& env, CheckTask<storm::logic::EventuallyFormula, ConstantType> const& checkTask);
            virtual void specifyReachabilityRewardFormula(Environment const& env, CheckTask<storm::logic::EventuallyFormula, ConstantType> const& checkTask);
            virtual void specifyCumulativeRewardFormula(Environment const& env, CheckTask<storm::logic::CumulativeRewardFormula, ConstantType> const& checkTask);

            
            virtual storm::modelchecker::SparseInstantiationModelChecker<SparseModelType, ConstantType>& getInstantiationChecker() = 0;
            
            virtual std::unique_ptr<CheckResult> computeQuantitativeValues(Environment const& env, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, storm::solver::OptimizationDirection const& dirForParameters, storm::analysis::Order* reachabilityOrder = nullptr, std::shared_ptr<storm::analysis::LocalMonotonicityResult<typename RegionModelChecker<typename SparseModelType::ValueType>::VariableType>> localMonotonicityResult = nullptr) = 0;
            
            
            std::shared_ptr<SparseModelType> parametricModel;
            std::unique_ptr<CheckTask<storm::logic::Formula, ConstantType>> currentCheckTask;

        private:
            // store the current formula. Note that currentCheckTask only stores a reference to the formula.
            std::shared_ptr<storm::logic::Formula const> currentFormula;
            
        };
    }
}
