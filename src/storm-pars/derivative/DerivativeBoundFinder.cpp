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
        void DerivativeBoundFinder<FunctionType, ConstantType>::liftingTest(Environment const& env) {
            std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> lower;
            std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> upper;
            carl::Variable parameter;
            for (VariableType<FunctionType> loopParam : storm::models::sparse::getAllParameters(model)) {
                lower[loopParam] = storm::utility::convertNumber<CoefficientType<FunctionType>>(0.1);
                upper[loopParam] = storm::utility::convertNumber<CoefficientType<FunctionType>>(0.9);
                parameter = loopParam;
            }

            storage::ParameterRegion<FunctionType> region(lower, upper);

            this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(model), *this->currentCheckTaskNoBound);
            auto min = liftingModelChecker->getBound(env, region, OptimizationDirection::Minimize)->template asExplicitQuantitativeCheckResult<ConstantType>();
            auto max = liftingModelChecker->getBound(env, region, OptimizationDirection::Maximize)->template asExplicitQuantitativeCheckResult<ConstantType>();

            std::cout << "w.r.t. " << parameter << std::endl;
            const storage::SparseMatrix<FunctionType> transitionMatrix = model.getTransitionMatrix();
            models::sparse::Dtmc<FunctionType> modelCopy = model;

            std::vector<FunctionType> stateRewardsMax(transitionMatrix.getRowCount());
            std::vector<FunctionType> stateRewardsMin(transitionMatrix.getRowCount());

            for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) {
                for (storage::MatrixEntry<uint_fast64_t, FunctionType> const& entry : transitionMatrix.getRow(i)) {
                    FunctionType derivative = entry.getValue().derivative(parameter);
                    if (derivative != utility::zero<FunctionType>()) {
                        ConstantType derivative = utility::convertNumber<ConstantType>(derivative);
                        ConstantType extremalValueMin = 0;
                        ConstantType extremalValueMax = 0;
                        if (derivative < 0) {
                            extremalValueMax = min[i];
                            extremalValueMin = max[i];
                        } else if (derivative > 0) {
                            extremalValueMax = max[i];
                            extremalValueMin = min[i];
                        }
                        // TODO only works for probs
                        stateRewardsMax[i] = utility::convertNumber<FunctionType>(derivative) * utility::convertNumber<FunctionType>(extremalValueMax);
                        stateRewardsMin[i] = utility::convertNumber<FunctionType>(derivative) * utility::convertNumber<FunctionType>(extremalValueMin);
                    }
                }
            }

            models::sparse::StandardRewardModel<FunctionType> rewardModelMax(std::move(stateRewardsMax));
            models::sparse::StandardRewardModel<FunctionType> rewardModelMin(std::move(stateRewardsMin));

            modelCopy.addRewardModel("derivative-max", rewardModelMax);
            modelCopy.addRewardModel("derivative-min", rewardModelMin);
            
            /* std::shared_ptr<const storm::logic::Formula> formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula, std::string("derivative-max"), this->formulaOperatorInformation, logic::RewardMeasureType::Expectation)->asRewardOperatorFormula().asSharedPointer(); */
            auto formulaMax = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula, std::string("derivative-max"), this->formulaOperatorInformation);
            auto formulaMin = std::make_shared<storm::logic::RewardOperatorFormula>(this->subformula, std::string("derivative-min"), this->formulaOperatorInformation);
            auto checkTaskMin = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);
            auto checkTaskMax = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*formulaMax);
            std::cout << checkTaskMax->getFormula().isInFragment(storm::logic::reachability().setRewardOperatorsAllowed(true).setReachabilityRewardFormulasAllowed(true).setBoundedUntilFormulasAllowed(true).setCumulativeRewardFormulasAllowed(true).setStepBoundedCumulativeRewardFormulasAllowed(true).setTimeBoundedCumulativeRewardFormulasAllowed(true).setTimeBoundedUntilFormulasAllowed(true).setStepBoundedUntilFormulasAllowed(true).setTimeBoundedUntilFormulasAllowed(true)) << std::endl;

            std::cout << this->currentCheckTaskNoBound->getFormula() << std::endl;
            std::cout << modelCopy.getTransitionMatrix() << std::endl;

            this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMax);
            auto derivativeMax = liftingModelChecker->getBound(env, region, OptimizationDirection::Maximize)->template asExplicitQuantitativeCheckResult<ConstantType>();
            ConstantType derMaxAtInit = (derivativeMax)[model.getInitialStates()[0]];
            std::cout << "derivative max: " << derMaxAtInit << std::endl;

            this->liftingModelChecker->specify(env, std::make_shared<storm::models::sparse::Dtmc<FunctionType>>(modelCopy), *checkTaskMin);
            auto derivativeMin = liftingModelChecker->getBound(env, region, OptimizationDirection::Minimize)->template asExplicitQuantitativeCheckResult<ConstantType>();

            ConstantType derMinAtInit = (derivativeMin)[model.getInitialStates()[0]];

            std::cout << "derivative min: " << derMinAtInit << std::endl;
        }

        template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
        template class DerivativeBoundFinder<RationalFunction, double>;
		}
}
