#ifndef STORM_DERIVATIVEBOUNDFINDER_H
#define STORM_DERIVATIVEBOUNDFINDER_H

#include "environment/Environment.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "storm-pars/utility/parametric.h"
namespace storm {
    namespace derivative {

        template<typename FunctionType, typename ConstantType>
        class DerivativeBoundFinder {
            public:
                DerivativeBoundFinder<FunctionType, ConstantType>(
                        storm::models::sparse::Dtmc<FunctionType> const model
                ) : model(model) {
                }

            void specifyFormula(Environment const& env, modelchecker::CheckTask<logic::Formula, FunctionType> const& checkTask) {
                this->currentFormula = checkTask.getFormula().asSharedPointer();
                this->currentCheckTask = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(checkTask.substituteFormula(*currentFormula).template convertValueType<FunctionType>());

                if (!checkTask.getFormula().isRewardOperatorFormula()) {
                    this->currentFormulaNoBound = std::make_shared<storm::logic::ProbabilityOperatorFormula>(
                        checkTask.getFormula().asProbabilityOperatorFormula().getSubformula().asSharedPointer(), storm::logic::OperatorInformation(boost::none, boost::none));
                } else {
                    // No worries, this works as intended, the API is just weird.
                    this->currentFormulaNoBound = std::make_shared<storm::logic::RewardOperatorFormula>(
                        checkTask.getFormula().asRewardOperatorFormula().getSubformula().asSharedPointer());
                }
                this->currentCheckTaskNoBound = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, FunctionType>>(*currentFormulaNoBound);
                this->currentCheckTaskNoBoundConstantType = std::make_unique<storm::modelchecker::CheckTask<storm::logic::Formula, ConstantType>>(*currentFormulaNoBound);

                this->parameters = storm::models::sparse::getProbabilityParameters(model);
                if (checkTask.getFormula().isRewardOperatorFormula()) {
                    for (auto const& rewardParameter : storm::models::sparse::getRewardParameters(model)) {
                        this->parameters.insert(rewardParameter);
                    }
                }
            }

            void liftingTest();

            private:
                std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, FunctionType>> currentCheckTask;
                std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, FunctionType>> currentCheckTaskNoBound;
                std::unique_ptr<modelchecker::CheckTask<storm::logic::Formula, ConstantType>> currentCheckTaskNoBoundConstantType;
                std::shared_ptr<storm::logic::Formula const> currentFormula;
                std::shared_ptr<storm::logic::Formula const> currentFormulaNoBound;
                const models::sparse::Dtmc<FunctionType> model;
                std::set<typename utility::parametric::VariableType<FunctionType>::type> parameters;
        };

    }
}

#endif
