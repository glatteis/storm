#include "DerivativeBoundFinder.h"
#include "storm-pars/modelchecker/instantiation/SparseDtmcInstantiationModelChecker.h"
#include "storm-pars/utility/parametric.h"

namespace storm {
    namespace derivative {

        template<typename FunctionType>
        using VariableType = typename utility::parametric::VariableType<FunctionType>::type;
        template<typename FunctionType>
        using CoefficientType = typename utility::parametric::CoefficientType<FunctionType>::type;
            
        template<typename FunctionType, typename ConstantType>
        void DerivativeBoundFinder<FunctionType, ConstantType>::liftingTest() {
            models::sparse::Dtmc<FunctionType> modelCopy = model;
            std::cout << "Reward models:" << std::endl;
            std::cout << model.getNumberOfRewardModels() << std::endl;
            for (auto const& rewardModel : model.getRewardModels()) {
                std::cout << rewardModel.first << std::endl;
                std::cout << rewardModel.second << std::endl;
                std::cout << typeid(rewardModel.second).name() << std::endl;
            }
            /* std::vector<FunctionType> stateRewardVector; */
            /* for (uint_fast64_t i = 0; i < model.getNumberOfStates(); i++) { */
            /*     /1* stateRewardVector.push_back(utility::convertNumber<FunctionType>(-5.0)); *1/ */
            /*     if (modelCopy.getUniqueRewardModel().getStateActionReward(i) != utility::convertNumber<FunctionType>(0.0)) { */
            /*         modelCopy.getUniqueRewardModel().setStateActionReward(i, utility::convertNumber<FunctionType>(-5.0)); */
            /*     } */
            /* } */
            /* modelCopy.getUniqueRewardModel().setStateActionReward(0, utility::convertNumber<FunctionType>(-5.0)); */
            /* storm::models::sparse::StandardRewardModel<FunctionType> rewardModel(stateRewardVector); */
            /* rewardModel.set */
            /* modelCopy.addRewardModel(modelCopy.getUniqueRewardModelName(), rewardModel); */

            std::map<VariableType<FunctionType>, CoefficientType<FunctionType>> valuation;
            for (VariableType<FunctionType> parameter : storm::models::sparse::getAllParameters(modelCopy)) {
                valuation[parameter] = storm::utility::convertNumber<CoefficientType<FunctionType>>(0.5);
            }

            modelchecker::SparseDtmcInstantiationModelChecker<models::sparse::Dtmc<FunctionType>, ConstantType> modelChecker(modelCopy);
            modelChecker.specifyFormula(*currentFormula);
            modelChecker.check(Environment(), valuation);
        }

        template class DerivativeBoundFinder<RationalFunction, RationalNumber>;
        template class DerivativeBoundFinder<RationalFunction, double>;
		}
}
