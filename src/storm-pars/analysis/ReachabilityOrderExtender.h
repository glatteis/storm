#pragma once
#include "storm-pars/analysis/OrderExtender.h"

namespace storm {
    namespace analysis {
        template<typename ValueType, typename ConstantType>
        class ReachabilityOrderExtender : public OrderExtender<ValueType, ConstantType> {

           public:
            typedef typename utility::parametric::VariableType<ValueType>::type VariableType;

            // Used to call the constructor of OrderExtender
            ReachabilityOrderExtender(std::shared_ptr<models::sparse::Model<ValueType>> model, std::shared_ptr<logic::Formula const> formula, bool useAssumptions = true);

            // Used to call the constructor of OrderExtender
            ReachabilityOrderExtender(storm::storage::BitVector* topStates,  storm::storage::BitVector* bottomStates, storm::storage::SparseMatrix<ValueType> matrix, bool useAssumptions = true);

           protected:
                // Override methods from OrderExtender
                void handleOneSuccessor(std::shared_ptr<Order> order, uint_fast64_t currentState, uint_fast64_t successor) override;
                std::shared_ptr<Order> getInitialOrder() override;
                virtual void addInitialStatesMinMax(std::shared_ptr<Order> order) = 0;

                std::pair<uint_fast64_t, uint_fast64_t> extendByBackwardReasoning(std::shared_ptr<Order> order, storm::storage::ParameterRegion<ValueType> region, uint_fast64_t currentState) override;
                std::pair<uint_fast64_t, uint_fast64_t> extendByForwardReasoning(std::shared_ptr<Order> order, storm::storage::ParameterRegion<ValueType> region, uint_fast64_t currentState) override;

        };
    }
}