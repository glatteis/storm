#include <storage/StronglyConnectedComponentDecomposition.h>
#include "storm-pars/analysis/ReachabilityOrderExtender.h"


namespace storm {
    namespace analysis {

        template<typename ValueType, typename ConstantType>
        ReachabilityOrderExtender<ValueType, ConstantType>::ReachabilityOrderExtender(std::shared_ptr<models::sparse::Model<ValueType>> model, std::shared_ptr<logic::Formula const> formula) : OrderExtender<ValueType, ConstantType>(model, formula) {
            // intentionally left empty
        }

        template<typename ValueType, typename ConstantType>
        ReachabilityOrderExtender<ValueType, ConstantType>::ReachabilityOrderExtender(storm::storage::BitVector* topStates,  storm::storage::BitVector* bottomStates, storm::storage::SparseMatrix<ValueType> matrix) : OrderExtender<ValueType, ConstantType>(topStates, bottomStates, matrix) {
            // intentionally left empty
        }

        template <typename ValueType, typename ConstantType>
        void ReachabilityOrderExtender<ValueType, ConstantType>::checkParOnStateMonRes(uint_fast64_t s, std::shared_ptr<Order> order, typename OrderExtender<ValueType, ConstantType>::VariableType param, std::shared_ptr<MonotonicityResult<VariableType>> monResult) {
            monResult->updateMonotonicityResult(param, this->monotonicityChecker.checkLocalMonotonicity(order, s, param, this->region));
        }

        template <typename ValueType, typename ConstantType>
        std::tuple<std::shared_ptr<Order>, uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::toOrder(storage::ParameterRegion<ValueType> region, std::shared_ptr<MonotonicityResult<VariableType>> monRes) {
            return this->extendOrder(nullptr, region, monRes, nullptr);
        }

        template<typename ValueType, typename ConstantType>
        std::tuple<std::shared_ptr<Order>, uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendOrder(std::shared_ptr<Order> order, storm::storage::ParameterRegion<ValueType> region, std::shared_ptr<MonotonicityResult<VariableType>> monRes, std::shared_ptr<expressions::BinaryRelationExpression> assumption) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This function is only implemented in sub-classes");
        }

        template<typename ValueType, typename ConstantType>
        std::tuple<std::shared_ptr<Order>, uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendOrder(std::shared_ptr<Order> order, std::shared_ptr<MonotonicityResult<VariableType>> monRes, std::shared_ptr<expressions::BinaryRelationExpression> assumption) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This function is only implemented in sub-classes");
        }

        template<typename ValueType, typename ConstantType>
        std::pair<uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendByBackwardReasoning(std::shared_ptr<Order> order, uint_fast64_t currentState) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This function is only implemented in sub-classes");
        }

        template<typename ValueType, typename ConstantType>
        std::pair<uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendNormal(std::shared_ptr<Order> order, uint_fast64_t currentState, const std::vector<uint_fast64_t> &successors, bool allowMerge)  {
            // when it is cyclic and the current state is part of an SCC we do forwardreasoning
            if (this->cyclic && !order->isTrivial(currentState) && order->contains(currentState)) {
                // Try to extend the order for this scc
                return  extendByForwardReasoning(order, currentState, successors, allowMerge);
            } else {
                assert (order->isTrivial(currentState) || !order->contains(currentState));
                // Do backward reasoning, all successor states must be in the order
                return  extendByBackwardReasoning(order, currentState, successors, allowMerge);
            }
        }

        template <typename ValueType, typename ConstantType>
        std::pair<uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendByBackwardReasoning(std::shared_ptr<Order> order, uint_fast64_t currentState, std::vector<uint_fast64_t> const& successors, bool allowMerge) {
            assert (!order->isOnlyBottomTopOrder());
            assert (successors.size() > 1);

            bool pla = (this->usePLA.find(order) != this->usePLA.end() && this->usePLA.at(order));
            std::vector<uint_fast64_t> sortedSuccs;

            if (pla && (this->continueExtending.find(order) == this->continueExtending.end() || this->continueExtending.at(order))) {
                for (auto& state1 : successors) {
                    if (sortedSuccs.size() == 0) {
                        sortedSuccs.push_back(state1);
                    } else {
                        bool added = false;
                        for (auto itr = sortedSuccs.begin(); itr != sortedSuccs.end(); ++itr) {
                            auto& state2 = *itr;
                            auto compareRes = order->compareFast(state1, state2);
                            if (compareRes == Order::NodeComparison::UNKNOWN) {
                                compareRes = this->addStatesBasedOnMinMax(order, state1, state2);
                            }
                            if (compareRes == Order::NodeComparison::UNKNOWN) {
                                // If fast comparison did not help, we continue by checking "slow" comparison
                                compareRes = order->compare(state1, state2);
                            }
                            if (compareRes == Order::NodeComparison::ABOVE ||
                                compareRes == Order::NodeComparison::SAME) {
                                // insert at current pointer (while keeping other values)
                                sortedSuccs.insert(itr, state1);
                                added = true;
                                break;
                            } else if (compareRes == Order::NodeComparison::UNKNOWN) {
                                this->continueExtending[order] = false;
                                return {state1, state2};
                            }
                        }
                        if (!added) {
                            sortedSuccs.push_back(state1);
                        }
                    }
                }
            } else {
                auto temp = order->sortStatesUnorderedPair(&successors);
                if (temp.first.first != this->numberOfStates) {
                    return temp.first;
                } else {
                    sortedSuccs = std::move(temp.second);
                }
            }

            if (order->compare(sortedSuccs[0], sortedSuccs[sortedSuccs.size() - 1]) == Order::SAME) {
                if (!order->contains(currentState)) {
                    order->addToNode(currentState, order->getNode(sortedSuccs[0]));
                } else {
                    order->merge(currentState, sortedSuccs[0]);
                }
            } else {
                if (!order->contains(sortedSuccs[0])) {
                    assert (order->isBottomState(sortedSuccs[sortedSuccs.size() - 1]));
                    assert (sortedSuccs.size() == 2);
                    order->addAbove(sortedSuccs[0], order->getBottom());
                }
                if (!order->contains(sortedSuccs[sortedSuccs.size() - 1])) {
                    assert (order->isTopState(sortedSuccs[0]));
                    assert (sortedSuccs.size() == 2);
                    order->addBelow(sortedSuccs[sortedSuccs.size() - 1], order->getTop());
                }
                // sortedSuccs[0] is highest
                if (!order->contains(currentState)) {
                    order->addBetween(currentState, sortedSuccs[0], sortedSuccs[sortedSuccs.size()-1]);
                } else {
                    order->addRelation(sortedSuccs[0], currentState, allowMerge);
                    order->addRelation(currentState, sortedSuccs[sortedSuccs.size() - 1], allowMerge);
                }

            }
            assert (order->contains(currentState) && order->compare(order->getNode(currentState), order->getBottom()) == Order::ABOVE && order->compare(order->getNode(currentState), order->getTop()) == Order::BELOW);
            return {this->numberOfStates, this->numberOfStates};
        }

        template <typename ValueType, typename ConstantType>
        std::pair<uint_fast64_t, uint_fast64_t> ReachabilityOrderExtender<ValueType, ConstantType>::extendByForwardReasoning(std::shared_ptr<Order> order, uint_fast64_t currentState, std::vector<uint_fast64_t> const& successors, bool allowMerge)  {
            assert (successors.size() > 1);
            assert (order->contains(currentState));
            assert (this->cyclic);

            std::vector<uint_fast64_t> statesSorted;
            statesSorted.push_back(currentState);
            bool pla = (this->usePLA.find(order) != this->usePLA.end() && this->usePLA.at(order));
            // Go over all states
            bool oneUnknown = false;
            bool unknown = false;
            uint_fast64_t s1 = this->numberOfStates;
            uint_fast64_t s2 = this->numberOfStates;
            for (auto& state : successors) {
                unknown = false;
                bool added = false;
                for (auto itr = statesSorted.begin(); itr != statesSorted.end(); ++itr) {
                    auto compareRes = order->compareFast(state, (*itr));
                    if (pla && compareRes == Order::NodeComparison::UNKNOWN) {
                        compareRes = this->addStatesBasedOnMinMax(order, state, (*itr));
                    }
                    if (compareRes == Order::NodeComparison::UNKNOWN) {
                        compareRes = order->compare(state, *itr);
                    }
                    if (compareRes == Order::NodeComparison::ABOVE || compareRes == Order::NodeComparison::SAME) {
                        if (!order->contains(state) && compareRes == Order::NodeComparison::ABOVE) {
                            order->add(state);
                            order->addStateToHandle(state);
                        }
                        added = true;
                        // insert at current pointer (while keeping other values)
                        statesSorted.insert(itr, state);
                        break;
                    } else if (compareRes == Order::NodeComparison::UNKNOWN && !oneUnknown) {
                        // We miss state in the result.
                        s1 = state;
                        s2 = *itr;
                        oneUnknown = true;
                        added = true;
                        break;
                    } else if (compareRes == Order::NodeComparison::UNKNOWN && oneUnknown) {
                        unknown = true;
                        added = true;
                        break;
                    }
                }
                if (!(unknown && oneUnknown) && !added ) {
                    // State will be last in the list
                    statesSorted.push_back(state);
                }
                if (unknown && oneUnknown) {
                    break;
                }
            }
            if (!unknown && oneUnknown) {
                assert (statesSorted.size() == successors.size());
                s2 = this->numberOfStates;
            }

            if (s1 == this->numberOfStates) {
                assert (statesSorted.size() == successors.size() + 1);
                // all could be sorted, no need to do anything
            } else if (s2 == this->numberOfStates) {
                if (!order->contains(s1)) {
                    order->add(s1);
                }

                if (statesSorted[0] == currentState) {
                    order->addRelation(s1, statesSorted[0], allowMerge);
                    auto res = order->compare(s1, statesSorted[0]);
                    assert ((order->compare(s1, statesSorted[0]) == Order::ABOVE) || (allowMerge && (order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::SAME)));
                    order->addRelation(s1, statesSorted[statesSorted.size() - 1], allowMerge);
                    assert ((order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::ABOVE) || (allowMerge && (order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::SAME)));
                    order->addStateToHandle(s1);
                } else if (statesSorted[statesSorted.size() - 1] == currentState) {
                    order->addRelation(statesSorted[0], s1, allowMerge);
                    assert ((order->compare(s1, statesSorted[0]) == Order::BELOW) || (allowMerge && (order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::SAME)));
                    order->addRelation(statesSorted[statesSorted.size() - 1], s1, allowMerge);
                    assert ((order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::BELOW) || (allowMerge && (order->compare(s1, statesSorted[statesSorted.size() - 1]) == Order::SAME)));
                    order->addStateToHandle(s1);
                } else {
                    bool continueSearch = true;
                    for (auto& entry :  this->matrix.getRow(currentState)) {
                        if (entry.getColumn() == s1) {
                            if (entry.getValue().isConstant()) {
                                continueSearch = false;
                            }
                        }
                    }
                    if (continueSearch) {
                        for (auto& i : statesSorted) {
                            if (order->compare(i, s1) == Order::UNKNOWN) {
                                return {i, s1};
                            }
                        }
                    }
                }
            } else {
                return {s1, s2};
            }
            assert (order->contains(currentState) && order->compare(order->getNode(currentState), order->getBottom()) == Order::ABOVE && order->compare(order->getNode(currentState), order->getTop()) == Order::BELOW);
            return {this->numberOfStates, this->numberOfStates};
        }

        template<typename ValueType, typename ConstantType>
        bool ReachabilityOrderExtender<ValueType, ConstantType>::extendByAssumption(std::shared_ptr<Order> order, uint_fast64_t stateSucc1, uint_fast64_t stateSucc2) {
            bool usePLANow = this->usePLA.find(order) != this->usePLA.end() && this->usePLA[order];
            assert (order->compare(stateSucc1, stateSucc2) == Order::UNKNOWN);
            auto assumptions = usePLANow ? this->assumptionMaker->createAndCheckAssumptions(stateSucc1, stateSucc2,  order, this->region, this->minValues[order], this->maxValues[order]) : this->assumptionMaker->createAndCheckAssumptions(stateSucc1, stateSucc2, order, this->region);
            if (assumptions.size() == 1 && assumptions.begin()->second == AssumptionStatus::VALID) {
                this->handleAssumption(order, assumptions.begin()->first);
                // Assumptions worked, we continue
                return true;
            }
            return false;
        }

        template<typename ValueType, typename ConstantType>
        void ReachabilityOrderExtender<ValueType, ConstantType>::handleOneSuccessor(std::shared_ptr<Order> order, uint_fast64_t currentState, uint_fast64_t successor) {
            assert (order->contains(successor));
            if (currentState != successor) {
                if (order->contains(currentState)) {
                    order->merge(currentState, successor);
                } else {
                    order->addToNode(currentState, order->getNode(successor));
                }
            }
        }

        template <typename ValueType, typename ConstantType>
        std::shared_ptr<Order> ReachabilityOrderExtender<ValueType, ConstantType>::getBottomTopOrder() {
            if (this->bottomTopOrder == nullptr) {
                assert (this->model != nullptr);
                STORM_LOG_THROW(this->matrix.getRowCount() == this->matrix.getColumnCount(), exceptions::NotSupportedException,"Creating order not supported for non-square matrix");
                modelchecker::SparsePropositionalModelChecker<models::sparse::Model<ValueType>> propositionalChecker(*(this->model));
                storage::BitVector phiStates;
                storage::BitVector psiStates;
                assert (this->formula->isProbabilityOperatorFormula());
                if (this->formula->asProbabilityOperatorFormula().getSubformula().isUntilFormula()) {
                    phiStates = propositionalChecker.check(
                            this->formula->asProbabilityOperatorFormula().getSubformula().asUntilFormula().getLeftSubformula())->asExplicitQualitativeCheckResult().getTruthValuesVector();
                    psiStates = propositionalChecker.check(
                            this->formula->asProbabilityOperatorFormula().getSubformula().asUntilFormula().getRightSubformula())->asExplicitQualitativeCheckResult().getTruthValuesVector();
                } else {
                    assert (this->formula->asProbabilityOperatorFormula().getSubformula().isEventuallyFormula());
                    phiStates = storage::BitVector(this->numberOfStates, true);
                    psiStates = propositionalChecker.check(
                            this->formula->asProbabilityOperatorFormula().getSubformula().asEventuallyFormula().getSubformula())->asExplicitQualitativeCheckResult().getTruthValuesVector();
                }
                // Get the maybeStates
                std::pair<storage::BitVector, storage::BitVector> statesWithProbability01 = utility::graph::performProb01(this->model->getBackwardTransitions(), phiStates, psiStates);
                storage::BitVector topStates = statesWithProbability01.second;
                storage::BitVector bottomStates = statesWithProbability01.first;

                STORM_LOG_THROW(topStates.begin() != topStates.end(), exceptions::NotSupportedException,"Formula yields to no 1 states");
                STORM_LOG_THROW(bottomStates.begin() != bottomStates.end(), exceptions::NotSupportedException,"Formula yields to no zero states");
                auto& matrix = this->model->getTransitionMatrix();
                std::vector<uint64_t> firstStates;

                storm::storage::BitVector subStates (topStates.size(), true);
                for (auto state : topStates) {
                    firstStates.push_back(state);
                    subStates.set(state, false);
                }
                for (auto state : bottomStates) {
                    firstStates.push_back(state);
                    subStates.set(state, false);
                }
                this->cyclic = storm::utility::graph::hasCycle(matrix, subStates);
                storm::storage::StronglyConnectedComponentDecomposition<ValueType> decomposition;
                if (this->cyclic) {
                    storm::storage::StronglyConnectedComponentDecompositionOptions options;
                    options.forceTopologicalSort();
                    decomposition = storm::storage::StronglyConnectedComponentDecomposition<ValueType>(matrix, options);
                }
                auto statesSorted = storm::utility::graph::getTopologicalSort(matrix.transpose(), firstStates);
                this->bottomTopOrder = std::shared_ptr<Order>(new Order(&topStates, &bottomStates, this->numberOfStates, std::move(decomposition), std::move(statesSorted)));

                // Build stateMap
                this->buildStateMap(subStates);
            }

            if (this->minValuesInit && this->maxValuesInit) {
                this->continueExtending[this->bottomTopOrder] = true;
                this->usePLA[this->bottomTopOrder] = true;
                this->minValues[this->bottomTopOrder] = std::move(this->minValuesInit.get());
                this->maxValues[this->bottomTopOrder] = std::move(this->maxValuesInit.get());
            } else {
                this->usePLA[this->bottomTopOrder] = false;
            }
            return this->bottomTopOrder;
        }

        template class ReachabilityOrderExtender<RationalFunction, double>;
        template class ReachabilityOrderExtender<RationalFunction, RationalNumber>;

    }
}