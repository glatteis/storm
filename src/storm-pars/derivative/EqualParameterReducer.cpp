#include "EqualParameterReducer.h"
#include <carl/core/Variable.h>
#include <carl/core/VariablePool.h>
#include <cstdint>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include "adapters/RationalFunctionAdapter.h"
#include "adapters/RationalNumberAdapter.h"
#include "logic/UntilFormula.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "models/sparse/StandardRewardModel.h"
#include "models/sparse/StateLabeling.h"
#include "storage/BitVector.h"
#include "storage/FlexibleSparseMatrix.h"
#include "storage/SparseMatrix.h"
#include "storm-pars/transformer/SparseParametricDtmcSimplifier.h"
#include "utility/constants.h"
#include "utility/graph.h"
#include "utility/macros.h"

namespace storm {
namespace derivative {

models::sparse::Dtmc<RationalFunction> EqualParameterReducer::minimizeEqualParameters(
    models::sparse::Dtmc<RationalFunction> dtmc, modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask) {
    // dtmc.writeDotToStream(std::cout);
    storage::SparseMatrix<RationalFunction> transitionMatrix = dtmc.getTransitionMatrix();

    STORM_LOG_ASSERT(dtmc.getTransitionMatrix().isProbabilistic(), "Matrix not probabilistic!");
    // Repeat this algorithm until we can't mimimize anything anymore
    bool somethingChanged = true;
    auto allParameters = storm::models::sparse::getAllParameters(dtmc);

    models::sparse::StateLabeling newLabels(dtmc.getStateLabeling());

    // Check the reward model - do not touch states with rewards
    boost::optional<std::vector<RationalFunction>> stateRewardVector;
    boost::optional<std::string> stateRewardName;
    if (checkTask.getFormula().isRewardOperatorFormula()) {
        if (checkTask.isRewardModelSet()) {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel(checkTask.getRewardModel()).getStateRewardVector();
            stateRewardName = checkTask.getRewardModel();
        } else {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel("").getStateRewardVector();
            stateRewardName = dtmc.getUniqueRewardModelName();
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
        std::map<storm::RationalFunctionVariable, RationalFunction> pFunctions;

        uint_fast64_t numberOfSearchingTransitions = 0;
        // for (auto const& parameter : allParameters) {
        //     alreadyVisitedStates[parameter] = std::map<uint_fast64_t, std::vector<uint_fast64_t>>();
        //     talliedUpProbabilities[parameter] = std::map<uint_fast64_t, std::vector<RationalNumber>>();
        //     increasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
        //     decreasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
        //     joiningLabels[parameter] = std::map<uint_fast64_t, std::set<std::string>>();
        // }
        // Search for all occurences of parameters and fill up the maps
        for (uint_fast64_t row = 0; row < transitionMatrix.getRowCount(); row++) {
            for (auto const& entry : transitionMatrix.getRow(row)) {
                if (entry.getValue().isConstant()) {
                    continue;
                }
                STORM_LOG_ASSERT(entry.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
                auto parameter = *entry.getValue().gatherVariables().begin();
                auto cache = entry.getValue().nominatorAsPolynomial().pCache();
                STORM_LOG_ASSERT(entry.getValue().denominator().isOne() && entry.getValue().nominator().isUnivariate() &&
                                     entry.getValue().nominator().getSingleVariable() == parameter && entry.getValue().nominator().factorization().size() == 1,
                                 "Flip minimization only supports simple pMCs.");
                pFunctions[parameter] = entry.getValue();

                if (utility::isOne(entry.getValue().derivative(entry.getValue().nominator().getSingleVariable()))) {
                    numberOfSearchingTransitions++;
                    auto parameter = entry.getValue().nominatorAsPolynomial().getSingleVariable();
                    STORM_LOG_ASSERT(alreadyVisitedStates[parameter].count(row) == 0, "Flip minimization only supports simple pMCs.");
                    alreadyVisitedStates[parameter][row].push_back(row);
                    talliedUpProbabilities[parameter][row].push_back(utility::one<RationalNumber>());
                    increasingSuccessors[parameter][row] = entry.getColumn();
                    doneSearching[parameter][row] = false;
                    joiningLabels[parameter][row] = newLabels.getLabelsOfState(row);
                } else if (utility::isOne(-entry.getValue().derivative(entry.getValue().nominator().getSingleVariable()))) {
                    decreasingSuccessors[parameter][row] = entry.getColumn();
                } else {
                    STORM_LOG_ASSERT(false, "Flip minimization only supports transitions with a single parameter.");
                }
            }
        }
        // Step by step, do a backwards search
        uint_fast64_t transitionsDoneSearching = 0;
        auto backwardsTransitions = transitionMatrix.transpose(true);
        // Final results of the search: Sets of sets of states that we can transform together
        // The maps are <parameter> -> <state where the paths join> -> <set of states that leads there, with this probability>
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
                    if (doneSearching.at(parameter).at(firstState)) {
                        continue;
                    }
                    auto row = backwardsTransitions.getRow(visitedStates.back());

                    storage::MatrixEntry<uint_fast64_t, RationalFunction> entry;
                    bool entryFound = false;
                    for (auto const& loopEntry : row) {
                        if (loopEntry.getValue().isZero()) {
                            continue;
                        }
                        // We only support one backlink per state for now, so if we already found an entry, abort
                        if (entryFound) {
                            std::cout << "Abort because of two backlinks" << std::endl;
                            entryFound = false;
                            break;
                        }
                        entry = loopEntry;
                        entryFound = true;
                    }
                    bool doStep = true;
                    if (!entryFound || !entry.getValue().isConstant()) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                        doStep = false;
                    }

                    uint_fast64_t currentState = visitedStates.back();
                    if (doStep) {
                        // If the last label is different, we stop at the next state (and still join if possible)
                        // If the state reward vector is not zero, we stop at the next state (and still join if possible)
                        if (joiningLabels.at(parameter).at(firstState) != newLabels.getLabelsOfState(entry.getColumn()) ||
                            (stateRewardVector && !(*stateRewardVector)[visitedStates.back()].isZero())) {
                            transitionsDoneSearching++;
                            doneSearching[parameter][firstState] = true;
                        }
                        auto constantProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                        currentState = entry.getColumn();

                        // The next probability is the previous one times the considered tranistion
                        RationalNumber numberToPush = talliedUpProbabilities[parameter][firstState].back() * constantProbability;
                        // std::cout << "pushing " << numberToPush << std::endl;
                        talliedUpProbabilities[parameter][firstState].push_back(numberToPush);
                        alreadyVisitedStates[parameter][firstState].push_back(currentState);
                    }

                    // If we are not done searching, don't join anything
                    if (!doneSearching[parameter][firstState]) {
                        continue;
                    }

                    // Check if any search of the same parameter has already gotten here
                    for (auto const& pair3 : alreadyVisitedStates[parameter]) {
                        auto otherFirstState = pair3.first;
                        // Don't join a state with itself
                        if (otherFirstState == firstState) {
                            continue;
                        }
                        // If the labels don't match, this is not a match
                        if (joiningLabels.at(parameter).at(otherFirstState) != joiningLabels.at(parameter).at(firstState)) {
                            continue;
                        }
                        auto otherVisitedStates = pair3.second;

                        auto state = otherVisitedStates.back();

                        // doneSearching[parameter][otherFirstState] can only be false if the other path has gotten here,
                        // but has not noticed it is the end of the path yet. The easiest way to mitigate this is let
                        // that path join us when it notices.
                        if (state == currentState && doneSearching.at(parameter).at(otherFirstState)) {
                            // Wow!!! Another search of the same parameter has gotten here.
                            if (!joinedFirstStates.at(parameter).count(state)) {
                                joinedFirstStates[parameter][state] = std::map<uint_fast64_t, RationalNumber>();
                            }
                            RationalNumber p1 = talliedUpProbabilities.at(parameter).at(firstState).back();
                            joinedFirstStates[parameter][state][firstState] = p1;
                            RationalNumber p2 = talliedUpProbabilities.at(parameter).at(otherFirstState).back();
                            joinedFirstStates[parameter][state][otherFirstState] = p2;

                            break;
                        }
                    }
                }
            }
        }

        // Check if we found any parameters that we can join into one transition, and if so, do it
        for (auto const& parameter : allParameters) {
            if (joinedFirstStates.at(parameter).size() > 0) {
                somethingChanged = true;
                break;
            }
        }
        // States that are on parametric paths
        std::map<storm::RationalFunctionVariable, std::set<uint_fast64_t>> statesOnParametricPaths;
        for (auto const& parameter : allParameters) {
            for (auto const& joinedFirstStatesPair : joinedFirstStates.at(parameter)) {
                for (auto const& firstStatePair : joinedFirstStatesPair.second) {
                    auto firstState = firstStatePair.first;
                    for (auto const& visitedState : alreadyVisitedStates.at(parameter).at(firstState)) {
                        statesOnParametricPaths[parameter].emplace(visitedState);
                    }
                }
            }
        }
        if (!somethingChanged) {
            break;
        }

        storage::SparseMatrix<RationalFunction> wipMatrix(transitionMatrix);
        for (auto const& parameter : allParameters) {
            // We add three new states to the DTMC for this procedure for every state where parameters join
            if (joinedFirstStates.at(parameter).size() == 0) {
                continue;
            }
            std::cout << "joining " << joinedFirstStates.at(parameter).size() << " cases for " << parameter << std::endl;
            uint_fast64_t runningCounter = wipMatrix.getRowCount();

            std::map<uint_fast64_t, std::vector<uint_fast64_t>> newStates;
            for (auto const& pair : joinedFirstStates.at(parameter)) {
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

            storage::SparseMatrixBuilder<RationalFunction> builder;
            for (uint_fast64_t row = 0; row < wipMatrix.getRowCount(); row++) {
                // If from row, we are starting with a transition we will flip...
                if (joinedFirstStates.at(parameter).count(row)) {
                    std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>,
                                        std::greater<std::pair<uint_fast64_t, RationalFunction>>>
                        insertInOrder;
                    // For checking correctness
                    RationalFunction rowSum = utility::zero<RationalFunction>();
                    // We will join until these transitions
                    auto statesBeforeParams = joinedFirstStates.at(parameter).at(row);

                    std::set<uint_fast64_t> alreadyConsidered;
                    // We go, in reverse, through the visited states of every state
                    // We choose all constant-only transitions in the process to add to the first state in the end
                    // For example
                    // S1 -> 15/16 -> S2
                    // S2 -> 1/15 -> S3 -> .. (non-parametric)
                    // S2 -> 1/15 -> S4 -> .. (non-parametric)
                    // S2 -> 13/15 -> S5 -> p
                    // S1 -> 1/16 -> S6 -> p
                    // we find all the non-parametric transitions and add them to S1
                    // S1 -> 15/16 * 1/15 -> S3
                    // S1 -> 15/16 * 1/15 -> S4
                    // S1 -> SN1
                    // ...
                    for (auto const& startingState : statesBeforeParams) {
                        auto visitedStates = alreadyVisitedStates.at(parameter).at(startingState.first);
                        // The probability to reach the starting state. It will be divided as we step back
                        RationalNumber probability = talliedUpProbabilities.at(parameter).at(startingState.first).back();
                        for (uint_fast64_t i = 1; i < visitedStates.size(); i++) {
                            // Previous (next in reality) visited state is i-1. Divide the probability
                            auto const& visitedState = visitedStates.at(i);
                            if (alreadyConsidered.count(visitedState)) {
                                break;
                            }
                            alreadyConsidered.emplace(visitedState);
                            auto const& nextVisitedState = visitedStates[i - 1];
                            // get the transition probability between visitedState and nextVisitedState
                            boost::optional<RationalNumber> transitionProbability;
                            for (auto const& entry : wipMatrix.getRow(visitedState)) {
                                if (entry.getColumn() == nextVisitedState) {
                                    transitionProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                                    break;
                                }
                            }
                            STORM_LOG_ASSERT(transitionProbability, "Very bad internal error: could not trace back found path");
                            // std::cout << probability << std::endl;
                            // std::cout << *transitionProbability << std::endl;
                            probability = probability / *transitionProbability;
                            // Go through the matrix entries to see what other stuff is branching here
                            for (auto const& entry : wipMatrix.getRow(visitedState)) {
                                auto state = entry.getColumn();
                                if (!statesOnParametricPaths.at(parameter).count(state)) {
                                    // This state is not on a parametric path, thus it needs to be added to the original state
                                    RationalFunction value = utility::convertNumber<RationalFunction>(probability) * entry.getValue();
                                    insertInOrder.emplace(std::make_pair(state, value));
                                    rowSum += value;
                                }
                            }
                        }
                    }

                    RationalNumber summedProbability = 0;
                    for (auto const& entry2 : statesBeforeParams) {
                        summedProbability += entry2.second;
                    }
                    // Add transition to the state that will have the parametric decision
                    insertInOrder.emplace(std::make_pair(newStates.at(row)[0], utility::convertNumber<RationalFunction>(summedProbability)));
                    rowSum += utility::convertNumber<RationalFunction>(summedProbability);

                    STORM_LOG_ASSERT(utility::isOne(rowSum), "Internal error: Row sum of row " << row << " is " << rowSum << ", not one!");
                    while (!insertInOrder.empty()) {
                        auto pair = insertInOrder.top();
                        insertInOrder.pop();
                        builder.addNextValue(row, pair.first, pair.second);
                    }
                } else {  // Else, just add the stuff back into the matrix
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        builder.addNextValue(row, entry.getColumn(), entry.getValue());
                    }
                }
            }

            // RationalFunction parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
            RationalFunction parameterAsFunction = pFunctions.at(parameter);

            for (auto const& pair : joinedFirstStates.at(parameter)) {
                auto row = pair.first;
                RationalNumber summedProbability = 0;
                for (auto const& entry2 : joinedFirstStates.at(parameter).at(row)) {
                    summedProbability += entry2.second;
                }
                // Add parametric transitions
                builder.addNextValue(newStates.at(row)[0], newStates.at(row)[1], parameterAsFunction);
                builder.addNextValue(newStates.at(row)[0], newStates.at(row)[2], utility::one<RationalFunction>() - parameterAsFunction);

                // Add transitions from p and 1-p states to actual successors
                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>,
                                    std::greater<std::pair<uint_fast64_t, RationalFunction>>>
                    insertInOrder1;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto increasingSuccessor = increasingSuccessors.at(parameter).at(id);
                    insertInOrder1.emplace(std::make_pair(increasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder1.empty()) {
                    auto pair = insertInOrder1.top();
                    insertInOrder1.pop();
                    builder.addNextValue(newStates.at(row)[1], pair.first, pair.second);
                }

                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>,
                                    std::greater<std::pair<uint_fast64_t, RationalFunction>>>
                    insertInOrder2;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto decreasingSuccessor = decreasingSuccessors.at(parameter).at(id);
                    insertInOrder2.emplace(std::make_pair(decreasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder2.empty()) {
                    auto pair = insertInOrder2.top();
                    insertInOrder2.pop();
                    builder.addNextValue(newStates.at(row)[2], pair.first, pair.second);
                }

                auto someJoinedState = joinedFirstStates.at(parameter).at(row).begin()->first;
                for (auto const& label : joiningLabels.at(parameter).at(someJoinedState)) {
                    if (!newLabels.getLabels().count(label)) {
                        newLabels.addLabel(label);
                    }
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        newLabels.addLabelToState(label, newStates.at(row)[i]);
                    }
                }
                if (stateRewardVector) {
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        stateRewardVector->push_back(storm::utility::zero<RationalFunction>());
                    }
                }
            }
            wipMatrix = builder.build();
        }

        // Remove the states that are not needed anymore
        // First we mark which states to delete, then we complement the BitVector
        storage::BitVector statesToKeep(wipMatrix.getRowCount());
        for (auto const& parameter : allParameters) {
            for (auto const& pair : joinedFirstStates.at(parameter)) {
                auto joiningState = pair.first;
                for (auto const& pair2 : pair.second) {
                    // The identifier - the state the search started from
                    auto firstState = pair2.first;
                    auto visitedStates = alreadyVisitedStates.at(parameter).at(firstState);

                    // Delete everything up to the point where the joining state, which is not deleted, is met
                    for (auto const& visitedState : visitedStates) {
                        if (visitedState == joiningState) {
                            break;
                        }
                        // std::cout << "Deleting " << visitedState << std::endl;
                        statesToKeep.set(visitedState, true);
                    }
                }
            }
        }

        // Don't forget :)
        statesToKeep.complement();
        wipMatrix = wipMatrix.getSubmatrix(false, statesToKeep, statesToKeep, true);
        models::sparse::StateLabeling copyOfNewLabels = newLabels;
        newLabels = copyOfNewLabels.getSubLabeling(statesToKeep);
        if (stateRewardVector) {
            std::vector<RationalFunction> newStateRewards;
            for (uint_fast64_t i = 0; i < stateRewardVector->size(); i++) {
                if (statesToKeep[i]) {
                    newStateRewards.push_back((*stateRewardVector)[i]);
                }
            }
            stateRewardVector = newStateRewards;
        }
        storage::SparseMatrix<RationalFunction> copyOfWipMatrix = wipMatrix;
        transitionMatrix = copyOfWipMatrix;
    }
    models::sparse::Dtmc<RationalFunction> newDTMC(transitionMatrix, newLabels);

    if (stateRewardVector) {
        models::sparse::StandardRewardModel<RationalFunction> newRewardModel(*stateRewardVector);
        newDTMC.addRewardModel(*stateRewardName, newRewardModel);
    }

    STORM_LOG_ASSERT(newDTMC.getTransitionMatrix().isProbabilistic(), "Internal error: resulting matrix not probabilistic!");

    // Use model simplifier to get rid of unnecessary states
    // storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<RationalFunction>> simplifier(newDTMC);
    // STORM_LOG_ASSERT(simplifier.simplify(checkTask.getFormula()), "Could not simplify model.");
    // newDTMC = *simplifier.getSimplifiedModel()->template as<models::sparse::Dtmc<RationalFunction>>();
    // newDTMC.writeDotToStream(std::cout);
    return newDTMC;
}

models::sparse::Dtmc<RationalFunction> EqualParameterReducer::timeTravel(models::sparse::Dtmc<RationalFunction> dtmc,
                                                                         modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask) {
    carl::freshVariable(carl::VariableType::VT_REAL);
    storage::SparseMatrix<RationalFunction> transitionMatrix = dtmc.getTransitionMatrix();
    uint_fast64_t initialState = dtmc.getInitialStates().getNextSetIndex(0);

    STORM_LOG_ASSERT(dtmc.getTransitionMatrix().isProbabilistic(), "Matrix not probabilistic!");

    auto allParameters = storm::models::sparse::getAllParameters(dtmc);

    models::sparse::StateLabeling newLabels(dtmc.getStateLabeling());

    // Check the reward model - do not touch states with rewards
    boost::optional<std::vector<RationalFunction>> stateRewardVector;
    boost::optional<std::string> stateRewardName;
    if (checkTask.getFormula().isRewardOperatorFormula()) {
        if (checkTask.isRewardModelSet()) {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel(checkTask.getRewardModel()).getStateRewardVector();
            stateRewardName = checkTask.getRewardModel();
        } else {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel("").getStateRewardVector();
            stateRewardName = dtmc.getUniqueRewardModelName();
        }
    }

    auto const constantVariable = carl::VariablePool::getInstance().getFreshPersistentVariable();
    auto topologicalOrdering = utility::graph::getTopologicalSort<RationalFunction>(transitionMatrix, {initialState});
    auto flexibleMatrix = storage::FlexibleSparseMatrix<RationalFunction>(transitionMatrix);
    std::map<uint_fast64_t, bool> alreadyVisited;
    
    std::stack<uint_fast64_t> topologicalOrderingStack;
    for (auto rit = topologicalOrdering.begin(); rit != topologicalOrdering.end(); ++rit) {
        topologicalOrderingStack.push(*rit);
    }
    
    
    // Initialize counting
    std::map<RationalFunctionVariable, std::map<uint_fast64_t, std::set<uint_fast64_t>>> treeStates;
    std::map<RationalFunctionVariable, std::set<uint_fast64_t>> workingSets;
    // TODO you don't need to count anew every time if you're smort
    // 
    // 
    // dtmc.writeDotToStream(std::cout);

    
    auto backwardsTransitions = flexibleMatrix.createSparseMatrix().transpose(true);
    // Count number of parameter occurences per state
    for (uint_fast64_t row = 0; row < flexibleMatrix.getRowCount(); row++) {
        for (auto const& entry : flexibleMatrix.getRow(row)) {
            if (entry.getValue().isConstant()) {
                continue;
            }
            STORM_LOG_ERROR_COND(entry.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
            auto parameter = *entry.getValue().gatherVariables().begin();
            auto cache = entry.getValue().nominatorAsPolynomial().pCache();
            STORM_LOG_ERROR_COND(entry.getValue().denominator().isOne() && entry.getValue().nominator().isUnivariate() &&
                                     entry.getValue().nominator().getSingleVariable() == parameter &&
                                     entry.getValue().nominator().factorization().size() == 1,
                                 "Flip minimization only supports simple pMCs.");

            STORM_LOG_ERROR_COND(flexibleMatrix.getRow(row).size() == 2,
                                 "Flip minimization only supports transitions with a single parameter.");
            workingSets[parameter].emplace(row);
            treeStates[parameter][row].emplace(row);
            if (utility::isOne(entry.getValue().derivative(entry.getValue().nominator().getSingleVariable()))) {
            } else if (utility::isOne(-entry.getValue().derivative(entry.getValue().nominator().getSingleVariable()))) {
            } else {
                STORM_LOG_ERROR_COND(false, "Flip minimization only supports transitions with a single parameter.");
            }
        }
    }

    updateTreeStates(treeStates, workingSets, flexibleMatrix, allParameters);
    while (!topologicalOrderingStack.empty()) {
        auto state = topologicalOrderingStack.top();
        topologicalOrderingStack.pop();
        // Check if we can reach more than one var from here (by the original matrix)
        bool moreThanOneVarReachable = false;
        for (auto const& parameter : allParameters) {
            if (!treeStates.at(parameter).count(state)) {
                continue;
            }
            auto const entry = treeStates.at(parameter).at(state);
            if (entry.size() >= 2) {
                moreThanOneVarReachable = true;
            }
        }
        if (!moreThanOneVarReachable) {
            continue;
        }
        jipConvert(state, flexibleMatrix, alreadyVisited, treeStates, allParameters);
        // models::sparse::Dtmc<RationalFunction> newnewDTMC(flexibleMatrix.createSparseMatrix(), newLabels);
        // newnewDTMC.writeDotToStream(std::cout);
        
        // Now our matrix is in Jip normal form. Now just re-order
        std::map<storm::RationalFunctionVariable, std::set<uint_fast64_t>> parameterBuckets;
        std::map<storm::RationalFunctionVariable, RationalFunction> cumulativeProbabilities;

        std::map<uint_fast64_t, uint_fast64_t> pTransitions;
        std::map<uint_fast64_t, uint_fast64_t> oneMinusPTransitions;

        std::map<uint_fast64_t, RationalFunction> directProbs;

        std::map<storm::RationalFunctionVariable, RationalFunction> pRationalFunctions;
        std::map<storm::RationalFunctionVariable, RationalFunction> oneMinusPRationalFunctions;
        
        for (auto const& entry : flexibleMatrix.getRow(state)) {
            // Identify parameter of successor (or constant)
            RationalFunctionVariable parameterOfSuccessor;
            for (auto const& entry2 : flexibleMatrix.getRow(entry.getColumn())) {
                if (entry2.getValue().isZero()) {
                    continue;
                }
                if (entry2.getValue().isConstant()) {
                    parameterOfSuccessor = constantVariable;
                    break;
                }
                
                STORM_LOG_ERROR_COND(entry2.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
                parameterOfSuccessor = *entry2.getValue().gatherVariables().begin();
                auto cache = entry2.getValue().nominatorAsPolynomial().pCache();
                STORM_LOG_ERROR_COND(entry2.getValue().denominator().isOne() && entry2.getValue().nominator().isUnivariate() &&
                                         entry2.getValue().nominator().getSingleVariable() == parameterOfSuccessor &&
                                         entry2.getValue().nominator().factorization().size() == 1,
                                     "Flip minimization only supports simple pMCs.");
                STORM_LOG_ERROR_COND(flexibleMatrix.getRow(entry.getColumn()).size() == 2,
                                     "Flip minimization only supports transitions with a single parameter.");
                if (utility::isOne(entry2.getValue().derivative(entry2.getValue().nominator().getSingleVariable()))) {
                    pRationalFunctions[parameterOfSuccessor] = entry2.getValue();
                    pTransitions[entry.getColumn()] = entry2.getColumn();
                } else if (utility::isOne(-entry2.getValue().derivative(entry2.getValue().nominator().getSingleVariable()))) {
                    oneMinusPRationalFunctions[parameterOfSuccessor] = entry2.getValue();
                    oneMinusPTransitions[entry.getColumn()] = entry2.getColumn();
                } else {
                    STORM_LOG_ERROR_COND(false, "Flip minimization only supports transitions with a single parameter.");
                }
            }
            parameterBuckets[parameterOfSuccessor].emplace(entry.getColumn());
            cumulativeProbabilities[parameterOfSuccessor] += entry.getValue();
            directProbs[entry.getColumn()] = entry.getValue();
        }
        
        // TODO slow could be done better if flexible matrix had ability to add states
        uint_fast64_t newMatrixSize = flexibleMatrix.getRowCount() + 3 * parameterBuckets.size();
        if (parameterBuckets.count(constantVariable)) {
            newMatrixSize -= 2;
        }
        storage::SparseMatrixBuilder<RationalFunction> builder;
        storage::FlexibleSparseMatrix<RationalFunction> matrixWithAdditionalStates(builder.build(newMatrixSize, newMatrixSize, 0));
        for (uint_fast64_t row = 0; row < flexibleMatrix.getRowCount(); row++) {
            matrixWithAdditionalStates.getRow(row) = flexibleMatrix.getRow(row);
        }

        std::map<RationalFunctionVariable, std::set<uint_fast64_t>> countingWorkingSets;

        uint_fast64_t newStateIndex = flexibleMatrix.getRowCount();
        matrixWithAdditionalStates.getRow(state).clear();
        for (auto const& entry : parameterBuckets) {
            matrixWithAdditionalStates.getRow(state).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(newStateIndex, cumulativeProbabilities.at(entry.first)));
            std::cout << "Reorder: " << state << " -> " << newStateIndex << std::endl;
            
            if (entry.first == constantVariable) {
                for (auto const& successor : entry.second) {
                    matrixWithAdditionalStates.getRow(newStateIndex).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(
                        successor, directProbs.at(successor) / cumulativeProbabilities.at(entry.first)
                    ));
                }
                // Issue: multiple transitions can go to a single state, not allowed
                // Solution: Join them
                matrixWithAdditionalStates.getRow(newStateIndex) = joinDuplicateTransitions(matrixWithAdditionalStates.getRow(newStateIndex));
                newStateIndex += 1;
            } else {
                matrixWithAdditionalStates.getRow(newStateIndex).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(newStateIndex + 1, pRationalFunctions.at(entry.first)));
                matrixWithAdditionalStates.getRow(newStateIndex).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(newStateIndex + 2, oneMinusPRationalFunctions.at(entry.first)));
                
                for (auto const& successor : entry.second) {
                    // Remove transition from being counted (for now, we will re-add it below)
                    for (auto& state : treeStates.at(entry.first)) {
                        if (state.first != successor) {
                            state.second.erase(successor);
                        }
                    }
                    // If it's still needed, re-count it
                    workingSets.at(entry.first).emplace(successor);

                    matrixWithAdditionalStates.getRow(newStateIndex + 1).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(
                        pTransitions.at(successor), directProbs.at(successor) / cumulativeProbabilities.at(entry.first)
                    ));
                    matrixWithAdditionalStates.getRow(newStateIndex + 2).push_back(storage::MatrixEntry<uint_fast64_t, RationalFunction>(
                        oneMinusPTransitions.at(successor), directProbs.at(successor) / cumulativeProbabilities.at(entry.first)
                    ));
                }
                // Issue: multiple transitions can go to a single state, not allowed
                // Solution: Join them
                matrixWithAdditionalStates.getRow(newStateIndex + 1) = joinDuplicateTransitions(matrixWithAdditionalStates.getRow(newStateIndex + 1));
                matrixWithAdditionalStates.getRow(newStateIndex + 2) = joinDuplicateTransitions(matrixWithAdditionalStates.getRow(newStateIndex + 2));

                treeStates.at(entry.first)[newStateIndex].emplace(newStateIndex);
                workingSets.at(entry.first).emplace(newStateIndex);

                newStateIndex += 3;
            }
        }

        std::set<RationalFunctionVariable> paramsToUpdate;
        
        for (uint_fast64_t i = flexibleMatrix.getRowCount(); i < newMatrixSize; i++) {
            topologicalOrderingStack.push(i);
        }
        
        updateTreeStates(treeStates, countingWorkingSets, matrixWithAdditionalStates, paramsToUpdate);
        
        
        // Extend labeling to more states (Found no better way to do it)
        models::sparse::StateLabeling nextNewLabels(newMatrixSize);
        for (auto const& label : newLabels.getLabels()) {
            nextNewLabels.addLabel(label);
        }
        for (uint_fast64_t state = 0; state < flexibleMatrix.getRowCount(); state++) {
            for (auto const& label : newLabels.getLabelsOfState(state)) {
                nextNewLabels.addLabelToState(label, state);
            }
        }
        newLabels = nextNewLabels;
        flexibleMatrix = matrixWithAdditionalStates;
        // models::sparse::Dtmc<RationalFunction> newnewnewDTMC(flexibleMatrix.createSparseMatrix(), newLabels);
        // newnewnewDTMC.writeDotToStream(std::cout);
    }

    auto newTransitionMatrix = flexibleMatrix.createSparseMatrix();
    transitionMatrix = newTransitionMatrix;

    models::sparse::Dtmc<RationalFunction> newDTMC(transitionMatrix, newLabels);

    if (stateRewardVector) {
        models::sparse::StandardRewardModel<RationalFunction> newRewardModel(*stateRewardVector);
        newDTMC.addRewardModel(*stateRewardName, newRewardModel);
    }

    // 
    // storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<RationalFunction>> simplifier(newDTMC);
    // STORM_LOG_ASSERT(simplifier.simplify(checkTask.getFormula()), "Could not simplify model.");
    // newDTMC = *simplifier.getSimplifiedModel()->template as<models::sparse::Dtmc<RationalFunction>>();
    // newDTMC.writeDotToStream(std::cout);
    STORM_LOG_ASSERT(newDTMC.getTransitionMatrix().isProbabilistic(), "Internal error: resulting matrix not probabilistic!");

    return newDTMC;
}

std::vector<storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> EqualParameterReducer::joinDuplicateTransitions(std::vector<storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> const& entries) {
    std::vector<uint_fast64_t> keyOrder;
    std::map<uint_fast64_t, storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> existingEntries;
    for (auto const& entry : entries) {
        if (existingEntries.count(entry.getColumn())) {
            existingEntries.at(entry.getColumn()).setValue(existingEntries.at(entry.getColumn()).getValue() + entry.getValue());
        } else {
            existingEntries[entry.getColumn()] = entry;
            keyOrder.push_back(entry.getColumn());
        }
    }
    std::vector<storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> newEntries;
    for (uint_fast64_t key : keyOrder) {
        newEntries.push_back(existingEntries.at(key));
    }
    return newEntries;
}

void EqualParameterReducer::updateTreeStates(
    std::map<RationalFunctionVariable, std::map<uint_fast64_t, std::set<uint_fast64_t>>>& treeStates,
    std::map<RationalFunctionVariable, std::set<uint_fast64_t>>& workingSets,
    storage::FlexibleSparseMatrix<RationalFunction>& flexibleMatrix,
    const std::set<carl::Variable>& allParameters
) {
    auto backwardsTransitions = flexibleMatrix.createSparseMatrix().transpose(true);
    for (auto const& parameter : allParameters) {
        std::set<uint_fast64_t> workingSet = workingSets.at(parameter);
        while (!workingSet.empty()) {
            std::set<uint_fast64_t> newWorkingSet;
            for (uint_fast64_t row : workingSet) {
                for (auto const& entry : backwardsTransitions.getRow(row)) {
                    if (entry.getValue().isConstant()) {
                        // If the set of tree states at the current position is a subset of the set of
                        // tree states of the parent state, we've reached some loop. Then we can stop.
                        bool isSubset = true;
                        for (auto const& state : treeStates.at(parameter).at(row)) {
                            if (!treeStates.at(parameter)[entry.getColumn()].count(state)) {
                                isSubset = false;
                                break;
                            }
                        }
                        if (isSubset) {
                            continue;
                        }
                        for (auto const& state : treeStates.at(parameter).at(row)) {
                            treeStates.at(parameter).at(entry.getColumn()).emplace(state);
                        }
                        newWorkingSet.emplace(entry.getColumn());
                    }
                }
            }
            workingSet = newWorkingSet;
        }
    }
}



bool EqualParameterReducer::jipConvert(uint_fast64_t state, storage::FlexibleSparseMatrix<RationalFunction>& matrix, std::map<uint_fast64_t, bool>& alreadyVisited,
                                       const std::map<RationalFunctionVariable, std::map<uint_fast64_t, std::set<uint_fast64_t>>>& treeStates,
                                       const std::set<carl::Variable>& allParameters) {
    auto copiedRow = matrix.getRow(state);
    bool firstIteration = true;
    for (auto const& entry : copiedRow) {
        // ignore zero-entries
        if (entry.getValue().isZero()) {
            continue;
        }
        // if this is a parameteric transition, for now this means returning and ending our preprocessing
        if (!entry.getValue().isConstant()) {
            return false;
        }
        uint_fast64_t nextState = entry.getColumn();
        bool constantTransition;
        if (alreadyVisited.count(nextState)) {
            constantTransition = alreadyVisited.at(nextState);
        } else {
            alreadyVisited[nextState] = false;
            constantTransition = jipConvert(nextState, matrix, alreadyVisited, treeStates, allParameters);
            alreadyVisited[nextState] = constantTransition;
        }
        RationalFunction probability = entry.getValue();
        if (firstIteration) {
            matrix.getRow(state).clear();
            firstIteration = false;
        }
        if (constantTransition) {
            for (auto const& successor : matrix.getRow(nextState)) {
                RationalFunction succProbability = successor.getValue();
                storm::storage::MatrixEntry<uint_fast64_t, RationalFunction> newEntry(successor.getColumn(), probability * succProbability);
                std::cout << "JipConvert: " << state << " -> " << successor.getColumn() << " w/ " << probability * succProbability << std::endl;
                matrix.getRow(state).push_back(newEntry);
            }
        } else {
            matrix.getRow(state).push_back(entry);
        }
    }
    std::sort(matrix.getRow(state).begin(), matrix.getRow(state).end(),
              [](const storage::MatrixEntry<uint_fast64_t, RationalFunction>& a, const storage::MatrixEntry<uint_fast64_t, RationalFunction>& b) -> bool {
                  return a.getColumn() > b.getColumn();
              });
    return true;
}

class EqualParameterReducer;
}  // namespace derivative
}  // namespace storm