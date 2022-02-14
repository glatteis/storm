#include "EqualParameterReducer.h"
#include <cstdint>
#include <map>
#include <queue>
#include <string>
#include "adapters/RationalFunctionAdapter.h"
#include "adapters/RationalNumberAdapter.h"
#include "logic/UntilFormula.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "models/sparse/StandardRewardModel.h"
#include "storage/BitVector.h"
#include "storage/SparseMatrix.h"
#include "storm-pars/transformer/SparseParametricDtmcSimplifier.h"
#include "utility/constants.h"
#include "utility/macros.h"

namespace storm {
namespace derivative {

models::sparse::Dtmc<RationalFunction> EqualParameterReducer::minimizeEqualParameters(models::sparse::Dtmc<RationalFunction> dtmc, modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask) {
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
        for (auto const& parameter : allParameters) {
            alreadyVisitedStates[parameter] = std::map<uint_fast64_t, std::vector<uint_fast64_t>>();
            talliedUpProbabilities[parameter] = std::map<uint_fast64_t, std::vector<RationalNumber>>();
            increasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            decreasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            joiningLabels[parameter] = std::map<uint_fast64_t, std::set<std::string>>();
        }
        // Search for all occurences of parameters and fill up the maps
        for (uint_fast64_t row = 0; row < transitionMatrix.getRowCount(); row++) {
            for (auto const& entry : transitionMatrix.getRow(row)) {
                if (entry.getValue().isConstant()) {
                    continue;
                }
                STORM_LOG_ASSERT(entry.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
                auto parameter = *entry.getValue().gatherVariables().begin();
                auto cache = entry.getValue().nominatorAsPolynomial().pCache();
                RationalFunction parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
                if (entry.getValue() == parameterAsFunction) {
                    pFunctions[parameter] = entry.getValue();
                    numberOfSearchingTransitions++;
                    auto parameter = entry.getValue().nominatorAsPolynomial().getSingleVariable();
                    STORM_LOG_ASSERT(alreadyVisitedStates[parameter].count(row) == 0, "Flip minimization only supports simple pMCs.");
                    alreadyVisitedStates[parameter][row] = std::vector<uint_fast64_t>();
                    alreadyVisitedStates[parameter][row].push_back(row);
                    talliedUpProbabilities[parameter][row] = std::vector<RationalNumber>();
                    talliedUpProbabilities[parameter][row].push_back(utility::one<RationalNumber>());
                    increasingSuccessors[parameter][row] = entry.getColumn();
                    doneSearching[parameter][row] = false;
                    joiningLabels[parameter][row] = newLabels.getLabelsOfState(row);
                } else if (utility::one<RationalFunction>() - entry.getValue() == parameterAsFunction) {
                    decreasingSuccessors[parameter][row] = entry.getColumn();
                } else {
                    STORM_LOG_ASSERT(false, "Flip minimization only supports simple pMCs.");
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
                    if (doneSearching[parameter][firstState]) {
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
                        if (joiningLabels[parameter][firstState] != newLabels.getLabelsOfState(entry.getColumn()) ||
                            (stateRewardVector && !(*stateRewardVector)[visitedStates.back()].isZero())) {
                            transitionsDoneSearching++;
                            doneSearching[parameter][firstState] = true;
                        }
                        auto constantProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                        currentState = entry.getColumn();
                        
                        // The next probability is the previous one times the considered tranistion
                        talliedUpProbabilities[parameter][firstState].push_back(talliedUpProbabilities[parameter][firstState].back() * constantProbability);
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
                        if (joiningLabels[parameter][otherFirstState] != joiningLabels[parameter][firstState]) {
                            continue;
                        }
                        auto otherVisitedStates = pair3.second;

                        auto state = otherVisitedStates.back();

                        // doneSearching[parameter][otherFirstState] can only be false if the other path has gotten here,
                        // but has not noticed it is the end of the path yet. The easiest way to mitigate this is let
                        // that path join us when it notices.
                        if (state == currentState && doneSearching[parameter][otherFirstState]) {
                            // Wow!!! Another search of the same parameter has gotten here.
                            if (!joinedFirstStates[parameter].count(state)) {
                                joinedFirstStates[parameter][state] = std::map<uint_fast64_t, RationalNumber>();
                            }
                            RationalNumber p1 = talliedUpProbabilities[parameter][firstState].back();
                            joinedFirstStates[parameter][state][firstState] = p1;
                            RationalNumber p2 = talliedUpProbabilities[parameter][otherFirstState].back();
                            joinedFirstStates[parameter][state][otherFirstState] = p2;
                            
                            break;
                        }
                    }
                }
            }
        }
        
        // Check if we found any parameters that we can join into one transition, and if so, do it
        for (auto const& parameter : allParameters) {
            if (joinedFirstStates[parameter].size() > 0) {
                somethingChanged = true;
                break;
            }
        }
        // States that are on parametric paths
        std::map<storm::RationalFunctionVariable, std::set<uint_fast64_t>> statesOnParametricPaths;
        for (auto const& parameter : allParameters) {
            for (auto const& joinedFirstStatesPair : joinedFirstStates[parameter]) {
                for (auto const& firstStatePair : joinedFirstStatesPair.second) {
                    auto firstState = firstStatePair.first;
                    for (auto const& visitedState : alreadyVisitedStates[parameter][firstState]) {
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
            if (joinedFirstStates[parameter].size() == 0) {
                continue;
            }
            std::cout << "joining " << joinedFirstStates[parameter].size() << " cases for " << parameter << std::endl;
            uint_fast64_t runningCounter = wipMatrix.getRowCount();

            std::map<uint_fast64_t, std::vector<uint_fast64_t>> newStates;
            for (auto const& pair : joinedFirstStates[parameter]) {
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
                if (joinedFirstStates[parameter].count(row)) {
                    std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>, std::greater<std::pair<uint_fast64_t, RationalFunction>>> insertInOrder;
                    // For checking correctness
                    RationalFunction rowSum = utility::zero<RationalFunction>();
                    // We will join until these transitions
                    auto statesBeforeParams = joinedFirstStates[parameter][row];
                    
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
                        auto visitedStates = alreadyVisitedStates[parameter][startingState.first];
                        // The probability to reach the starting state. It will be divided as we step back
                        RationalNumber probability = talliedUpProbabilities[parameter][startingState.first].back();
                        for (uint_fast64_t i = 1; i < visitedStates.size(); i++) {
                            // Previous (next in reality) visited state is i-1. Divide the probability
                            auto const& visitedState = visitedStates[i];
                            if (alreadyConsidered.count(visitedState)) {
                                break;
                            }
                            alreadyConsidered.emplace(visitedState);
                            auto const& nextVisitedState = visitedStates[i-1];
                            // get the transition probability between visitedState and nextVisitedState
                            boost::optional<RationalNumber> transitionProbability;
                            for (auto const& entry : wipMatrix.getRow(visitedState)) {
                                if (entry.getColumn() == nextVisitedState) {
                                    transitionProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                                    break;
                                }
                            }
                            STORM_LOG_ASSERT(transitionProbability, "Very bad internal error: could not trace back found path");
                            probability = probability / *transitionProbability;
                            // Go through the matrix entries to see what other stuff is branching here
                            for (auto const& entry : wipMatrix.getRow(visitedState)) {
                                auto state = entry.getColumn();
                                if (!statesOnParametricPaths[parameter].count(state)) {
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
                    insertInOrder.emplace(std::make_pair(newStates[row][0], utility::convertNumber<RationalFunction>(summedProbability)));
                    rowSum += utility::convertNumber<RationalFunction>(summedProbability);

                    STORM_LOG_ASSERT(utility::isOne(rowSum), "Internal error: Row sum of row " << row << " is " << rowSum << ", not one!");
                    while (!insertInOrder.empty()) {
                        auto pair = insertInOrder.top();
                        insertInOrder.pop();
                        builder.addNextValue(row, pair.first, pair.second);
                    }
                } else { // Else, just add the stuff back into the matrix
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        builder.addNextValue(row, entry.getColumn(), entry.getValue());
                    }
                }
            }
            
            // RationalFunction parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
            RationalFunction parameterAsFunction = pFunctions[parameter];

            for (auto const& pair : joinedFirstStates[parameter]) {
                auto row = pair.first;
                RationalNumber summedProbability = 0;
                for (auto const& entry2 : joinedFirstStates[parameter][row]) {
                    summedProbability += entry2.second;
                }
                // Add parametric transitions 
                builder.addNextValue(newStates[row][0], newStates[row][1], parameterAsFunction);
                builder.addNextValue(newStates[row][0], newStates[row][2], utility::one<RationalFunction>() - parameterAsFunction);
                
                // Add transitions from p and 1-p states to actual successors
                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>, std::greater<std::pair<uint_fast64_t, RationalFunction>>> insertInOrder1;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto increasingSuccessor = increasingSuccessors[parameter][id];
                    insertInOrder1.emplace(std::make_pair(increasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder1.empty()) {
                    auto pair = insertInOrder1.top();
                    insertInOrder1.pop();
                    builder.addNextValue(newStates[row][1], pair.first, pair.second);
                }
                
                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>, std::greater<std::pair<uint_fast64_t, RationalFunction>>> insertInOrder2;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto decreasingSuccessor = decreasingSuccessors[parameter][id];
                    insertInOrder2.emplace(std::make_pair(decreasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder2.empty()) {
                    auto pair = insertInOrder2.top();
                    insertInOrder2.pop();
                    builder.addNextValue(newStates[row][2], pair.first, pair.second);
                }
                
                auto someJoinedState = joinedFirstStates[parameter][row].begin()->first;
                for (auto const& label : joiningLabels[parameter][someJoinedState]) {
                    if (!newLabels.getLabels().count(label)) {
                        newLabels.addLabel(label);
                    }
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        newLabels.addLabelToState(label, newStates[row][i]);
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
            for (auto const& pair : joinedFirstStates[parameter]) {
                auto joiningState = pair.first;
                for (auto const& pair2 : pair.second) {
                    // The identifier - the state the search started from
                    auto firstState = pair2.first;
                    auto visitedStates = alreadyVisitedStates[parameter][firstState];

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
        newLabels = newLabels.getSubLabeling(statesToKeep);
        if (stateRewardVector) {
            std::vector<RationalFunction> newStateRewards;
            for (uint_fast64_t i = 0; i < stateRewardVector->size(); i++) {
                if (statesToKeep[i]) {
                    newStateRewards.push_back((*stateRewardVector)[i]);
                }
            }
            stateRewardVector = newStateRewards;
        }
        transitionMatrix = wipMatrix;
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

class EqualParameterReducer;
}
}